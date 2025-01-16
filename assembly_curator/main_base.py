import os
import sys
import shutil
import json
import logging
import tempfile
from typing import List, Type
import multiprocessing as mp
import importlib.resources as pkg_resources

from assembly_curator.ContigGroup import ContigGroup
from assembly_curator.dotplots_minimap2 import process_cluster
from assembly_curator.utils import AssemblyFailedException, rgb_array_to_css, css_escape
from assembly_curator.ani_dendrogram import ani_clustermap, add_cluster_info_to_assemblies
from assembly_curator.Assembly import Assembly
from assembly_curator.AssemblyImporter import AssemblyImporter

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembly_curator', 'templates'),
    autoescape=select_autoescape(['html']),
    cache_size=0 if '--debug' in sys.argv else -1
)
env.filters['css_escape'] = css_escape
template_overview = env.get_template('index.html.jinja2')
template_assemblies = env.get_template('assemblies.html.jinja2')
template_assemblies_css = env.get_template('assemblies_dynamic.css.jinja2')

# Load default GC-content thresholds (os.environ)
GC_LOW = float(os.environ.get('GC_LOW', '25')) / 100
GC_HIGH = float(os.environ.get('GC_HIGH', '65')) / 100


def process_sample(
        sample: str,
        sample_dir: str,
        importers: List[Type[AssemblyImporter]],
        raise_error: bool = False,
        force_rerun: bool = False
) -> [Assembly]:
    template_assemblies = env.get_template('assemblies.html.jinja2')

    if force_rerun:
        shutil.rmtree(f"{sample_dir}/assembly-curator")
    elif os.path.exists(f"{sample_dir}/assembly-curator"):
        msg = f"Sample {sample} is already being processed!"
        if raise_error:
            raise AssertionError(msg)
        else:
            logging.info(msg)
            return []

    os.makedirs(f"{sample_dir}/assembly-curator")

    try:
        assemblies, messages = load_assemblies(sample, sample_dir, importers)
    except AssemblyFailedException as e:
        logging.warning(str(e))
        template_assemblies.stream(
            messages=[e], sample=sample,
            assemblies=[], cluster_to_color={},
        ).dump(f"{sample_dir}/assemblies.html")
        return []

    similarity_matrix, cluster_to_color, cg_to_cluster = ani_clustermap(
        assemblies=assemblies,
        fname=f"{sample_dir}/assembly-curator/ani_clustermap.svg"
    )
    similarity_matrix.to_csv(f'{sample_dir}/assembly-curator/assemblies_pyskani_similarity_matrix.tsv', sep='\t')

    add_cluster_info_to_assemblies(assemblies, cluster_to_color, cg_to_cluster)

    # cluster_to_color =
    create_all_dotplots(assemblies, sample_dir)

    json_data = {assembly.assembler: assembly.to_json() for assembly in assemblies}

    with open(f"{sample_dir}/assembly-curator/assemblies.json", 'w') as f:
        json.dump(json_data, f, indent=2)

    template_assemblies.stream(
        messages=messages,
        sample=sample,
        assemblies=assemblies,
        cluster_to_color={cluster_id: rgb_array_to_css(color) for cluster_id, color in cluster_to_color.items()},
        cg_to_cluster=cg_to_cluster,
    ).dump(f"{sample_dir}/assemblies.html")

    template_assemblies_css.stream(
        assemblies=assemblies
    ).dump(f"{sample_dir}/assembly-curator/assemblies_dynamic.css")

    return assemblies


def load_assemblies(
        sample: str,
        sample_dir: str,
        importers: List[Type[AssemblyImporter]],
        only_one: bool = False
) -> ([Assembly], [str]):
    print(f"Processing {sample}")

    assemblies: [Assembly] = []
    messages: [str] = []
    for importer_class in importers:
        try:
            importer = importer_class(sample_dir)
            print(f"  -> loading {sample} with {importer.name}")
            assembly = importer.load_assembly()
            assemblies.append(assembly)
            assembly.pprint()
            if only_one:
                return [assembly]
        except AssemblyFailedException as e:
            logging.warning(str(e))
            messages.append(e)

    if not assemblies:
        err = f"Failed to load any assemblies for {sample}!"
        if messages:
            err += f"<br><br>{'<br>'.join(str(m) for m in messages)}"
        raise AssemblyFailedException(err)
    else:
        print(f"Loaded {len(assemblies)} assemblies for {sample}")

    # Add warning if in GC-content is below GC_LOW or above GC_HIGH
    for assembly in assemblies:
        for cg in assembly.contig_groups:
            for contig in cg.contigs:
                if contig.gc_rel < GC_LOW:
                    messages.append(AssemblyFailedException(
                        f"GC content below {GC_LOW * 100:.2f} ({contig.gc_rel * 100:.2f}%) for {contig.id}"))
                elif contig.gc_rel > GC_HIGH:
                    messages.append(AssemblyFailedException(
                        f"High GC content above {GC_HIGH * 100:.2f} ({contig.gc_rel * 100:.2f}%) for {contig.id}"))

    return assemblies, messages


def create_all_dotplots(assemblies, sample_dir: str):
    dotplot_outdir = os.path.join(sample_dir, 'assembly-curator', 'dotplots')
    os.makedirs(dotplot_outdir, exist_ok=True)
    cgs = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}

    cluster_to_color = {cg.cluster_id: cg.cluster_color for cg in cgs.values()}

    tmpdirs = {}
    cluster_to_cgs = {}
    for cluster_id in cluster_to_color:
        cluster_cgs = [cg for cg in cgs.values() if cg.cluster_id == cluster_id]
        cluster_to_cgs[cluster_id] = cluster_cgs
        tmpdir = tempfile.TemporaryDirectory()
        for i, cg in enumerate(cluster_cgs):
            with open(os.path.join(tmpdir.name, f'cg_{i}.json'), 'w') as f:
                json.dump(cg.to_json(), f, indent=4)
            with open(os.path.join(tmpdir.name, f'cg_{i}.fasta'), 'w') as f:
                f.write(f'>{cg.id}\n')
                for c in cg.contigs:
                    f.write(c.sequence)
        tmpdirs[cluster_id] = tmpdir

        # start_mm2 = time.time()
        # create_dotplots_minimap2(cluster_cgs, output=os.path.join(dotplot_outdir, f'{cluster_id}.svg'))
        # delta_mm2 = time.time() - start_mm2

    MULTIPROCESSING = os.environ.get('MULTIPROCESSING_DOTPLOTS', 'FALSE').lower() == 'true'
    if MULTIPROCESSING:
        with mp.Pool(mp.cpu_count()) as pool:
            pool.starmap(
                func=process_cluster,
                iterable=[(cluster_id, tmpdirs[cluster_id].name, dotplot_outdir) for cluster_id in cluster_to_color])
    else:
        for cluster_id in cluster_to_color:
            process_cluster(cluster_id, tmpdirs[cluster_id].name, dotplot_outdir)

    # cleanup
    for tmpdir in tmpdirs.values():
        tmpdir.cleanup()

    return cluster_to_color


def prepare_website(samples_dir: str, link: bool = True):
    def copy_or_link(src, dst):
        if os.path.isfile(dst):
            os.remove(dst)
        if link:
            os.link(src, dst)
        else:
            shutil.copy(src, dst)

    files_to_copy = ['dotplot.js', 'assemblies.css', 'assemblies.js']
    for file_name in files_to_copy:
        with pkg_resources.path('assembly_curator.templates', file_name) as src:
            dst = os.path.join(samples_dir, file_name)
            copy_or_link(src, dst)
