import os
import shutil
import json
import logging
import tempfile
import time
from textwrap import indent
from typing import List, Type
import multiprocessing as mp

from assembler_tools.ContigGroup import ContigGroup
from assembler_tools.dotplots_minimap2 import process_cluster
from assembler_tools.utils import AssemblyFailedException, rgb_array_to_css, css_escape
from assembler_tools.ani_dendrogram import ani_clustermap, add_cluster_info_to_assemblies
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembler_tools', 'templates'),
    autoescape=select_autoescape(['html'])
)
env.filters['css_escape'] = css_escape
template_overview = env.get_template('overview.html')
template_assemblies = env.get_template('assemblies.html.jinja2')
template_assemblies_css = env.get_template('assemblies_dynamic.css.jinja2')


def process_sample(sample: str, sample_dir: str, importers: List[Type[AssemblyImporter]]):
    try:
        assemblies, messages = load_assemblies(sample, sample_dir, importers)
    except AssemblyFailedException as e:
        logging.warning(str(e))
        template_assemblies.stream(
            messages=[e], sample=sample,
            assemblies=[], sample_to_cluster={}, cluster_to_color={},
        ).dump(f"{sample_dir}/assemblies.html")
        return []

    similarity_matrix, sample_to_cluster = ani_clustermap(
        assemblies=assemblies,
        fname=f"{sample_dir}/ani_clustermap.svg"
    )
    similarity_matrix.to_csv(f'{sample_dir}/assemblies_pyskani_similarity_matrix.tsv', sep='\t')

    add_cluster_info_to_assemblies(assemblies, sample_to_cluster)

    cluster_to_color = create_all_dotplots(assemblies, sample_dir)

    json_data = {assembly.assembler: assembly.to_json() for assembly in assemblies}

    with open(f"{sample_dir}/assemblies.json", 'w') as f:
        json.dump(json_data, f, indent=2)

    for contig in sample_to_cluster:
        sample_to_cluster[contig]['color'] = rgb_array_to_css(sample_to_cluster[contig]['color'])

    template_assemblies.stream(
        messages=messages,
        sample=sample,
        assemblies=assemblies,
        sample_to_cluster=sample_to_cluster,
        cluster_to_color={k: rgb_array_to_css(v) for k, v in cluster_to_color.items()},
    ).dump(f"{sample_dir}/assemblies.html")

    template_assemblies_css.stream(
        assemblies=assemblies
    ).dump(f"{sample_dir}/assemblies_dynamic.css")

    return assemblies


def load_assemblies(sample: str, sample_dir: str, importers: List[Type[AssemblyImporter]]) -> ([Assembly], [str]):
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
    return assemblies, messages


def create_all_dotplots(assemblies, sample_dir: str):
    dotplot_outdir = os.path.join(sample_dir, 'dotplots')
    os.makedirs(dotplot_outdir, exist_ok=True)
    cgs = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}

    cluster_to_color = {cg.cluster_id: cg.cluster_color for cg in cgs.values()}
    cluster_to_color = dict(sorted(cluster_to_color.items(), key=lambda item: item[0]))

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

    # multiprocessing
    with mp.Pool(mp.cpu_count()) as pool:
        pool.starmap(
            func=process_cluster,
            iterable=[(cluster_id, tmpdirs[cluster_id].name, dotplot_outdir) for cluster_id in cluster_to_color])

    # # no multiprocessing
    # for cluster_id in cluster_to_color:
    #     process_cluster(cluster_id, tmpdirs[cluster_id].name, dotplot_outdir)

    # cleanup
    for tmpdir in tmpdirs.values():
        tmpdir.cleanup()

    return cluster_to_color


def prepare_website(samples_dir: str, link: bool = True):
    def link(src, dst):
        if os.path.isfile(dst):
            os.remove(dst)
        os.link(src, dst)

    copy_or_link = link if link else shutil.copy

    # Add dotplot.js
    copy_or_link(
        os.path.join('/home/thomas/PycharmProjects/dotplot.js/dotplot.js'),
        f"{samples_dir}/dotplot.js"
    )
    # use link instead

    # Add assemblies.css
    copy_or_link(
        os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.css'),
        f"{samples_dir}/assemblies.css"
    )

    # Add assemblies.js
    copy_or_link(
        os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.js'),
        f"{samples_dir}/assemblies.js"
    )
