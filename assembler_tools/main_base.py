import os
import shutil
import json
import logging
from typing import List, Type

from assembler_tools.ContigGroup import ContigGroup
from assembler_tools.utils import AssemblyFailedException, rgb_array_to_css, css_escape
from assembler_tools.ani_dendrogram import ani_clustermap, add_cluster_info_to_assemblies
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter
from assembler_tools.dotplots import DotplotConfig, create_dotplots

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
    assemblies, messages = load_assemblies(sample, sample_dir, importers)
    similarity_matrix, sample_to_cluster = ani_clustermap(
        assemblies=assemblies,
        fname=f"{sample_dir}/ani_clustermap.svg"
    )
    add_cluster_info_to_assemblies(assemblies, sample_to_cluster)
    if similarity_matrix is not None:
        similarity_matrix.to_csv(f'{sample_dir}/assemblies_pyskani_similarity_matrix.tsv', sep='\t')

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
        msg = f"Failed to load any assemblies for {sample}"
        logging.warning(msg)
        messages.insert(0, msg)
    else:
        print(f"Loaded {len(assemblies)} assemblies for {sample}")
    return assemblies, messages


def create_all_dotplots(assemblies, sample_dir: str):
    dotplot_outdir = os.path.join(sample_dir, 'dotplots')
    os.makedirs(dotplot_outdir, exist_ok=True)
    cgs = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}

    cluster_to_color = {cg.cluster_id: cg.cluster_color for cg in cgs.values()}
    # sort by data['cluster']
    cluster_to_color = dict(sorted(cluster_to_color.items(), key=lambda item: item[0]))

    for cluster_id, color in cluster_to_color.items():
        cluster_cgs = [cg for cg in cgs.values() if cg.cluster_id == cluster_id]
        length = sum(len(c) for c in cluster_cgs)
        if length > 1_000_000:
            print(f"Skipping cluster {cluster_id} with total length {length}")
            continue
        create_dotplots(cluster_cgs, output=os.path.join(dotplot_outdir, f'{cluster_id}.svg'))

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
