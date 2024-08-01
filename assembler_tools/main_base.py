import os
import shutil
import json
import logging
from typing import List, Type

from assembler_tools.utils import AssemblyFailedException
from assembler_tools.ani_dendrogram import ani_clustermap, group_by_linkage
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembler_tools', 'templates'),
    autoescape=select_autoescape(['html'])
)
template_overview = env.get_template('overview.html')
template_assemblies = env.get_template('assemblies.html')


def process_sample(sample: str, sample_dir: str, importers: List[Type[AssemblyImporter]]):
    assemblies, messages = load_assemblies(sample, sample_dir, importers)
    similarity_matrix, ani_html, linkage_matrix = ani_clustermap(assemblies=assemblies)

    dotplot_tables = create_all_dotplots(assemblies, similarity_matrix, linkage_matrix, sample_dir)

    if similarity_matrix is not None:
        similarity_matrix.to_csv(f'{sample_dir}/assemblies_pyskani_similarity_matrix.tsv', sep='\t')

    with open(f"{sample_dir}/assemblies.json", 'w') as f:
        json.dump({assembly.assembler: assembly.to_json() for assembly in assemblies}, f, indent=2)

    template_assemblies.stream(
        messages=messages,
        sample=sample,
        assemblies=assemblies,
        ani_html=ani_html,
        dotplot_tables=dotplot_tables,
    ).dump(f"{sample_dir}/assemblies.html")

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


from assembler_tools.dotplots import create_dotplots, add_separators, DotplotConfig

config = DotplotConfig(kmer=40, background_colour='white',
                       color_fwd='mediumblue', color_rev='firebrick', line_color='black')


def create_dotplot(cluster_id, contig_group_1, contig_group_2, output_dir: str, res=1000):
    if len(contig_group_2) > 1_000_000:
        print('too large')
        return ''

    # concatenate sequences
    seq_a = ''.join(contig.sequence for contig in contig_group_1.contigs)
    seq_b = ''.join(contig.sequence for contig in contig_group_2.contigs)

    # calculate bp_per_pixel: make sure that the dotplot is not too large
    bp_per_pixel = max(len(seq_a), len(seq_b)) / res
    config.bp_per_pixel = bp_per_pixel

    # add separators
    separators_a = [0] + [len(contig.sequence) for contig in contig_group_1.contigs]
    separators_b = [0] + [len(contig.sequence) for contig in contig_group_2.contigs]

    # create dotplot
    image, draw = create_dotplots(seq_a, seq_b, config)
    add_separators(image, draw, separators_a, separators_b, config)

    filename = f"{cluster_id}-{contig_group_1.id}-{contig_group_2.id}.png"
    output_file = os.path.join(output_dir, filename)
    image.save(output_file)
    return filename


def create_all_dotplots(assemblies, similarity_matrix, linkage_matrix, sample_dir: str):
    dotplot_outdir = os.path.join(sample_dir, 'dotplots')
    os.makedirs(dotplot_outdir, exist_ok=True)
    cgs = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}
    clusters = group_by_linkage(similarity_matrix, linkage_matrix)

    tables = {}
    for cluster, contig_groups in clusters.items():
        # compare all pairs of contig groups in the cluster, only once though
        tables[cluster] = create_cluster_dotplots(cluster, contig_groups, dotplot_outdir, cgs)
        # break

    print(tables)
    return tables


def create_cluster_dotplots(cluster, contig_groups: [str], dotplot_outdir: str, cgs):
    table = {}
    dotplots = {}
    for i, cg_name_i in enumerate(contig_groups):
        table[cg_name_i] = {}
        for j, cg_name_j in enumerate(contig_groups):
            cg_i, cg_j = cgs[cg_name_i], cgs[cg_name_j]

            if len(cg_j) < len(cg_i):  # faster
                cg_i, cg_j = cg_j, cg_i

            if (cg_i, cg_j) not in dotplots:
                dotplots[(cg_i.id, cg_j.id)] = create_dotplot(i, cg_i, cg_j, output_dir=dotplot_outdir)

            table[cg_name_i][cg_name_j] = {
                'dotplot': dotplots[(cg_i.id, cg_j.id)],
                'is_inverted': (cg_name_i, cg_name_j) not in dotplots
            }

    return table


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
