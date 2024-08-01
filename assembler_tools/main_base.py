import os
import shutil
import json
import logging
from typing import List, Type

from assembler_tools.utils import AssemblyFailedException
from assembler_tools.ani_dendrogram import ani_clustermap
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
    distance_matrix, ani_html = ani_clustermap(assemblies=assemblies)

    if distance_matrix is not None:
        distance_matrix.to_csv(f'{sample_dir}/assemblies_pyskani_similarity_matrix.tsv', sep='\t')

    with open(f"{sample_dir}/assemblies.json", 'w') as f:
        json.dump({assembly.assembler: assembly.to_json() for assembly in assemblies}, f, indent=2)

    template_assemblies.stream(
        messages=messages,
        sample=sample,
        assemblies=assemblies,
        ani_html=ani_html,
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
