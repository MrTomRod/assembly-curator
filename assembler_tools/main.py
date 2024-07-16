import json
import os
import shutil
import logging
from typing import List, Type

from assembler_tools.ani_dendrogram import ani_clustermap
from assembler_tools.utils import load_plugins, get_relative_path, AssemblyFailedException
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembler_tools', 'templates'),
    autoescape=select_autoescape(['html'])
)
template_overview = env.get_template('overview.html')
template_assemblies = env.get_template('assemblies.html')


def curate_assemblies(assemblies: [Assembly]) -> Assembly:
    print(f"Curating {len(assemblies)} assemblies")
    return assemblies[0]


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


def process_samples(importers: [Type[AssemblyImporter]], samples_dir: str, overview_file: str = 'overview.html'):
    samples = [s for s in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, s))]
    print(f"Processing {len(samples)} samples in {samples_dir}")

    template_overview.stream(
        samples=samples,
        relpath=get_relative_path(overview_file, samples_dir)
    ).dump(overview_file)

    for sample in samples:
        sample_dir = os.path.join(samples_dir, sample)
        if not os.path.isdir(sample_dir):
            logging.info(f"Skipping {sample} as it is not a directory")
            continue
        assemblies, messages = load_assemblies(sample, sample_dir, importers)
        ani_html = ani_clustermap(
            assemblies=assemblies,
            output_file=f'{samples_dir}/gendiscal_distance_matrix.tsv'
        )

        with open(f"{sample_dir}/assemblies.json", 'w') as f:
            json.dump({assembly.assembler: assembly.to_json() for assembly in assemblies}, f, indent=2)

        template_assemblies.stream(
            messages=messages,
            sample=sample,
            assemblies=assemblies,
            ani_html=ani_html,
        ).dump(f"{sample_dir}/assemblies.html")

        # Add dotplot.js
        shutil.copy(
            os.path.join('/home/thomas/PycharmProjects/dotplot.js/dotplot.js'),
            f"{samples_dir}/dotplot.js"
        )

        # Add assemblies.css
        shutil.copy(
            os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.css'),
            f"{samples_dir}/assemblies.css"
        )

        # Add assemblies.js
        shutil.copy(
            os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.js'),
            f"{samples_dir}/assemblies.js"
        )

        # final_assembly = curate_assemblies(assemblies)
        # print(f"Final assembly for {sample}: {final_assembly}")


def cli(samples_dir: str = './data', plugin_dir: str = None):
    if plugin_dir is None:
        plugin_dir = os.environ.get('PLUGIN_DIR', './plugins')
    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    importers = load_plugins(plugin_dir)

    process_samples(importers, samples_dir)


def main():
    from fire import Fire

    Fire(cli)


if __name__ == '__main__':
    main()
