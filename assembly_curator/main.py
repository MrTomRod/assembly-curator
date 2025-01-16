import json
import os
import shutil
import logging
from typing import List, Type
from tempfile import TemporaryDirectory

from assembly_curator.main_base import prepare_website, process_sample, load_assemblies
from assembly_curator.phylogenetic_tree import calculate_phylogenetic_tree
from assembly_curator.utils import load_importers, get_relative_path
from assembly_curator.Assembly import Assembly
from assembly_curator.AssemblyImporter import AssemblyImporter

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembly_curator', 'templates'),
    autoescape=select_autoescape(['html'])
)
template_index = env.get_template('index.html.jinja2')

s_to_n = lambda sample: {'name': sample, 'status': 'not started', 'btn_cls': 'secondary', 'icon': 'bi-pause-circle'}


def process_samples(
        importers: [Type[AssemblyImporter]],
        samples_dir: str,
        overview_file: str = None,
        calculate_tree: bool = True,
        force_rerun: bool = False
):
    if overview_file is None:
        overview_file = f'{samples_dir}/index.html'
    folders = [s for s in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, s))]
    print(f"Processing {len(folders)} samples in {samples_dir}")

    if calculate_tree:
        samples = calculate_phylogenetic_tree(importers, samples_dir, force_rerun=force_rerun)
        folders = [f for f in folders if f not in samples]
    else:
        samples, folders = folders, []

    template_index.stream(
        title=f'Index of: {samples_dir}',
        relpath=get_relative_path(overview_file, samples_dir),
        samples=[s_to_n(s) for s in samples],
        folders=sorted(folders),
        files=sorted([]),
        links=sorted([]),
    ).dump(overview_file)

    prepare_website(samples_dir)

    for sample in samples:
        sample_dir = os.path.join(samples_dir, sample)

        assemblies = process_sample(sample, sample_dir, importers, force_rerun=force_rerun)

        # final_assembly = export_assemblies(assemblies)
        # print(f"Final assembly for {sample}: {final_assembly}")


# DATADIR = ''
DATADIR = '-nanopore'


def cli(samples_dir: str = f'./data{DATADIR}', plugin_dir: str = f'./plugins{DATADIR}', force_rerun: bool = False):
    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    if plugin_dir is None:
        plugin_dir = os.environ.get('PLUGIN_DIR', './plugins')
    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    importers = load_importers(plugin_dir)

    process_samples(importers, samples_dir, force_rerun=force_rerun)


def main():
    from fire import Fire

    Fire(cli)


if __name__ == '__main__':
    main()
