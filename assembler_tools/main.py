import json
import os
import shutil
import logging
from typing import List, Type

from assembler_tools.main_base import prepare_website, process_sample
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


def process_samples(importers: [Type[AssemblyImporter]], samples_dir: str, overview_file: str = None):
    if overview_file is None:
        overview_file = f'{samples_dir}/overview.html'
    samples = [s for s in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, s))]
    print(f"Processing {len(samples)} samples in {samples_dir}")

    template_overview.stream(
        samples=samples,
        relpath=get_relative_path(overview_file, samples_dir)
    ).dump(overview_file)

    prepare_website(samples_dir)

    import time
    all_time = time.time()

    for sample in samples:
        sample_dir = os.path.join(samples_dir, sample)
        if not os.path.isdir(sample_dir):
            logging.info(f"Skipping {sample} as it is not a directory")
            continue

        process_sample(sample, sample_dir, importers)
        print('XXXXXXXXXXXXXXXX')
        print(f"{time.time() - all_time:.2f}sec.png")
        exit(0)

        # final_assembly = curate_assemblies(assemblies)
        # print(f"Final assembly for {sample}: {final_assembly}")


def cli(samples_dir: str = './data', plugin_dir: str = None):
    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

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
