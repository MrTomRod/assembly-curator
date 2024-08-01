import dill
from bottle import route, run, static_file, template, redirect, request, response

import os
import json
import logging
from typing import List, Type

from assembler_tools.main_base import load_assemblies, prepare_website, process_sample
from assembler_tools.ani_dendrogram import ani_clustermap
from assembler_tools.utils import load_plugins, get_relative_path, AssemblyFailedException
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter

samples_directory: str = None
importers: List[Type[AssemblyImporter]] = None
assenblies: [Assembly] = None

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembler_tools', 'templates'),
    autoescape=select_autoescape(['html'])
)


def list_directory(samples_directory, rel_path, is_root=False):
    cwd = os.path.join(samples_directory, rel_path)
    folders = []
    files = []
    links = []
    for path in os.listdir(cwd):
        abs_path = os.path.join(cwd, path)
        assert os.path.exists(abs_path), f"Path {abs_path} does not exist"
        if os.path.isdir(abs_path):
            folders.append(path)
        elif os.path.isfile(abs_path):
            files.append(path)
        elif os.path.islink(abs_path):
            links.append(dict(name=path, url=os.readlink(abs_path)))
        else:
            logging.warning(f"Unknown type for {path}")

    if is_root:
        return folders, files, links

    context = {
        'title': f'Index of: {rel_path}',
        'cwd': cwd, 'rel_path': rel_path,
        'folders': sorted(folders), 'files': sorted(files), 'links': sorted(links)
    }

    template_index = env.get_template('index.html')
    return template_index.render(context)


@route('/')
def serve_root():
    # create a list for each sample:
    # [
    #   { name: 'sample1', complete: True },
    #   ...
    # ]
    # name: dirname, complete: whether folder dirname/curated_assembly exists

    def is_complete(dirname):
        return os.path.isdir(os.path.join(samples_directory, dirname, 'curated_assembly'))

    samples, files, links = list_directory(samples_directory, '.', is_root=True)
    samples = [{'name': path, 'complete': is_complete(path)} for path in samples]

    template_index = env.get_template('index.html')
    return template_index.render(
        title='Overview',
        samples=samples,
        files=files,
        relpath=get_relative_path(f'{samples_directory}/overview.html', samples_directory)
    )


def serve_assembly(samples_directory, filepath):
    full_path = os.path.join(samples_directory, filepath)
    sample_dir = full_path.removesuffix('/assemblies.html')
    sample = os.path.basename(sample_dir)
    print('serving assembly', sample_dir)

    assemblies = process_sample(sample, sample_dir, importers, samples_directory)

    dill.dump(assemblies, open(f"{sample_dir}/assemblies.pkl", 'wb'))

    return static_file(filepath, root=samples_directory)


@route('/<filepath:path>')
def serve_path(filepath):
    full_path = os.path.join(samples_directory, filepath)

    # If it is a sample directory (samples_directory/<sample>/assemblies.html, serve_assembly
    if full_path.endswith('/assemblies.html'):
        return serve_assembly(samples_directory, filepath)

    if os.path.isfile(full_path):
        return static_file(filepath, mimetype='auto', download=False, root=samples_directory)
    elif os.path.isdir(full_path):
        # If full_path does not end with /, redirect to the same path with /
        if not filepath.endswith('/'):
            return redirect(f'/{filepath}/')
        # Otherwise, list the directory
        return list_directory(samples_directory, filepath)
    else:
        return "404 Not Found :("


@route('/curate', method='POST')
def curate():
    data = request.json
    if data is None:
        response.status = 400
        return "Invalid JSON"

    sample = data['sample']
    contigs = data['contigs']
    sample_dir = os.path.join(samples_directory, sample)

    assert os.path.isfile(f"{sample_dir}/assemblies.pkl")
    assemblies = dill.load(open(f"{sample_dir}/assemblies.pkl", 'rb'))
    print(f"Loaded {len(assemblies)} assemblies for {sample}")


def run_server(samples_dir: str = 'data-pb', plugin_dir: str = None, port=8080, address='localhost'):
    global samples_directory
    samples_directory = samples_dir

    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    if plugin_dir is None:
        plugin_dir = os.environ.get('PLUGIN_DIR', './plugins-pb')
    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    global importers
    importers = load_plugins(plugin_dir)

    run(host=address, port=port)


def main():
    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    from fire import Fire
    Fire(run_server)


if __name__ == '__main__':
    main()

exit(0)


def run_server(samples_dir: str = 'data-pb', plugin_dir: str = 'plugins-pb', port=8010, address='localhost'):
    @route('/hello/<name>')
    def index(name):
        return template('<b>Hello {{name}}</b>!', name=name)

    run(host='localhost', port=8080)


def curate_assemblies(assemblies: [Assembly]) -> Assembly:
    print(f"Curating {len(assemblies)} assemblies")
    return assemblies[0]


def main():
    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    from fire import Fire

    Fire(run_server)


if __name__ == '__main__':
    main()


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


def process_samples(importers: [Type[AssemblyImporter]], samples_dir: str, overview_file: str = None):
    if overview_file is None:
        overview_file = f'{samples_dir}/overview.html'
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
            output_file=f'{samples_dir}/pyskani_similarity_matrix.tsv'
        )

        with open(f"{sample_dir}/assemblies.json", 'w') as f:
            json.dump({assembly.assembler: assembly.to_json() for assembly in assemblies}, f, indent=2)

        template_assemblies = env.get_template('assemblies.html')
        template_assemblies.stream(
            messages=messages,
            sample=sample,
            assemblies=assemblies,
            ani_html=ani_html,
        ).dump(f"{sample_dir}/assemblies.html")

        def link(src, dst):
            if os.path.isfile(dst):
                os.remove(dst)
            os.link(src, dst)

        copy_or_link = link  # shutil.copy
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

        # final_assembly = curate_assemblies(assemblies)
        # print(f"Final assembly for {sample}: {final_assembly}")
