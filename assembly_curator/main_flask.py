import os
import shutil
import json
import logging
from typing import List, Type

import socket

from assembly_curator.phylogenetic_tree import calculate_phylogenetic_tree

socket.setdefaulttimeout(1000)  # seconds
from werkzeug.serving import WSGIRequestHandler

import dill
from flask import Flask, send_from_directory, render_template, redirect, request, jsonify
from pywebio.input import input, FLOAT, SELECT
from pywebio.output import put_text, put_html
from pywebio.platform.flask import webio_view
from pywebio.session import eval_js

from urllib.parse import urlparse, parse_qs

from assembly_curator.ContigGroup import ContigGroup
from assembly_curator.main_base import process_sample, prepare_website
from assembly_curator.utils import load_importers, load_get_custom_html, get_relative_path, detach_process
from assembly_curator.Assembly import Assembly
from assembly_curator.AssemblyImporter import AssemblyImporter
from assembly_curator.contig_curator import *

samples_directory: str = None
importers: List[Type[AssemblyImporter]] = None
get_custom_html = None

process_assembly_huey = None
assemblies: [Assembly] = None

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembly_curator', 'templates'),
    autoescape=select_autoescape(['html'])
)

app = Flask(__name__)
pywebio_html = """
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet"
      integrity="sha384-QWTKZyjpPEjISv5WaRU9OFeRpok6YctnYmDr5pNlyT2bRjXh0JMhjY6hW+ALEwIH" crossorigin="anonymous">
<script>
    document.documentElement.setAttribute('data-bs-theme', (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'))
</script>
<style>
#input-container, .footer {
    background: unset !important;
}
"""


def get_status(path):
    res = {
        'name': path,
        'custom_html': get_custom_html(path)
    }
    try:
        files = os.listdir(os.path.join(samples_directory, path, 'assembly-curator'))
    except FileNotFoundError:
        return res | {'status': 'not started', 'btn_cls': 'secondary', 'icon': 'bi-pause-circle'}
    if 'hybrid.fasta' in files:
        return res | {'status': 'finished', 'btn_cls': 'success', 'icon': 'bi-check-circle'}
    elif 'assemblies.pkl' in files:
        return res | {'status': 'preprocessed', 'btn_cls': 'primary', 'icon': 'bi-play-circle'}
    elif 'failed' in files:
        return res | {'status': 'failed', 'btn_cls': 'danger', 'icon': 'bi-x-circle'}
    else:
        return res | {'status': 'processing', 'btn_cls': 'warning', 'icon': 'bi-hourglass-split'}


def list_directory(samples_directory, rel_path, is_root=False):
    cwd = os.path.join(samples_directory, rel_path)
    folders = []
    files = []
    links = []

    with os.scandir(cwd) as it:
        for entry in it:
            if entry.is_dir():
                folders.append(entry.name)
            elif entry.is_symlink():
                links.append(dict(name=entry.name, url=os.readlink(entry.path)))
            elif entry.is_file():
                files.append(entry.name)
            else:
                logging.warning(f"Unknown type for {entry.name}")

    if is_root:
        return folders, files, links

    context = {
        'title': f'Index of: {rel_path}',
        'cwd': cwd, 'rel_path': rel_path,
        'folders': sorted(folders), 'files': sorted(files), 'links': sorted(links)
    }

    template_index = env.get_template('index.html.jinja2')
    return template_index.render(context)


@app.route('/')
def serve_root():
    print('>>>>>>>>>>>>>>> serving overview')

    folders, files, links = list_directory(samples_directory, '.', is_root=True)

    calculate_tree = True
    if calculate_tree:
        samples = calculate_phylogenetic_tree(importers, samples_directory, force_rerun=False)
        folders = [f for f in folders if f not in samples]
    else:
        samples, folders = folders, []

    samples = [get_status(path) for path in samples]

    template_index = env.get_template('index.html.jinja2')
    return template_index.render(
        title='Overview',
        samples=samples,
        folders=sorted(folders),
        files=sorted(files),
        links=sorted(links),
        relpath=get_relative_path(f'{samples_directory}/index.html.jinja2', samples_directory)
    )


def process_assembly(sample, sample_dir):
    print(f'>>>>>>>>>>>>>>> processing assembly {sample}')
    assemblies = process_sample(sample, sample_dir, importers)

    if not assemblies:
        # touch /assembly-curator/fail
        with open(f"{sample_dir}/assembly-curator/failed", 'w') as f:
            f.write('Assembly failed.')
        return

    dill.dump(assemblies, open(f"{sample_dir}/assembly-curator/assemblies.pkl", 'wb'))


def serve_assembly(samples_directory, sample):
    sample_dir = os.path.join(samples_directory, sample)

    if not os.path.isfile(f"{sample_dir}/assembly-curator/assemblies.pkl"):
        process_assembly(sample, sample_dir)
    else:
        print(f'>>>>>>>>>>>>>>> serving assembly {sample}')

    return send_from_directory(
        directory=os.path.abspath(sample_dir),
        path='assemblies.html'
    )


@app.route('/reset_sample', methods=['GET', 'POST'])
def reset_sample():
    sample_name = request.json.get('sample_name')
    path = os.path.join(samples_directory, sample_name, 'assembly-curator')
    print(f'>>>>>>>>>>>>>>> Resetting {sample_name} ({path})')

    if os.path.isdir(path):
        shutil.rmtree(path)

    return jsonify({'status': 'success'}), 200


@app.route('/reset_all_samples', methods=['GET', 'POST'])
def reset_all_samples():
    print('>>>>>>>>>>>>>>> Resetting all samples')
    samples, _, _ = list_directory(samples_directory, '.', is_root=True)
    samples = [get_status(path) for path in samples]

    for sample in samples:
        if sample['status'] != 'finished':
            path = os.path.join(samples_directory, sample['name'], 'assembly-curator')
            if os.path.isdir(path):
                shutil.rmtree(path)

    return jsonify({'status': 'success'}), 200


@app.route('/dispatch_all_not_started_samples', methods=['POST'])
def dispatch_all_not_started_samples():
    print('>>>>>>>>>>>>>>> Dispatching all not started samples')
    samples, _, _ = list_directory(samples_directory, '.', is_root=True)
    samples = [get_status(path) for path in samples]

    for sample in samples:
        if sample['status'] == 'not started':
            process_assembly_huey(sample['name'], os.path.join(samples_directory, sample['name']))

    return jsonify({'status': 'success'}), 200


@app.route('/<path:filepath>')
def serve_path(filepath):
    print(f">>>>>>>>>>>>>>> Requested {filepath}")
    full_path = os.path.abspath(os.path.join(samples_directory, filepath))
    dirname, basename = os.path.dirname(full_path), os.path.basename(full_path)
    action = request.args.get('action', 'auto')

    # Import and serve assembly
    if basename == 'assemblies.html':
        sample = os.path.basename(dirname)
        return serve_assembly(samples_directory, sample)

    # Serve files
    if os.path.isfile(full_path):
        extension = full_path.rsplit('.', 1)[-1]
        mimetype, as_attachment = None, False

        if extension in ['fasta', 'tsv', 'gfa', 'gv']:
            mimetype = 'text/plain'

        if action == 'view':
            mimetype = 'text/plain'
        elif action == 'download':
            as_attachment = True

        print(full_path, as_attachment, mimetype)

        # Special case for favicon
        if filepath == 'favicon.ico':
            # return assembly_curator/templates/assembly-curator.webp
            return send_from_directory('templates', 'assembly-curator.webp')

        return send_from_directory(dirname, basename, as_attachment=as_attachment, mimetype=mimetype)

    # Serve directories
    elif os.path.isdir(full_path):
        if not filepath.endswith('/'):
            return redirect(f'/{filepath}/')
        return list_directory(samples_directory, filepath)

    else:
        return "404 Not Found :(", 404


# PyWebIO application function
def pywebio_export():
    put_html(pywebio_html)

    url = eval_js('window.location.href')
    print('URL:', url)
    # get sample and contig_groups from URL
    params = parse_qs(urlparse(url).query)
    sample = params.get('sample', [None])[0]
    contig_groups = params.get('contig_groups', [None])[0]

    if sample is None or contig_groups is None:
        put_text(f"Sample: {sample}")
        put_text(f"Contig groups: {contig_groups}")
        put_text("Aborting. Please provide a sample and contig groups.")
        return

    contig_groups = contig_groups.split(',')
    put_text(f">>>>>>>>>>>>>>> Curating {sample} with {len(contig_groups)} contig groups")

    sample_dir = os.path.join(samples_directory, sample)
    assert os.path.isfile(f"{sample_dir}/assembly-curator/assemblies.pkl")
    assemblies = dill.load(open(f"{sample_dir}/assembly-curator/assemblies.pkl", 'rb'))
    put_text(f"Loaded {len(assemblies)} assemblies for {sample}")
    cgs = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}

    missing = [cg for cg in contig_groups if cg not in cgs]
    if missing:
        put_text(f"Aborting. Missing contig groups: {missing}")
        return

    contig_groups: [ContigGroup] = [cgs[cg] for cg in contig_groups]

    # sort by length
    contig_groups.sort(key=lambda cg: len(cg), reverse=True)
    headers = create_headers(sample, contig_groups)
    export(f"{sample_dir}/assembly-curator/hybrid.fasta", headers)


# Add the PyWebIO endpoint
app.add_url_rule('/export', 'webio_view', webio_view(pywebio_export), methods=['GET', 'POST', 'OPTIONS'])


def run_server(
        samples_dir: str = '/data',
        plugin_dir: str = '/plugins',
        port=8080,
        address='localhost',
        debug: bool = False,
        n_workers: int = 0
):
    assert os.path.isdir(samples_dir), f"Samples directory {samples_dir} does not exist"
    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    global samples_directory, process_assembly_huey, huey
    samples_directory = samples_dir

    os.environ['HUEY_DB_PATH'] = os.path.join(samples_dir, 'huey.db')
    from assembly_curator.huey_tasks import process_assembly as process_assembly_huey, huey

    if n_workers > 0:
        from assembly_curator.huey_main import run_huey
        detach_process(run_huey, samples_dir=samples_dir, plugin_dir=plugin_dir, n_workers=n_workers)

    os.environ['MULTIPROCESSING_DOTPLOTS'] = 'TRUE'

    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    global importers, get_custom_html
    importers = load_importers(plugin_dir)
    get_custom_html = load_get_custom_html(plugin_dir)

    prepare_website(samples_dir, link=False)

    app.run(host=address, port=port, debug=debug)


def main():
    from fire import Fire

    Fire(run_server)


if __name__ == '__main__':
    main()
