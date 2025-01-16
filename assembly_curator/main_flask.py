import os
import sys
import shutil
import logging
from glob import glob
from typing import List, Type

import socket

from assembly_curator.phylogenetic_tree import calculate_phylogenetic_tree

socket.setdefaulttimeout(1000)  # seconds

import dill
from flask import Flask, send_from_directory, redirect, request, jsonify, abort
from pywebio.output import put_text, put_html
from pywebio.input import select, SELECT
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
    autoescape=select_autoescape(['html']),
    cache_size=0 if '--debug' in sys.argv else -1
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
    try:
        with open(os.path.join(samples_directory, path, 'note.md')) as f:
            note = f.read()
    except FileNotFoundError:
        note = None

    res = {
        'name': path,
        'custom_html': get_custom_html(path),
        'note': note
    }

    try:
        files = os.listdir(os.path.join(samples_directory, path, 'assembly-curator'))
    except FileNotFoundError:
        return res | {'status': 'not started', 'btn_cls': 'secondary', 'icon': 'bi-pause-circle'}
    if 'hybrid.fasta' in files:
        return res | {'status': 'finished', 'btn_cls': 'success', 'icon': 'bi-check-circle'}
    elif 'failed' in files:
        return res | {'status': 'failed', 'btn_cls': 'danger', 'icon': 'bi-x-circle'}
    elif 'assemblies.pkl' in files:
        return res | {'status': 'preprocessed', 'btn_cls': 'primary', 'icon': 'bi-play-circle'}
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
    dirs, files, links = list_directory(samples_directory, '.', is_root=True)

    calculate_tree = True
    if calculate_tree:
        samples = calculate_phylogenetic_tree(importers, samples_directory, force_rerun=False)
        folders = [f for f in dirs if f not in samples]
    else:
        samples, folders = dirs, []

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

    return send_from_directory(
        directory=os.path.abspath(sample_dir),
        path='assemblies.html'
    )


@app.route('/reset_sample', methods=['GET', 'POST'])
def reset_sample():
    sample_name = request.json.get('sample_name')
    path = os.path.join(samples_directory, sample_name, 'assembly-curator')

    if os.path.isdir(path):
        shutil.rmtree(path)

    return jsonify({'status': 'success'}), 200


@app.route('/reset_all_samples', methods=['GET', 'POST'])
def reset_all_samples():
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
    samples, _, _ = list_directory(samples_directory, '.', is_root=True)
    samples = [get_status(path) for path in samples]

    for sample in samples:
        if sample['status'] == 'not started':
            process_assembly_huey(sample['name'], os.path.join(samples_directory, sample['name']))

    return jsonify({'status': 'success'}), 200


@app.route('/get_status', methods=['POST'])
def get_status_endpoint():
    sample_name = request.json.get('sample_name')
    if not sample_name:
        return jsonify({"error": "sample_name parameter is required"}), 400
    # Call the get_status function and get the result
    status = get_status(sample_name)
    return jsonify(status)


@app.route('/toggle_failed', methods=['POST'])
def toggle_failed():
    sample_name = request.json.get('sample_name')
    if not sample_name:
        return jsonify({"error": "sample_name parameter is required"}), 400
    sample_dir = os.path.join(samples_directory, sample_name)
    failed_file = f"{sample_dir}/assembly-curator/failed"
    hybrid_file = f"{sample_dir}/assembly-curator/hybrid.fasta"
    if os.path.isfile(failed_file):
        os.remove(failed_file)
    else:
        try:
            os.remove(hybrid_file)
        except FileNotFoundError:
            pass
        with open(failed_file, 'w') as f:
            f.write('Assembly failed.')
    return jsonify({'status': 'success'}), 200


@app.route('/<path:filepath>', methods=['PUT'])
def put_path(filepath):
    full_path = os.path.abspath(os.path.join(samples_directory, filepath))
    dirname, basename = os.path.dirname(full_path), os.path.basename(full_path)
    if not os.path.exists(dirname):
        logging.warning(f"Directory {dirname} does not exist, cannot write {full_path}")
        return abort(404)

    data = request.data.decode()
    if data:
        with open(full_path, 'w') as f:
            f.write(request.data.decode())
    else:
        os.remove(full_path)

    return jsonify({
        'status': 'success',
        'message': f'Wrote {full_path}' if data else f'Deleted {full_path}'}
    ), 200


@app.route('/<path:filepath>', methods=['GET'])
def serve_path(filepath):
    full_path = os.path.abspath(os.path.join(samples_directory, filepath))
    dirname, basename = os.path.dirname(full_path), os.path.basename(full_path)
    action = request.args.get('action', 'auto')

    # Import and serve assembly
    if basename == 'assemblies.html':
        sample = os.path.basename(dirname)
        return serve_assembly(samples_directory, sample)

    # If the path does not exist, try to use glob to get a path
    if not os.path.exists(full_path):
        _globresult = glob(full_path)
        if len(_globresult) == 1:
            full_path = _globresult[0]
            dirname, basename = os.path.dirname(full_path), os.path.basename(full_path)
        else:
            logging.warning(f"Path {full_path} does not exist")
            return abort(404)

    # Serve files
    if os.path.isfile(full_path):
        # Special case for favicon
        if filepath == 'favicon.ico':
            # return assembly_curator/templates/assembly-curator.webp
            return send_from_directory('templates', 'assembly-curator.webp')

        extension = full_path.rsplit('.', 1)[-1]
        mimetype, as_attachment = None, False

        if extension in ['fasta', 'tsv', 'gfa', 'gv']:
            mimetype = 'text/plain'

        if action == 'view':
            mimetype = 'text/plain'
        elif action == 'download':
            as_attachment = True

        return send_from_directory(dirname, basename, as_attachment=as_attachment, mimetype=mimetype)

    # Serve directories
    elif os.path.isdir(full_path):
        if not filepath.endswith('/'):
            return redirect(f'/{filepath}/')
        return list_directory(samples_directory, filepath)

    else:
        logging.warning(f"Path {full_path} is not a file or directory")
        return abort(404)


# PyWebIO application function
def pywebio_export():
    put_html(pywebio_html)

    url = eval_js('window.location.href')
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
    put_text(f"Curating {sample} with {len(contig_groups)} contig groups...")

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

    # show warning if multiple chromosome contigs
    n_chromosomes = [cg.location for cg in contig_groups].count('chromosome')

    if n_chromosomes != 1:
        if select(f'Warning: {n_chromosomes} contig groups were marked as chromosome. Continue?',
                  type=SELECT, options=['Yes', 'No']) == 'No':
            put_text("Aborting.")
            return

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
