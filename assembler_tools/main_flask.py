import os
import json
import logging
from typing import List, Type

import socket

socket.setdefaulttimeout(1000)  # seconds
from werkzeug.serving import WSGIRequestHandler


import dill
from flask import Flask, send_from_directory, render_template, redirect, request, jsonify
from pywebio.input import input, FLOAT, SELECT
from pywebio.output import put_text, put_html
from pywebio.platform.flask import webio_view
from pywebio.session import eval_js

from urllib.parse import urlparse, parse_qs

from assembler_tools.ContigGroup import ContigGroup
from assembler_tools.main_base import process_sample
from assembler_tools.utils import load_plugins, get_relative_path, AssemblyFailedException
from assembler_tools.Assembly import Assembly
from assembler_tools.AssemblyImporter import AssemblyImporter
from assembler_tools.contig_curator import *

samples_directory: str = None
importers: List[Type[AssemblyImporter]] = None
assemblies: [Assembly] = None

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader('assembler_tools', 'templates'),
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


@app.route('/')
def serve_root():
    print('>>>>>>>>>>>>>>> serving root')

    def is_complete(dirname):
        return os.path.isfile(os.path.join(samples_directory, dirname, 'curated.fasta'))

    samples, files, links = list_directory(samples_directory, '.', is_root=True)
    samples = [{'name': path, 'complete': is_complete(path)} for path in samples]

    template_index = env.get_template('index.html')
    return template_index.render(
        title='Overview',
        samples=samples,
        files=files,
        relpath=get_relative_path(f'{samples_directory}/overview.html', samples_directory)
    )


def serve_assembly(samples_directory, sample):
    print(f'>>>>>>>>>>>>>>> serving assembly {sample}')
    sample_dir = os.path.join(samples_directory, sample)

    assemblies = process_sample(sample, sample_dir, importers)

    dill.dump(assemblies, open(f"{sample_dir}/assemblies.pkl", 'wb'))

    return send_from_directory(
        directory=os.path.abspath(sample_dir),
        path='assemblies.html'
    )


@app.route('/<path:filepath>')
def serve_path(filepath):
    print(f">>>>>>>>>>>>>>> Requested {filepath}")
    full_path = os.path.abspath(os.path.join(samples_directory, filepath))
    dirname, basename = os.path.dirname(full_path), os.path.basename(full_path)

    # Import and serve assembly
    if basename == 'assemblies.html':
        sample = os.path.basename(dirname)
        return serve_assembly(samples_directory, sample)

    # Serve files
    if os.path.isfile(full_path):
        extension = full_path.rsplit('.', 1)[-1]
        mimetype = None
        if extension in ['fasta', 'tsv', 'gfa', 'gv']:
            mimetype = 'text/plain'
        return send_from_directory(dirname, basename, as_attachment=False, mimetype=mimetype)

    # Serve directories
    elif os.path.isdir(full_path):
        if not filepath.endswith('/'):
            return redirect(f'/{filepath}/')
        return list_directory(samples_directory, filepath)

    else:
        return "404 Not Found :(", 404


# PyWebIO application function
def pywebio_curate():
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
    assert os.path.isfile(f"{sample_dir}/assemblies.pkl")
    assemblies = dill.load(open(f"{sample_dir}/assemblies.pkl", 'rb'))
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
    save_curated(f"{sample_dir}/curated.fasta", headers)


# Add the PyWebIO endpoint
app.add_url_rule('/curate', 'webio_view', webio_view(pywebio_curate), methods=['GET', 'POST', 'OPTIONS'])


def run_server(
        samples_dir: str = 'data-pb',
        plugin_dir: str = None,
        port=8080,
        address='localhost',
        debug: bool = False
):
    global samples_directory
    samples_directory = samples_dir

    LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=LOGLEVEL)

    if plugin_dir is None:
        plugin_dir = os.environ.get('PLUGIN_DIR', './plugins-pb')
    assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"

    global importers
    importers = load_plugins(plugin_dir)


    handler = WSGIRequestHandler
    app.run(host=address, port=port, debug=debug)


def main():
    from fire import Fire
    Fire(run_server)


if __name__ == '__main__':
    main()
