import os
import re
import json
import logging
from io import StringIO

import pyskani
from typing import Type
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib as mpl

mpl.use('SVG')
import matplotlib.pyplot as plt

from assembly_curator.AssemblyImporter import AssemblyImporter


def get_assembly(sample_dir, assembly_importer: AssemblyImporter):
    fasta = os.path.join(sample_dir, assembly_importer.assembly_dir, assembly_importer.assembly)
    if os.path.isfile(fasta) and os.stat(fasta).st_size > 0:
        return fasta


def get_sequences(fasta: str):
    with open(fasta) as f:
        data = f.read()
    encoded_sequences = []
    for entry in data.split('>'):
        entry = entry.strip()
        if entry == '': continue  # start or end of FASTA file
        contig_header, sequence = entry.split('\n', 1)
        sequence = sequence.replace('\n', '')
        assert set(sequence) == {'A', 'T', 'C', 'G'}
        encoded_sequences.append(sequence.encode('ascii'))
    return encoded_sequences


def calculate_similarity_matrix(samples_dir, importers):
    logging.info(f"Calculating Phylogenetic tree...")

    pyskani_db = pyskani.Database()

    # create db
    samples = os.listdir(samples_dir)
    samples_to_sequence = {}
    for i, sample in enumerate(samples):
        logging.info(f"Creating pyskani db: {i}/{len(samples)}")
        for importer in importers:
            fasta = get_assembly(os.path.join(samples_dir, sample), importer)
            if fasta is not None:
                samples_to_sequence[sample] = get_sequences(fasta)
                pyskani_db.sketch(sample, *samples_to_sequence[sample])

    # initiate similarity_matrix
    similarity_matrix = pd.DataFrame(index=samples_to_sequence.keys(), columns=samples_to_sequence.keys(), data=0.0)

    # fill values
    for i, (sample, sequence) in enumerate(samples_to_sequence.items()):
        logging.info(f"Querying pyskani db: {i}/{len(samples)}")
        hits = pyskani_db.query(sample, *sequence)
        for hit in hits:
            similarity_matrix.at[sample, hit.reference_name] = hit.identity

    # This similarity matrix is not always symmetric. We enforce this to get a nice diagonal after clustering.
    similarity_matrix = (similarity_matrix + similarity_matrix.T) / 2

    return similarity_matrix


def calculate_phylogenetic_tree(
        importers: [Type[AssemblyImporter]],
        samples_dir: str,
        force_rerun: bool = False
):
    similarity_matrix_file = os.path.join(samples_dir, 'similarity_matrix.tsv')
    similarity_matrix_plot = os.path.join(samples_dir, 'similarity_matrix.svg')
    similarity_matrix_plot_order = os.path.join(samples_dir, 'similarity_matrix.json')
    if (not force_rerun
            and os.path.isfile(similarity_matrix_file)
            and os.path.isfile(similarity_matrix_plot)
            and os.path.isfile(similarity_matrix_plot_order)):
        logging.info(f"Skipping phylogenetic tree calculation.")
        with open(similarity_matrix_plot_order) as f:
            samples_sorted = json.load(f)
        return samples_sorted

    similarity_matrix = calculate_similarity_matrix(samples_dir, importers)
    # similarity_matrix = pd.read_csv(similarity_matrix_file, sep='\t', index_col=0)
    similarity_matrix.to_csv(similarity_matrix_file, sep='\t')

    low_score = similarity_matrix[similarity_matrix > 0].min().min()
    similarity_matrix.replace(0, max(0, low_score - 0.02), inplace=True)

    dendrogram_params = plot_dendrogram(similarity_matrix, similarity_matrix_plot)

    samples_sorted = list(reversed(dendrogram_params['ivl']))
    with open(similarity_matrix_plot_order, 'w') as f:
        json.dump(samples_sorted, f)

    return samples_sorted


def plot_dendrogram(similarity_matrix, fname: str):
    plt.rcParams['svg.fonttype'] = 'none'

    # calculate plot proportions
    content_height = len(similarity_matrix.index) / 3  # height dependent on number of compounds

    # create matplotlib figure
    plt.close()
    fig = plt.figure(figsize=(2, content_height))  # dpi irrelevant if mpl.use('SVG')
    ax = fig.add_axes([0, 0, .95, 1])  # [left, bottom, width, height]

    # Compute the linkage matrix manually so it can be returned
    linkage_matrix = sch.linkage(similarity_matrix, method='single')

    dendrogram_params = sch.dendrogram(
        linkage_matrix,
        orientation='left',
        labels=similarity_matrix.index,
        no_labels=True,
        ax=ax,
        color_threshold=0,
        above_threshold_color='k'  # greyscale
    )

    ax.grid(False)
    plt.box(False)
    ax.tick_params(
        axis='both', which='both',
        bottom=False, top=False, left=False, right=False,
        labelleft=False,
    )

    # plt.savefig(fname, format='svg')
    svg_buffer = StringIO()
    fig.savefig(svg_buffer, format='svg')
    plt.close()
    svg_content = svg_buffer.getvalue()
    svg_buffer.close()

    # Disable preserveAspectRatio (https://www.sarasoueidan.com/blog/svg-object-fit/)
    pattern = r'width="[^"]*"[^>]*viewBox'
    replacement = 'preserveAspectRatio="xMinYMin slice" viewBox'
    svg_content = re.sub(pattern, replacement, svg_content, count=1)

    with open(fname, 'w') as f:
        f.write(svg_content)

    return dendrogram_params
