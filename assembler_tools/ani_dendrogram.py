import logging
from io import StringIO
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.figure_factory as ff
import pyskani

from assembler_tools import ContigGroup
from assembler_tools.utils import AssemblyFailedException, MinorAssemblyException
from assembler_tools.Assembly import Assembly

matplotlib.use('agg')  # non-interactive backend that can only write to files


def ani_clustermap(assemblies: [Assembly], fname: str, cutoff: float = .9) -> (pd.DataFrame, str):
    assert len(assemblies) > 0, "No assemblies to cluster!"
    similarity_matrix = calculate_similarity_matrix(assemblies)

    length_matrix = calculate_length_matrix(similarity_matrix, assemblies, label_cutoff=cutoff)
    sample_to_cluster = plot_clustermap(
        similarity_matrix,
        annot=length_matrix, annot_kws={"fontsize": 6}, fmt='',
        vmin=cutoff, vmax=1,
        fname=fname
    )
    return similarity_matrix, sample_to_cluster


def add_cluster_info_to_assemblies(assemblies, sample_to_cluster):
    for assembly in assemblies:
        for contig_group in assembly.contig_groups:
            cluster_info = sample_to_cluster[contig_group.id]
            contig_group.cluster_id = cluster_info['cluster']
            contig_group.set_cluster_color(*cluster_info['color'])


def calculate_length_matrix(
        similarity_matrix: pd.DataFrame,
        assemblies: [Assembly],
        label_cutoff: float = .9
) -> pd.DataFrame:
    contig_groups = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}
    contig_group_to_len = {id: len(cg) for id, cg in contig_groups.items()}

    # Create a dataframe with same col/rows as similarity_matrix
    # for each contig group (i, j), store the len(i) / len(j), i.e. the ratio of the lengths
    length_matrix = pd.DataFrame(index=similarity_matrix.index, columns=similarity_matrix.columns)
    for i, cg_id_i in enumerate(similarity_matrix.index):
        for j, cg_id_j in enumerate(similarity_matrix.columns):
            if cg_id_i == cg_id_j:
                cg: ContigGroup = contig_groups[cg_id_i]
                length_matrix.iloc[i, j] = cg.topology_or_n_contigs(short=True)
                continue
            if similarity_matrix.iloc[i, j] > label_cutoff:
                frac = contig_group_to_len[cg_id_i] / contig_group_to_len[cg_id_j]
                length_matrix.iloc[i, j] = f'{frac:.2g}'
            else:
                length_matrix.iloc[i, j] = ''

    return length_matrix


def calculate_similarity_matrix(assemblies: [Assembly]):
    """Decision: create fasta for each contig or for each contig_group"""
    # Create a temporary directory for the pyskani sketches
    pyskani_db = pyskani.Database()

    # Import all contig groups into the skani database
    for assembly in assemblies:
        for contig_group in assembly.contig_groups:
            pyskani_db.sketch(contig_group.id, *contig_group.encode_sequences())

    contig_group_ids = [cg.id for assembly in assemblies for cg in assembly.contig_groups]
    similarity_matrix = pd.DataFrame(index=contig_group_ids, columns=contig_group_ids, data=0.0)

    for assembly in assemblies:
        for contig_group in assembly.contig_groups:
            hits = pyskani_db.query(contig_group.id, *contig_group.encode_sequences())
            for hit in hits:
                similarity_matrix.at[contig_group.id, hit.reference_name] = hit.identity

    # This similarity matrix is not always symmetric. We enforce this to get a nice diagonal after clustering.
    similarity_matrix = (similarity_matrix + similarity_matrix.T) / 2

    return similarity_matrix


def plot_clustermap(similarity_matrix: pd.DataFrame, fname: str, **kwargs) -> (str, np.ndarray):
    plt.rcParams['svg.fonttype'] = 'none'
    assert (similarity_matrix.index == similarity_matrix.columns).all(), \
        "Similarity matrix must be square:\nindex={similarity_matrix.index}\ncolumns={similarity_matrix.columns}"

    # Compute the linkage matrix manually so it can be returned
    linkage_matrix = sch.linkage(similarity_matrix, method='average')

    sample_to_cluster = group_by_linkage(similarity_matrix, linkage_matrix)
    row_col_colors = [sample_to_cluster[cg]['color'] for cg in similarity_matrix.index]

    # Generate the clustermap
    clustermap_fig = sns.clustermap(
        similarity_matrix,
        figsize=(10, 10),
        cmap='mako_r',
        linewidths=5,
        row_cluster=True,
        col_cluster=True,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        row_colors=row_col_colors,
        col_colors=row_col_colors,
        **kwargs
    )

    svg_buffer = StringIO()
    clustermap_fig.savefig(svg_buffer, format='svg')
    plt.close()
    svg_content = svg_buffer.getvalue()
    svg_buffer.close()

    # Convert to JSON, add to SVG
    json_data = similarity_matrix.to_json()  # Reordering not necessary because json uses keys, not order
    svg_content = svg_content.rsplit('</svg>', 1)
    assert svg_content[1].strip() == '', f"Unexpected content after closing svg tag: {svg_content[1]}"
    svg_content = f'{svg_content[0]}\n<script type="application/json" id="ani-matrix-data">{json_data}</script></svg>'

    with open(fname, 'w') as f:
        f.write(svg_content)

    return sample_to_cluster


def group_by_linkage(similarity_matrix, linkage_matrix: np.ndarray, threshold: float = .95) -> {int: [str]}:
    """Form flat clusters based on threshold"""
    cluster_group = sch.fcluster(linkage_matrix, t=threshold, criterion='distance')
    cluster_colors = sns.color_palette('tab20', n_colors=len(set(cluster_group)))

    return {
        i: {'cluster': int(g), 'color': cluster_colors[int(g) - 1]}
        for i, g
        in zip(similarity_matrix.index, cluster_group)
    }
