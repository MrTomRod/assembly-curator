import logging
from io import StringIO
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.figure_factory as ff
import pyskani

from assembler_tools import ContigGroup
from assembler_tools.utils import AssemblyFailedException, MinorAssemblyException
from assembler_tools.Assembly import Assembly


def ani_clustermap(assemblies: [Assembly], cutoff: float = .9) -> str:
    try:
        distance_matrix = calculate_similarity_matrix(assemblies)
        length_matrix = calculate_length_matrix(distance_matrix, assemblies, label_cutoff=cutoff)
        html = plot_clustermap(distance_matrix,
                               annot=length_matrix, annot_kws={"fontsize": 6}, fmt='',
                               vmin=cutoff, vmax=1)
    except MinorAssemblyException as e:
        logging.warning(f"Failed to compute ANI matrix: {e}")
        return None, f'<p class="error">{e}</p>'
    return distance_matrix, html


def calculate_length_matrix(
        distance_matrix: pd.DataFrame,
        assemblies: [Assembly],
        label_cutoff: float = .9
) -> pd.DataFrame:
    contig_groups = {cg.id: cg for assembly in assemblies for cg in assembly.contig_groups}
    contig_group_to_len = {id: len(cg) for id, cg in contig_groups.items()}

    # Create a dataframe with same col/rows as distance_matrix
    # for each contig group (i, j), store the len(i) / len(j), i.e. the ratio of the lengths
    length_matrix = pd.DataFrame(index=distance_matrix.index, columns=distance_matrix.columns)
    for i, cg_id_i in enumerate(distance_matrix.index):
        for j, cg_id_j in enumerate(distance_matrix.columns):
            if cg_id_i == cg_id_j:
                cg: ContigGroup = contig_groups[cg_id_i]
                length_matrix.iloc[i, j] = cg.topology_or_n_contigs(short=True)
                continue
            if distance_matrix.iloc[i, j] > label_cutoff:
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

    if len(similarity_matrix) < 2:
        raise MinorAssemblyException("Not enough contig groups to create a similarity matrix")

    return similarity_matrix


def plot_clustermap_with_seaborn(distance_matrix: pd.DataFrame, **kwargs) -> str:
    plt.rcParams['svg.fonttype'] = 'none'

    # Generate the clustermap
    clustermap_fig = sns.clustermap(
        distance_matrix,
        figsize=(10, 10),
        cmap='mako_r',
        linewidths=5,
        **kwargs
    )

    svg_buffer = StringIO()
    clustermap_fig.savefig(svg_buffer, format='svg')
    svg_content = svg_buffer.getvalue()
    svg_buffer.close()

    # Convert to JSON, add to SVG
    json_data = distance_matrix.to_json()  # Reordering not necessary because json uses keys, not order
    svg_content = svg_content.rsplit('</svg>', 1)
    assert svg_content[1].strip() == '', f"Unexpected content after closing svg tag: {svg_content[1]}"
    svg_content = f'{svg_content[0]}\n<script type="application/json" id="ani-matrix-data">{json_data}</script></svg>'

    return svg_content


plot_clustermap = plot_clustermap_with_seaborn


def plot_clustermap_plotly(distance_matrix: pd.DataFrame) -> str:
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(distance_matrix, orientation='bottom', labels=distance_matrix.columns)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(distance_matrix, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    heat_data = distance_matrix.values
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = go.Heatmap(
        x=dendro_leaves,
        y=dendro_leaves,
        z=heat_data,
        colorscale='Blues',
        showscale=False,
        customdata=np.indices(heat_data.shape),  # Add custom data for each cell
        hoverinfo='none'  # Disable hover info
    )

    # Align the heatmap and the dendrogram
    heatmap['x'] = fig['layout']['xaxis']['tickvals']
    heatmap['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    fig.add_trace(heatmap)

    # Edit Layout
    fig.update_layout({
        'width': 1000, 'height': 1000,
        'showlegend': False, 'hovermode': 'closest',
    })

    fig.update_layout(
        xaxis={  # xaxis of heatmap and dendrogram above
            'domain': [.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'tickangle': 90,  # easier to read
        },
        xaxis_autorange='reversed',
        yaxis={  # yaxis of heatmap and dendrogram on the left
            'domain': [0, .85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'tickvals': list(range(5, len(dendro_leaves) * 10, 10)),  # yaxis tick positions
            'ticktext': [distance_matrix.columns[i] for i in dendro_leaves],  # yaxis tick labels
            'tickangle': 0,  # fixed tickangle
        },
        yaxis_side='right',  # move yaxis ticks to the right
        xaxis2={  # xaxis of dendrogram on the left
            'domain': [0, .15],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        },
        yaxis2={  # yaxis of dendrogram above heatmap
            'domain': [.825, .975],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        })

    # fig.show()
    # exit(0)
    html = fig.to_html(full_html=False, div_id='ani-clustermap')
    return html
