from io import StringIO
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.figure_factory as ff
from gendiscalpy import GenDisCal

from assembler_tools.utils import AssemblyFailedException, MinorAssemblyException
from assembler_tools.Assembly import Assembly

gdc = GenDisCal()


def ani_clustermap(assemblies: [Assembly], output_file: str) -> str:
    try:
        distance_matrix = calculate_distance_matrix(assemblies, output_file)
        length_matrix = calculate_length_matrix(distance_matrix, assemblies)
        html = plot_clustermap(distance_matrix, annot=length_matrix, annot_kws={"fontsize": 6}, fmt='', )
    except MinorAssemblyException as e:
        return f'<p class="error">{e}</p>'
    return html


def calculate_length_matrix(distance_matrix: pd.DataFrame, assemblies: [Assembly]):
    # Create a dataframe with same col/rows as distance_matrix
    # for each contig group (i, j), store the len(i) / len(j), i.e. the ratio of the lengths
    contig_group_to_len = {cg.id: len(cg) for assembly in assemblies for cg in assembly.contig_groups}

    length_matrix = pd.DataFrame(index=distance_matrix.index, columns=distance_matrix.columns)
    for i, assembly_i in enumerate(distance_matrix.index):
        for j, assembly_j in enumerate(distance_matrix.columns):
            if distance_matrix.iloc[i, j] < .5:
                frac = contig_group_to_len[assembly_i] / contig_group_to_len[assembly_j]
                length_matrix.iloc[i, j] = f'{frac:.2g}'
            else:
                length_matrix.iloc[i, j] = ''

    return length_matrix


def calculate_distance_matrix(assemblies: [Assembly], output_file: str):
    # Decision: create fasta for each contig or for each contig_group?
    fasta_files = []
    tmpdir = TemporaryDirectory()
    for assembly in assemblies:
        for contig_group in assembly.contig_groups:
            fasta_file = f'{tmpdir.name}/{contig_group.id}'
            fasta_files.append(fasta_file)
            with open(fasta_file, 'w') as fasta:
                for contig in contig_group.contigs:
                    fasta.write(f'>{contig_group.id}\n{contig.sequence}\n')

    # Raise exception if there is no or only one fasta file
    if len(fasta_files) == 0:
        raise AssemblyFailedException("No contigs to compare", severity='danger')
    elif len(fasta_files) == 1:
        raise MinorAssemblyException("Only one contig, cannot create distance matrix")

    # Compute pairwise distance matrix
    try:
        df = gdc.run(*fasta_files, distance_matrix=True)
        df.to_csv(output_file, sep='\t')
    except Exception as e:
        print(f"Error: {e}")
        raise e

    tmpdir.cleanup()
    return df


def plot_clustermap_with_seaborn(distance_matrix: pd.DataFrame, **kwargs) -> str:
    plt.rcParams['svg.fonttype'] = 'none'

    # Generate the clustermap
    clustermap_fig = sns.clustermap(
        distance_matrix,
        figsize=(10, 10),
        cmap='bone_r', vmin=0, vmax=1,
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
