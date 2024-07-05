from tempfile import TemporaryDirectory
import pandas as pd
import plotly.graph_objects as go
import plotly.figure_factory as ff
from gendiscalpy import GenDisCal

from assembler_tools.utils import AssemblyFailedException, MinorAssemblyException
from assembler_tools.Assembly import Assembly

gdc = GenDisCal()


def ani_clustermap(assemblies: [Assembly], output_file: str) -> str:
    try:
        df = calculate_distance_matrix(assemblies, output_file)
    except MinorAssemblyException as e:
        return f'<p class="error">{e}</p>'
    html = plot_clustermap(df)
    return html


def calculate_distance_matrix(assemblies: [Assembly], output_file: str):
    # Decision: create fasta for each contig or for each contig_group?
    fasta_files = []
    tmpdir = TemporaryDirectory()
    for assembly in assemblies:
        for contig_group in assembly.contig_groups:
            fasta_file = f'{tmpdir.name}/{assembly.assembler}#{contig_group.id}'
            fasta_files.append(fasta_file)
            with open(fasta_file, 'w') as fasta:
                for contig in contig_group.contigs:
                    fasta.write(f'>{assembly.assembler}#{contig.original_id}\n{contig.sequence}\n')

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


def plot_clustermap(distance_matrix: pd.DataFrame):
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
        showscale=False
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

# def plot_heatmap(distance_matrix: pd.DataFrame):
#     # Plot heatmap using plotly
#     fig = go.Figure(data=go.Heatmap(
#         z=distance_matrix.values,
#         x=distance_matrix.columns,
#         y=distance_matrix.columns,
#         colorscale='Viridis',
#         showscale=False
#     ))
#     fig.update_layout(width=800, height=800)
#     fig.update_layout(yaxis_side='right')  # Move labels to the right
#     fig.show()
# def plot_dendrogram(distance_matrix: pd.DataFrame):
#     # Plot dendrogram heatmap using plotly
#     fig = ff.create_dendrogram(
#         X=distance_matrix.values,
#         labels=distance_matrix.columns,
#         orientation='left',
#     )
#     fig.update_layout(width=800, height=800)
#     fig.show()
