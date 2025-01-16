import json
import time
import os
import glob
import subprocess
import tempfile
from datetime import timedelta

import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, gridspec, ticker

from assembly_curator.utils import human_bp
from assembly_curator.Contig import Contig
from assembly_curator.ContigGroup import ContigGroup

matplotlib.use('SVG')
mnl = ticker.MaxNLocator(nbins=4, prune='upper')
fmt = ticker.FuncFormatter(lambda x, pos: human_bp(x, decimals=0, zero_val='0'))


def run_minimap(ref, qry, out, params=[]):
    # Run minimap2 with the -o option to specify the output file
    cmd = ['minimap2', *params, ref.fasta, qry.fasta, '-o', out]
    proc = subprocess.run(cmd, capture_output=True)

    if proc.returncode != 0:
        stderr = proc.stderr.decode('utf-8')
        stdout = proc.stdout.decode('utf-8')
        raise RuntimeError(f"minimap2 failed with return code {proc.returncode}:\n{stderr}\n{stdout}")

    return cmd


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def parse_minimap_output(minimap_output):
    # Define column names
    columns = ["queryID", "queryLen", "queryStart", "queryEnd", "strand", "refID", "refLen", "refStart", "refEnd",
               "numResidueMatches", "lenAln", "mapQ", "alnType"]

    # Read the minimap output into a pandas DataFrame
    if os.stat(minimap_output).st_size==0:
        return pd.DataFrame(columns=columns)
    df = pd.read_csv(minimap_output, sep='\t', header=None)

    # Remove the columns that are not needed
    df = df.iloc[:, :len(columns)]
    df.columns = columns

    # if strand is '-', swap queryStart and queryEnd
    df.loc[df.strand == '-', ['queryStart', 'queryEnd']] = df.loc[df.strand == '-', ['queryEnd', 'queryStart']].values

    return df


def create_dotplot(
        df,
        title=None,
        ax=None, ax2=None,
        cg_ref: ContigGroup = None, cg_qry: ContigGroup = None,
        figsize: tuple[int | float, int | float] = (10, 10),
        show=False
):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    def plot_one_ax(rs, re, qs, qe, color, linewidth):
        ax.plot(
            [rs, re], [qs, qe],
            color=color,
            linewidth=linewidth
        )

    def plot_two_ax(rs, re, qs, qe, color, linewidth):
        ax.plot(
            [rs, re], [qs, qe],
            color=color,
            linewidth=linewidth
        )
        ax2.plot(
            [qs, qe], [rs, re],
            color=color,
            linewidth=linewidth
        )

    plotting_fn = plot_one_ax if ax2 is None else plot_two_ax

    df = df.sort_values('alnType', ascending=False)  # show primary alignments on top

    # Plot alignments
    for _, row in df.iterrows():
        plotting_fn(
            row['refStart'], row['refEnd'], row['queryStart'], row['queryEnd'],
            color='blue' if row['strand'] == '+' else 'red',
            linewidth=2 if row.alnType == 'tp:A:P' else 1
        )

    # Add separators
    running_sum = 0
    for c in cg_ref.contigs:
        running_sum += len(c)
        ax.axvline(running_sum, color='black', linewidth=1)
        if ax2:
            ax2.axhline(running_sum, color='black', linewidth=1)
    running_sum = 0
    for c in cg_qry.contigs:
        running_sum += len(c)
        ax.axhline(running_sum, color='black', linewidth=1)
        if ax2:
            ax2.axvline(running_sum, color='black', linewidth=1)

    if title is not None:
        ax.set_xlabel('Reference Position')
        ax.set_ylabel('Query Position')
        ax.set_title(title)
    if show:
        plt.show()

    return ax


def dotplot_minimap2(ax, ax2, out, cg_ref: ContigGroup, cg_qry: ContigGroup, params) -> plt.Axes:
    ax.set_xlim(0, len(cg_ref))
    ax.set_ylim(0, len(cg_qry))
    ax.invert_yaxis()

    if ax2:
        ax2.set_xlim(0, len(cg_qry))
        ax2.set_ylim(0, len(cg_ref))
        ax2.invert_yaxis()

    cmd = run_minimap(ref=cg_ref, qry=cg_qry, params=params, out=out)
    df = parse_minimap_output(out)
    return create_dotplot(df, ax=ax, ax2=ax2, cg_ref=cg_ref, cg_qry=cg_qry)


def format_axis(ax, i, j, n_ctgs, x_label: str, y_label: str):
    ax.yaxis.tick_right()

    if i == 0:
        ax.set_title(x_label, fontsize=12)
    if j == 0:
        ax.set_ylabel(y_label, fontsize=12)
        ax.yaxis.set_label_position('left')
    if i != n_ctgs - 1:
        ax.tick_params(axis='x', labelbottom=False)
    if j != n_ctgs - 1:
        ax.tick_params(axis='y', labelright=False)

    ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True)


def format_ticks(ax, len_seq1, len_seq2):
    # Formatter: turn pixel into readable bp
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)

    # Set font size, orientation
    ax.tick_params(axis='x', rotation=-90, labelsize=8)
    ax.tick_params(axis='y', labelsize=8)


def create_dotplots(
        workdir: str,
        figsize: (int, int) = (10, 10),
        title: str = None,
        output: str = 'dotplots.png',
        params=None
):
    if params is None:
        params = [
            '-t', '1',
            '-k', '28',
            '-N', '1000000',
            '-p', '0.000001',
            '--no-long-join',
        ]
    plt.rcParams['svg.fonttype'] = 'none'

    cgs = []
    for cg_json in glob.glob(os.path.join(workdir, '*.json')):
        cg = ContigGroup.from_json(cg_json)
        cg.fasta = cg_json.replace('.json', '.fasta')
        cgs.append(cg)
    n_cgs = len(cgs)

    # Calculate the relative width of each subplot
    total_length = sum([len(c) for c in cgs])
    subplot_ratios = [len(c) / total_length for c in cgs]

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(n_cgs, n_cgs, width_ratios=subplot_ratios, height_ratios=subplot_ratios)

    for i in range(n_cgs):
        cg_i = cgs[i]
        for j in range(n_cgs):
            cg_j = cgs[j]
            if i == j:
                ax = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_j.id} - {cg_i.id}')
                # start = time.time()
                dotplot_minimap2(ax, None, os.path.join(workdir, f'{i}-{j}.mm2'), cg_i, cg_j, params)
                format_axis(ax, i, j, n_cgs, cg_i.id, cg_j.id)
                format_ticks(ax, len(cg_i), len(cg_j))
                if len(cg_j.contigs) == 1 and cg_j.contigs[0].topology == 'circular':
                    ax.set_facecolor('lightgreen')  # circular
                # diff = time.time() - start
            elif i < j:
                # Add the dotplot to lower triangle
                ax1 = fig.add_subplot(gs[j, i], gid=f'dotplot - {cg_j.id} - {cg_i.id}')
                ax2 = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_i.id} - {cg_j.id}')
                # start = time.time()
                dotplot_minimap2(ax1, ax2, os.path.join(workdir, f'{i}-{j}.mm2'), cg_i, cg_j, params)
                format_axis(ax1, j, i, n_cgs, cg_i.id, cg_j.id)
                format_ticks(ax1, len(cg_i), len(cg_j))
                format_axis(ax2, i, j, n_cgs, cg_j.id, cg_i.id)
                format_ticks(ax2, len(cg_j), len(cg_i))
                # diff = time.time() - start
            else:
                continue
            # print('TIMING-MM2:', cg_i.id, n_cgs, cg_j.id, len(cg_j), diff)

    if title:
        fig.suptitle(title, fontsize=16)

    padding, space = 0.05, 0.1
    plt.subplots_adjust(left=padding, right=1 - padding, top=1 - padding, bottom=padding, wspace=space, hspace=space)
    plt.savefig(output, format='svg')
    plt.close()


def process_cluster(cluster_id, workdir, dotplot_outdir):
    create_dotplots(workdir, output=os.path.join(dotplot_outdir, f'{cluster_id}.svg'))
    return cluster_id

# with open('data/15_N/lja/assembly.fasta') as f:
#     seq = f.read()
#     # remove headers (>...) and newlines
#     seq = ''.join(s for s in seq.split('\n') if not s.startswith('>'))
# import time
#
#
# def time_it(params, ref=seq, qry=seq + reverse_complement(seq) + seq):
#     start = time.time()
#     res, cmd = run_minimap(ref=ref, qry=qry, params=params)
#     timedelta_minimap = time.time() - start
#     print(timedelta_minimap)
#
#     start = time.time()
#     df = parse_minimap_output(res)
#     create_dotplot(df, title=' '.join(params))
#     timedelta_parse_plot = time.time() - start
#     print(timedelta_parse_plot, df.shape, params)
#     return timedelta_minimap, timedelta_parse_plot, df
#
#
# res, cmd = run_minimap(
#     ref=seq, qry=seq + reverse_complement(seq) + seq,
#     params=['-t', '8',
#             '-k', '28',
#             '-N', '100000',
#             '-p', '0.00001',
#             ])
#
# df = parse_minimap_output(res)
# create_dotplot(df, title=' '.join(cmd[7:-4]), show=True)
# print(df.shape)

# time_it(
#     ['-t', '1',
#      '-k', '28',
#      '-N', '1000000',
#      '-p', '0.000001',
#      ])

# for k in ['15', '28']:
#     for n, p in [('1000', '0.001'), ('10000', '0.0001'), ('100000', '0.00001'), ('1000000', '0.000001')]:
#         time_it(['-t', '8',
#                  '-k', k,
#                  '-N', n,
#                  '-p', p,
#                  ])
# 6.731383800506592 (267, 13)  ['-t', '8', '-k', '15', '-N', '1000', '-p', '0.001']
# 6.639932870864868 (2235, 13) ['-t', '8', '-k', '15', '-N', '10000', '-p', '0.0001']
# 6.660287618637085 (4500, 13) ['-t', '8', '-k', '15', '-N', '100000', '-p', '0.00001']
# 6.642897605895996 (4500, 13) ['-t', '8', '-k', '15', '-N', '1000000', '-p', '0.000001']
# 4.756570100784302 (267, 13)  ['-t', '8', '-k', '28', '-N', '1000', '-p', '0.001']
# 4.712647914886475 (2241, 13) ['-t', '8', '-k', '28', '-N', '10000', '-p', '0.0001']
# 4.758112907409668 (3960, 13) ['-t', '8', '-k', '28', '-N', '100000', '-p', '0.00001']
# 4.743515253067017 (3960, 13) ['-t', '8', '-k', '28', '-N', '1000000', '-p', '0.000001']


# times = list(range(1, 11))
# timedeltas_minimap, timedeltas_parse_plot = [], []
# for t in times:
#     timedelta_minimap, timedelta_parse_plot, df = time_it(
#         ref=seq, qry=seq,
#         params=['-t', str(t),
#                 '-k', '28',
#                 '-N', '1000000',
#                 '-p', '0.000001',
#                 ])
#     timedeltas_minimap.append(timedelta_minimap)
#     timedeltas_parse_plot.append(timedelta_parse_plot)
#
# # plot the times
# fig, ax = plt.subplots()
# ax.plot(times, timedeltas_minimap, label='minimap2')
# ax.plot(times, timedeltas_parse_plot, label='parse and plot')
# plt.show()
# # basically no effect of -t parameter on runtime
