import time
import os
import io
import subprocess
import tempfile
from datetime import timedelta

import matplotlib.pyplot as plt
import pandas as pd

from assembler_tools.utils import human_bp


def run_minimap(ref, qry, params=[]):
    with tempfile.TemporaryDirectory() as tempdir:
        ref_path = os.path.join(tempdir, 'ref.fna')
        qry_path = os.path.join(tempdir, 'qry.fna')
        outfile_path = os.path.join(tempdir, 'outfile.mm2')

        with open(ref_path, 'w') as ref_f, open(qry_path, 'w') as qry_f:
            ref_f.write('>ref\n' + ref)
            qry_f.write('>qry\n' + qry)

        # Run minimap2 with the -o option to specify the output file
        cmd = ['/home/thomas/PycharmProjects/assembler-tools/bin/minimap2',
               *params,
               ref_path, qry_path,
               '-o', outfile_path]
        minimap_out = subprocess.run(cmd, capture_output=True)

        if minimap_out.returncode != 0:
            stderr = minimap_out.stderr.decode('utf-8')
            stdout = minimap_out.stdout.decode('utf-8')
            raise RuntimeError(f"minimap2 failed with return code {minimap_out.returncode}:\n{stderr}\n{stdout}")

        # Read the contents of the output file
        with open(outfile_path, 'r') as outfile:
            result = outfile.read()

        return result, cmd


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def parse_minimap_output(minimap_output):
    # Define column names
    columns = ["queryID", "queryLen", "queryStart", "queryEnd", "strand", "refID", "refLen", "refStart", "refEnd",
               "numResidueMatches", "lenAln", "mapQ", "alnType"]

    # Read the minimap output into a pandas DataFrame
    df = pd.read_csv(io.StringIO(minimap_output), sep='\t', header=None)

    # Remove the columns that are not needed
    df = df.iloc[:, :len(columns)]
    df.columns = columns

    # if strand is '-', swap queryStart and queryEnd
    df.loc[df.strand == '-', ['queryStart', 'queryEnd']] = df.loc[df.strand == '-', ['queryEnd', 'queryStart']].values

    return df


def create_dotplot(df, title=None, ax=None, ax2=None, sep_ref=[], sep_qry=[], show=False):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

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
    for i in sep_ref:
        ax.axvline(i, color='black', linewidth=1)
        if ax2:
            ax2.axhline(i, color='black', linewidth=1)
    for i in sep_qry:
        ax.axhline(i, color='black', linewidth=1)
        if ax2:
            ax2.axvline(i, color='black', linewidth=1)

    if title is not None:
        ax.set_xlabel('Reference Position')
        ax.set_ylabel('Query Position')
        ax.set_title(title)
    if show:
        plt.show()

    return ax


def dotplot_minimap2(ax, ax2, ref, qry, sep_ref, sep_qry, params) -> plt.Axes:
    ax.set_xlim(0, len(ref))
    ax.set_ylim(0, len(qry))
    ax.invert_yaxis()

    if ax2:
        ax2.set_xlim(0, len(qry))
        ax2.set_ylim(0, len(ref))
        ax2.invert_yaxis()

    res, cmd = run_minimap(ref=ref, qry=qry, params=params)
    df = parse_minimap_output(res)
    return create_dotplot(df, ax=ax, ax2=ax2, sep_ref=sep_ref, sep_qry=sep_qry)


from assembler_tools.ContigGroup import ContigGroup
from matplotlib import pyplot as plt, gridspec, ticker

mnl = ticker.MaxNLocator(nbins=4, prune='upper')


def format_axis(ax, i, j, n_ctgs, x_label, y_label, bp_per_pixel):
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


def format_ticks(ax, seq1, seq2, bp_per_pixel):
    # Set ticks relative to bp, not pixels
    seq1_ticks = [t / bp_per_pixel for t in mnl.tick_values(0, len(seq1))]
    seq2_ticks = [t / bp_per_pixel for t in mnl.tick_values(0, len(seq2))]
    ax.set_xticks(seq1_ticks)
    ax.set_yticks(seq2_ticks)

    # Formatter: turn pixel into readable bp
    fmt = ticker.FuncFormatter(lambda x, pos: human_bp(x * bp_per_pixel, decimals=0, zero_val='0'))
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)

    # Set font size, orientation
    ax.tick_params(axis='x', rotation=-90, labelsize=8)
    ax.tick_params(axis='y', labelsize=8)


def create_dotplots(
        cgs: [ContigGroup],
        figsize: (int, int) = (10, 10),
        title: str = None,
        output: str = 'dotplots.png',
        kmer=12,
        bp_per_pixel: float = 17.651,
        background_colour='white',
        color_fwd='mediumblue', color_rev='firebrick',
        line_color='black',
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

    # Calculate the total length of all sequences
    seq_lengths = [sum(len(c.sequence) for c in cg.contigs) for cg in cgs]
    total_length = sum(seq_lengths)

    # Calculate the relative width of each subplot
    relative_widths = [length / total_length for length in seq_lengths]

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(len(cgs), len(cgs), width_ratios=relative_widths, height_ratios=relative_widths)

    for i, cg_i in enumerate(cgs):
        seq1 = ''.join(c.sequence for c in cg_i.contigs)
        sep1 = [0] + [len(contig.sequence) for contig in cg_i.contigs]

        for j, cg_j in enumerate(cgs):
            seq2 = ''.join(c.sequence for c in cg_j.contigs)
            sep2 = [0] + [len(contig.sequence) for contig in cg_j.contigs]

            if i == j:
                ax = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_i.id} - {cg_j.id}')
                start = time.time()
                dotplot_minimap2(ax, None, seq1, seq2, sep1, sep2, params)
                format_axis(ax, j, i, len(cgs), cg_i.id, cg_j.id, 1)
                format_ticks(ax, seq1, seq2, 1)
                diff = time.time() - start
            elif i < j:
                # Add the dotplot to lower triangle
                ax1 = fig.add_subplot(gs[j, i], gid=f'dotplot - {cg_j.id} - {cg_i.id}')
                ax2 = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_i.id} - {cg_j.id}')
                start = time.time()
                dotplot_minimap2(ax1, ax2, seq1, seq2, sep1, sep2, params)
                format_axis(ax1, j, i, len(cgs), cg_i.id, cg_j.id, 1)
                format_ticks(ax1, seq1, seq2, 1)
                format_axis(ax2, i, j, len(cgs), cg_j.id, cg_i.id, 1)
                format_ticks(ax2, seq2, seq1, 1)
                diff = time.time() - start
            else:
                continue
            print('TIMING-MM2:', cg_i.id, len(cg_i), cg_j.id, len(cg_j), diff)

    if title:
        fig.suptitle(title, fontsize=16)

    padding, space = 0.05, 0.1
    plt.subplots_adjust(left=padding, right=1 - padding, top=1 - padding, bottom=padding, wspace=space, hspace=space)
    plt.savefig(output, format='svg')
    plt.close()

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
