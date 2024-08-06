import logging
import collections
from PIL import Image, ImageDraw, ImageFont, ImageColor
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import pyplot as plt, gridspec, ticker

from assembler_tools.utils import human_bp
from assembler_tools.ContigGroup import ContigGroup

mnl = ticker.MaxNLocator(nbins=4, prune='upper')

REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def reverse_complement(seq):
    return ''.join([REV_COMP_DICT.get(base, 'N') for base in seq][::-1])


def extract_contig(fasta_path: str, contig_id: str) -> str:
    sequence = ''
    with open(fasta_path) as f:
        append = False
        for line in f:
            if line.startswith('>'):
                id = line.split(' ', 1)[0].strip('>\n')
                append = (id == contig_id)
                if append:
                    sequence = ''
                continue
            if append:
                sequence += line.strip()
    if not sequence:
        raise ValueError(f'Error: Contig {contig_id} not found in {fasta_path}')
    return sequence


def load_contig_sequence(fasta_and_contig: str) -> dict[str, str]:
    return extract_contig(*fasta_and_contig.rsplit(':', 1))


class DotplotConfig:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return '\n'.join([f'{k}: {v}' for k, v in self.__dict__.items()])


def add_separators(image, draw, sep_a, sep_b, config):
    for i in sep_b:
        i = int(round(i / config.bp_per_pixel))
        draw.line([(0, i), (image.width, i)], fill=config.line_color)
    for i in sep_a:
        i = int(round(i / config.bp_per_pixel))
        draw.line([(i, 0), (i, image.height)], fill=config.line_color)


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
        line_color='black'
):
    color_fwd = ImageColor.getrgb(color_fwd)
    color_rev = ImageColor.getrgb(color_rev)
    line_color = ImageColor.getrgb(line_color)
    config = DotplotConfig(**{k: v for k, v in locals().items()
                              if k not in ['fasta_paths_and_contigs', 'output', 'sep_a', 'sep_b']})

    plt.rcParams['svg.fonttype'] = 'none'

    fig = plt.figure(figsize=figsize)

    # Calculate the total length of all sequences
    seq_lengths = [sum(len(c.sequence) for c in cg.contigs) for cg in cgs]
    total_length = sum(seq_lengths)

    # Calculate the relative width of each subplot
    relative_widths = [length / total_length for length in seq_lengths]

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(len(cgs), len(cgs), width_ratios=relative_widths, height_ratios=relative_widths)

    for i, cg_i in enumerate(cgs):
        seq1 = ''.join(c.sequence for c in cg_i.contigs)
        sep1 = [0] + [len(contig.sequence) for contig in cg_i.contigs]

        a_forward_kmers, a_reverse_kmers = get_all_kmer_positions(config.kmer, seq1)
        for j, cg_j in enumerate(cgs):
            seq2 = ''.join(c.sequence for c in cg_j.contigs)
            sep2 = [0] + [len(contig.sequence) for contig in cg_j.contigs]

            if i == j:
                im, draw = create_dotplot(seq1, seq2, sep1, sep2, a_forward_kmers, a_reverse_kmers, config)
                ax = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_i.id} - {cg_j.id}')
                ax.imshow(im)
                format_axis(ax, j, i, len(cgs), cg_i.id, cg_j.id, bp_per_pixel)
                format_ticks(ax, seq1, seq2, bp_per_pixel)
            elif i < j:
                # Add the dotplot to upper triangle
                im, draw = create_dotplot(seq1, seq2, sep1, sep2, a_forward_kmers, a_reverse_kmers, config)
                ax = fig.add_subplot(gs[j, i], gid=f'dotplot - {cg_j.id} - {cg_i.id}')
                ax.imshow(im)
                format_axis(ax, j, i, len(cgs), cg_i.id, cg_j.id, bp_per_pixel)
                format_ticks(ax, seq1, seq2, bp_per_pixel)

                # Add the dotplot to lower triangle, simply rotating it 90 degrees
                ax = fig.add_subplot(gs[i, j], gid=f'dotplot - {cg_i.id} - {cg_j.id}')
                ax.imshow(im.transpose(Image.Transpose.ROTATE_90))
                format_axis(ax, i, j, len(cgs), cg_j.id, cg_i.id, bp_per_pixel)
                format_ticks(ax, seq2, seq1, bp_per_pixel)
            else:
                continue

    del a_forward_kmers, a_reverse_kmers

    if title:
        fig.suptitle(title, fontsize=16)

    padding, space = 0.05, 0.1
    plt.subplots_adjust(left=padding, right=1 - padding, top=1 - padding, bottom=padding, wspace=space, hspace=space)
    plt.savefig(output, format='svg')


def create_dotplot(
        seq_a: str, seq_b: str,
        sep_a, sep_b,
        a_forward_kmers, a_reverse_kmers,
        config: DotplotConfig
):
    image = Image.new(
        mode='RGB',
        size=[int(round(len(seq) / config.bp_per_pixel)) for seq in [seq_a, seq_b]],
        color=config.background_colour
    )
    draw = ImageDraw.Draw(image)

    def draw_dot(i, j, colour):
        draw.point((i, j), fill=colour)
        draw.point((i + 1, j), fill=colour)
        draw.point((i - 1, j), fill=colour)
        draw.point((i, j + 1), fill=colour)
        draw.point((i, j - 1), fill=colour)

    for j in range(len(seq_b) - config.kmer + 1):
        j_pixel = int(round(j / config.bp_per_pixel))
        k = seq_b[j:j + config.kmer]
        if k in a_reverse_kmers:
            for i in a_reverse_kmers[k]:
                i_pixel = int(round(i / config.bp_per_pixel))
                draw_dot(i_pixel, j_pixel, config.color_rev)

        if k in a_forward_kmers:
            for i in a_forward_kmers[k]:
                i_pixel = int(round(i / config.bp_per_pixel))
                draw_dot(i_pixel, j_pixel, config.color_fwd)

    add_separators(image, draw, sep_a, sep_b, config)

    return image, draw


def get_all_kmer_positions(kmer_size, seq):
    forward_kmers, reverse_kmers = collections.defaultdict(list), collections.defaultdict(list)
    rev_comp_seq = reverse_complement(seq)
    seq_len = len(seq) - kmer_size + 1
    for i in range(seq_len):
        k = seq[i:i + kmer_size]
        forward_kmers[k].append(i)
        k = rev_comp_seq[i:i + kmer_size]
        reverse_kmers[k].append(seq_len - i - 1)
    assert len(forward_kmers) < len(seq)
    assert len(reverse_kmers) < len(seq)
    return forward_kmers, reverse_kmers

# import logging
# import collections
# from PIL import Image, ImageDraw, ImageFont, ImageColor
#
# REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
#                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
#                  'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
#                  'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
#                  '?': '?'}
#
#
# def reverse_complement(seq):
#     return ''.join([REV_COMP_DICT.get(base, 'N') for base in seq][::-1])
#
#
# def extract_contig(fasta_path: str, contig_id: str) -> str:
#     sequence = ''
#     with open(fasta_path) as f:
#         append = False
#         for line in f:
#             if line.startswith('>'):
#                 id = line.split(' ', 1)[0].strip('>\n')
#                 append = (id == contig_id)
#                 if append:
#                     sequence = ''
#                 continue
#             if append:
#                 sequence += line.strip()
#     if not sequence:
#         raise ValueError(f'Error: Contig {contig_id} not found in {fasta_path}')
#     return sequence
#
#
# def load_contig_sequence(fasta_and_contig: str) -> dict[str, str]:
#     return extract_contig(*fasta_and_contig.rsplit(':', 1))
#
#
# class DotplotConfig:
#     def __init__(self, **kwargs):
#         self.__dict__.update(kwargs)
#
#     def __str__(self):
#         return '\n'.join([f'{k}: {v}' for k, v in self.__dict__.items()])
#
#
# def dotplot(
#         fasta_and_contig_1: str, fasta_and_contig_2: str,
#         sep_a: [int] = None, sep_b: [int] = None,
#         output: str = 'dotplots.png',
#         kmer=10,
#         bp_per_pixel: float = 20.,
#         background_colour='white',
#         color_fwd='mediumblue', color_rev='firebrick',
#         line_color='black'
# ):
#     color_fwd = ImageColor.getrgb(color_fwd)
#     color_rev = ImageColor.getrgb(color_rev)
#     line_color = ImageColor.getrgb(line_color)
#     if sep_a is None:
#         sep_a = []
#     if sep_b is None:
#         sep_b = []
#     config = DotplotConfig(**{k: v for k, v in locals().items()
#                               if k not in ['fasta_paths_and_contigs', 'output', 'sep_a', 'sep_b']})
#
#     seq_a = load_contig_sequence(fasta_and_contig_1)
#     seq_b = load_contig_sequence(fasta_and_contig_2)
#
#     image, draw = create_dotplots(seq_a, seq_b, config)
#
#     add_separators(image, draw, sep_a, sep_b, config)
#
#     image.save(output)
#     print(f'Saving dotplots to: {output}')
#
#
# def add_separators(image, draw, sep_a, sep_b, config):
#     for i in sep_a:
#         i = int(round(i / config.bp_per_pixel))
#         draw.line([(0, i), (image.width, i)], fill=config.line_color)
#     for i in sep_b:
#         i = int(round(i / config.bp_per_pixel))
#         draw.line([(i, 0), (i, image.height)], fill=config.line_color)
#
#
# def create_dotplots(
#         seq_a: str, seq_b: str,
#         config: DotplotConfig
# ):
#     image = Image.new(
#         mode='RGB',
#         size=[int(round(len(seq) / config.bp_per_pixel)) for seq in [seq_a, seq_b]],
#         color=config.background_colour
#     )
#     draw = ImageDraw.Draw(image)
#
#     def draw_dot(i, j, colour):
#         draw.point((i, j), fill=colour)
#         draw.point((i + 1, j), fill=colour)
#         draw.point((i - 1, j), fill=colour)
#         draw.point((i, j + 1), fill=colour)
#         draw.point((i, j - 1), fill=colour)
#
#     a_forward_kmers, a_reverse_kmers = get_all_kmer_positions(config.kmer, seq_a)
#
#     for j in range(len(seq_b) - config.kmer + 1):
#         j_pixel = int(round(j / config.bp_per_pixel))
#         k = seq_b[j:j + config.kmer]
#         if k in a_reverse_kmers:
#             for i in a_reverse_kmers[k]:
#                 i_pixel = int(round(i / config.bp_per_pixel))
#                 draw_dot(i_pixel, j_pixel, config.color_rev)
#
#         if k in a_forward_kmers:
#             for i in a_forward_kmers[k]:
#                 i_pixel = int(round(i / config.bp_per_pixel))
#                 draw_dot(i_pixel, j_pixel, config.color_fwd)
#
#     return image, draw
#
#
# def get_all_kmer_positions(kmer_size, seq):
#     forward_kmers, reverse_kmers = collections.defaultdict(list), collections.defaultdict(list)
#     rev_comp_seq = reverse_complement(seq)
#     seq_len = len(seq) - kmer_size + 1
#     for i in range(seq_len):
#         k = seq[i:i + kmer_size]
#         forward_kmers[k].append(i)
#         k = rev_comp_seq[i:i + kmer_size]
#         reverse_kmers[k].append(seq_len - i - 1)
#     assert len(forward_kmers) < len(seq)
#     assert len(reverse_kmers) < len(seq)
#     return forward_kmers, reverse_kmers
#
#
# if __name__ == '__main__':
#     from fire import Fire
#
#     Fire(dotplot)
