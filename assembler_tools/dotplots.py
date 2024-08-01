from rrwick_dotplots import *


def dotplot_api(
        seq_names,
        seqs,
        output: str,
        res=800, kmer=10,
        initial_top_left_gap=0.1, border_gap=0.015, between_seq_gap=0.0125, outline_width=0.0015,
        text_gap=0.005, max_font_size=0.025,
        background_colour='white', self_vs_self_colour='lightgrey', self_vs_other_colour='whitesmoke',
        text_colour='black', forward_strand_dot_colour='mediumblue', reverse_strand_dot_colour='firebrick'
):
    config = DotplotConfig(**{k: v for k, v in locals().items() if k not in ['seq_names', 'seqs', 'output']})
    image = create_dotplots(seq_names, seqs, config)
    image.save(output)

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
