"""
Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Trycycler

This file is part of Trycycler. Trycycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Trycycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Trycycler.
If not, see <http://www.gnu.org/licenses/>.
"""

import logging
import collections
from PIL import Image, ImageDraw, ImageFont, ImageColor

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


def load_contig_sequences(*fasta_paths_and_contigs: tuple[str]) -> dict[str, str]:
    """
    Load contig sequences from the given FASTA files and contig IDs.

    :param fasta_paths_and_contigs: Paths to FASTA files followed by contig IDs
        E.g. ('path/to/1.fasta:contig3 path/to/2.fasta:c5')
    :return: A dictionary of contig names to sequences
    """
    seqs = {}
    for fasta_path_and_contig in fasta_paths_and_contigs:
        fasta_path, contig_id = fasta_path_and_contig.rsplit(':', 1)
        seqs[contig_id] = extract_contig(fasta_path, contig_id)
    return seqs


log = logging.info
section_header = logging.info
explanation = logging.info


class DotplotConfig:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def dotplot(
        *fasta_paths_and_contigs: tuple[str], output: str = 'dotplots.png',
        res=800, kmer=10,
        initial_top_left_gap=0.1, border_gap=0.015, between_seq_gap=0.0125, outline_width=0.0015,
        text_gap=0.005, max_font_size=0.025,
        background_colour='white', self_vs_self_colour='lightgrey', self_vs_other_colour='whitesmoke',
        text_colour='black', forward_strand_dot_colour='mediumblue', reverse_strand_dot_colour='firebrick'
):
    """
    Draws all pairwise dotplots for a set of contigs in a cluster.

    This tool is heavily based on the dotplot tool from Trycycler:
    https://github.com/rrwick/Trycycler/blob/main/trycycler/dotplot.py

    :param fasta_paths_and_contigs: Paths to FASTA files followed by contig IDs, separated by a colon.
        E.g. ('path/to/1.fasta:contig3', 'path/to/2.fasta:c5')
    :param output: The output file path for the generated dotplots image.
    :param res: Size (in pixels) of each dot plot image (between 500 and 10000).
    :param kmer: K-mer size to use in dot plots.
    :param initial_top_left_gap: Initial gap size in the top left corner of the image.
    :param border_gap: Gap size between the border and the sequences.
    :param between_seq_gap: Gap size between sequences.
    :param outline_width: Width of the outline around each sequence box.
    :param text_gap: Gap size between the text labels and the sequence boxes.
    :param max_font_size: Maximum font size for the text labels.
    :param background_colour: Background color of the image.
    :param self_vs_self_colour: Color for self vs self comparisons.
    :param self_vs_other_colour: Color for self vs other comparisons.
    :param text_colour: Color of the text labels.
    :param forward_strand_dot_colour: Color of the dots for forward strand matches.
    :param reverse_strand_dot_colour: Color of the dots for reverse strand matches.
    """
    config = DotplotConfig(**{k: v for k, v in locals().items() if k not in ['fasta_paths_and_contigs', 'output']})
    welcome_message()
    assert 500 <= config.res <= 10000, 'Error: res size must be between 500 and 10000'
    seqs = load_contig_sequences(*fasta_paths_and_contigs)
    seq_names = sorted(seqs.keys())
    image = create_dotplots(seq_names, seqs, config)

    image.save(output)
    log(f'Saving dotplots to: {output}')
    finished_message()


def welcome_message():
    section_header('Starting Trycycler dotplot')
    explanation('Trycycler dotplot is a tool for drawing all pairwise dotplots for a set of '
                'contigs in a cluster. This step is optional but may help with interpretation '
                'of difficult clusters when performing reconciliation.')


def finished_message():
    section_header('Finished!')
    explanation('You can now examine the dotplots.png file to help understand how the contigs in '
                'this cluster relate to each other.')


def create_dotplots(seq_names, seqs, config: DotplotConfig):
    section_header('Creating dotplots')
    explanation('Each pairwise combination of sequences is now compared to generate a dotplot, '
                'all of which will be combined into a single image. Same-strand k-mer matches are '
                'drawn in blue, and opposite-strand k-mer matches are drawn in red.')
    seq_names = sorted(seq_names)

    # Set some sizes in terms of pixels
    top_left_gap = int(round(config.initial_top_left_gap * config.res))
    border_gap = max(2, int(round(config.border_gap * config.res)))
    between_seq_gap = max(2, int(round(config.between_seq_gap * config.res)))
    text_gap = max(1, int(round(config.text_gap * config.res)))
    outline_width = max(1, int(round(config.outline_width * config.res)))
    max_font_size = max(1, int(round(config.max_font_size * config.res)))

    # We create an initial image to test the label sizes.
    start_positions, end_positions, bp_per_pixel = \
        get_positions(seq_names, seqs, config.kmer, top_left_gap, border_gap, between_seq_gap, config.res)
    image = Image.new('RGB', (config.res, config.res), config.background_colour)
    min_font_size, max_text_height = \
        draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width,
                    max_font_size, config.background_colour, config.text_colour)

    # Now that we know the values for min_font_size and max_text_height, we start over, this time
    # limiting the font size to the minimum (so all text is the same size) and readjusting the
    # top-left gap (so it isn't bigger than necessary).
    top_left_gap = max_text_height + border_gap
    start_positions, end_positions, bp_per_pixel = \
        get_positions(seq_names, seqs, config.kmer, top_left_gap, border_gap, between_seq_gap, config.res)
    image = Image.new('RGB', (config.res, config.res), config.background_colour)
    draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width, config.self_vs_self_colour,
                        config.self_vs_other_colour)
    draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width,
                min_font_size, config.background_colour, config.text_colour)

    for name_a in seq_names:
        seq_a = seqs[name_a]
        for name_b in seq_names:
            seq_b = seqs[name_b]
            log(f'  {name_a} vs {name_b}')
            draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, config.kmer,
                      config.forward_strand_dot_colour, config.reverse_strand_dot_colour)

    # The boxes are drawn once more, this time with no fill. This is to overwrite any dots which
    # leaked into the outline, which would look messy.
    draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width, config.self_vs_self_colour,
                        config.self_vs_other_colour, fill=False)

    return image


def get_positions(seq_names, seqs, kmer_size, top_left_gap, bottom_right_gap,
                  between_seq_gap, res):
    """
    This function returns the image coordinates that start/end each sequence. Since the dot plot is
    symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    """
    seq_lengths = {n: len(seqs[n]) - kmer_size + 1 for n in seq_names}
    all_gaps = top_left_gap + bottom_right_gap + between_seq_gap * (len(seq_names) - 1)
    pixels_for_sequence = res - all_gaps
    total_seq_length = sum(seq_lengths[n] for n in seq_names)
    bp_per_pixel = total_seq_length / pixels_for_sequence

    start_positions, end_positions = {}, {}
    current_pos = top_left_gap
    for name in seq_names:
        start_positions[name] = current_pos
        rect_size = int(round(seq_lengths[name] / bp_per_pixel))
        current_pos += rect_size
        end_positions[name] = current_pos
        current_pos += between_seq_gap

    return start_positions, end_positions, bp_per_pixel


def draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width,
                        self_vs_self_colour, self_vs_other_colour, fill=True):
    """
    This function draws the box for each of the dot plots in the full image.
    """
    draw = ImageDraw.Draw(image)
    for name_a in seq_names:
        start_a = start_positions[name_a] - outline_width
        end_a = end_positions[name_a] + outline_width
        for name_b in seq_names:
            start_b = start_positions[name_b] - outline_width
            end_b = end_positions[name_b] + outline_width
            if fill:
                if name_a == name_b:
                    fill_colour = self_vs_self_colour
                else:
                    fill_colour = self_vs_other_colour
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               fill=fill_colour, outline='black', width=outline_width)
            else:
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               outline='black', width=outline_width)


def draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width, font_size,
                background_colour, text_colour):
    draw = ImageDraw.Draw(image)

    min_pos = min(p for p in start_positions.values())
    font_sizes, text_heights = [], []
    for name in seq_names:
        font, text_width, text_height, font_size = \
            get_font(draw, name, font_size, start_positions[name], end_positions[name])
        font_sizes.append(font_size)
        text_heights.append(text_height)

        # Horizontal labels on the top side.
        pos = min_pos - text_height - outline_width - text_gap
        draw.text((start_positions[name], pos), name, font=font, fill=text_colour)

        # Vertical labels on the left side.
        image_2 = Image.new('RGBA', (text_width, text_height), background_colour)
        draw_2 = ImageDraw.Draw(image_2)
        draw_2.text((0, 0), text=name, font=font, fill=text_colour)
        image_2 = image_2.rotate(90, expand=1)
        sx, sy = image_2.size
        image.paste(image_2, (pos, end_positions[name] - sy, pos + sx, end_positions[name]),
                    image_2)
    return min(font_sizes), max(text_heights)


def draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, kmer_size,
              forward_strand_dot_colour, reverse_strand_dot_colour):
    pixels = image.load()
    a_start_pos = start_positions[name_a]
    b_start_pos = start_positions[name_b]

    forward_colour = ImageColor.getrgb(forward_strand_dot_colour)
    reverse_colour = ImageColor.getrgb(reverse_strand_dot_colour)

    a_forward_kmers, a_reverse_kmers = get_all_kmer_positions(kmer_size, seq_a)

    for j in range(len(seq_b) - kmer_size + 1):
        j_pixel = int(round(j / bp_per_pixel)) + b_start_pos
        k = seq_b[j:j + kmer_size]
        if k in a_reverse_kmers:
            for i in a_reverse_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                draw_dot(pixels, i_pixel, j_pixel, reverse_colour)

        if k in a_forward_kmers:
            for i in a_forward_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                draw_dot(pixels, i_pixel, j_pixel, forward_colour)


def draw_dot(pixels, i, j, colour):
    try:
        pixels[i, j] = colour
        pixels[i + 1, j] = colour
        pixels[i - 1, j] = colour
        pixels[i, j + 1] = colour
        pixels[i, j - 1] = colour
    except IndexError:
        pass


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


def get_font(draw, label, font_size, start_position, end_position):
    font, is_default_font = load_font(font_size)

    # If we had to resort to the default font, then we can't size it.
    if is_default_font:
        left, top, right, bottom = font.getbbox(label)
        text_width, text_height = right - left, bottom - top
        return font, text_width, text_height, font_size

    # If we have a size-able font, then we adjust the size down until the text fits in the
    # available space.
    available_width = end_position - start_position
    left, top, right, bottom = font.getbbox(label)
    text_width, text_height = right - left, bottom - top
    while text_width > available_width:
        font_size -= 1
        font, _ = load_font(font_size)
        left, top, right, bottom = font.getbbox(label)
        text_width, text_height = right - left, bottom - top
    return font, text_width, text_height, font_size


def load_font(font_size):
    """
    This function loads a font, but it has to try a few different ones because different platforms
    have different fonts.
    """
    fonts = [
        'DejaVuSans.ttf', 'OpenSans-Regular.ttf', 'Arial.ttf', 'LiberationSans-Regular.ttf',
        'NimbusSans-Regular.otf', 'Verdana.ttf', 'Lato-Regular.ttf', 'FreeSans.ttf'
    ]

    for font in fonts:
        try:
            return ImageFont.truetype(font, font_size), False
        except OSError:
            continue

    return ImageFont.load_default(), True


if __name__ == '__main__':
    from fire import Fire

    Fire(dotplot)
