import logging
import os.path
from glob import glob
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from .utils import AssemblyFailedException, run_command
from .Assembly import Assembly
from .Contig import Contig
from .ContigGroup import ContigGroup


class AssemblyImporter(ABC):
    assembler: str = None
    assembly_dir: str = None  # The directory where the assembly is located, relative to sample_dir
    assembly: str = None  # The location of the assembly, relative to assembly_dir

    _sample_dir: str = None  # The directory where the sample is located
    _assembly_dir_abs: str = None  # The absolute path to the assembly directory

    def __init__(self, sample_dir: str):
        assert self.assembler is not None, f'{self.__class__.name} must define assembler'
        self._sample_dir = sample_dir
        self._assembly_dir_abs = os.path.join(sample_dir, self.assembly_dir)
        if not os.path.isdir(self._assembly_dir_abs):
            raise AssemblyFailedException(f'Folder {self._assembly_dir_abs} does not exist!', 'warning')

    @abstractmethod
    def load_assembly(self) -> Assembly:
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def load_fasta(self, fasta: str) -> {str: Contig}:
        matches = glob(fasta)
        if len(matches) != 1:
            raise AssemblyFailedException(
                f'{self.name}: Expected exactly one match for {fasta=}, but found {len(matches)}')
        fasta = matches[0]

        with open(fasta) as f:
            data = f.read().strip()
        if not data:
            raise AssemblyFailedException(f'{self.name}: FASTA file {fasta} is empty')

        contigs = {}
        for entry in data.split('>'):
            entry = entry.strip()
            if entry == '': continue  # start or end of FASTA file
            contig_header, sequence = entry.split('\n', 1)
            sequence = sequence.replace('\n', '')
            assert set(sequence) == {'A', 'T', 'C', 'G'}

            contig = Contig(
                importer=self,
                fasta_file=fasta,
                original_contig_header=contig_header,
                assembler=self.assembler,
                sequence=sequence,
            )
            contigs[contig.original_id] = contig

        return contigs

    def create_groups(self, connections) -> [[str]]:
        groups = [[segment] for segment in connections.keys()]
        for segment, connected_segments in connections.items():
            for connected_segment in connected_segments:
                # find and merge the groups that contain the connected segments
                group1 = None
                group2 = None
                for group in groups:
                    if segment in group:
                        group1 = group
                    if connected_segment in group:
                        group2 = group
                    if group1 and group2:
                        break
                if group1 is None:
                    raise AssertionError(f'{self.name}: {segment=} not in any group')
                if group2 is None:
                    raise AssertionError(f'{self.name}: {connected_segment=} not in any group')
                if group1 == group2:
                    continue
                group1.extend(group2)
                groups.remove(group2)

        return groups

    def create_assembly(self, groups: [[str]], contigs: {str: Contig}) -> Assembly:
        assembly = Assembly(assembler=self.assembler, assembly_dir=self.assembly_dir, assembly=self.assembly,
                            sample_dir=self._sample_dir)
        for group in groups:
            contig_group = ContigGroup()
            for segment in group:
                if segment in contigs:
                    contig_group.contigs.append(contigs.pop(segment))
                else:
                    logging.info(f'{self.name}: Segment {segment} not found in contigs')
            if contig_group.contigs:
                assembly.contig_groups.append(contig_group)
            else:
                logging.info(f'{self.name}: Empty contig group: {group=}')

        for name, contig in contigs.items():
            logging.info(f'{self.name}: Contig {name} not in any group')
            assembly.contig_groups.append(ContigGroup([contig]))

        # Sanity check: ensure no contigs were lost
        for name, contig in contigs.items():
            assert assembly.has_contig(contig), f'Contig {name} was lost! {contig=}'

        # Sort contig_groups by size
        assembly.sort()

        return assembly

    def declare_topology(self, contigs, circular):
        for contig in contigs.values():
            contig.topology = 'circular' if contig.original_id in circular else 'linear'

    def load_gfa(self, gfa: str):
        connections = {}  # {segment:  set(segment)}
        circular = set()  # {segment}

        def connect(segments1, segments2):
            for segment1 in segments1:
                for segment2 in segments2:
                    connections.setdefault(segment1, set()).add(segment2)
                    connections.setdefault(segment2, set()).add(segment1)
                    if segment1 == segment2:
                        circular.add(segment1)

        with open(gfa) as file:
            for line in file:
                parts = line.strip().split('\t')
                record_type = parts[0]

                if record_type == 'L':  # Link
                    from_name = parts[1]
                    to_name = parts[3]
                    connect([from_name], [to_name])
                elif record_type == 'P':  # Path
                    path_name = parts[1]
                    segment_names = [seg[:-1] for seg in parts[2].split(',')]
                    connect([path_name], segment_names)
                elif record_type in '#HAS':
                    continue
                else:
                    logging.warning(f'Unknown gfa record type: {record_type}. '
                                    f'Offending line (first 200 chars):\n{line[:200]}')

        return gfa, connections, circular

    def gfa_to_svg(self, gfa: str, overwrite: bool = True, params: [str] = ['--labels']):
        gfa_dirname = os.path.dirname(gfa)
        gfa_basename = os.path.basename(gfa)
        svg_basename = f'{gfa_basename}.svg'
        svg_path = os.path.join(gfa_dirname, svg_basename)

        if os.path.isfile(svg_path):
            if overwrite:
                logging.info(f'Overwriting {svg_basename}')
                os.remove(svg_path)
            else:
                logging.info(f'Skipping {svg_basename} as it already exists')
                return

        cmd = self._gfa_to_svg_cmd(gfa_basename, svg_basename, params)
        logging.info(f'Running: {cmd}')
        return_code = run_command(cmd, cwd=gfa_dirname)
        assert return_code == 0 and os.path.isfile(svg_path), f'Failed to create {svg_basename}'

    def _gfa_to_svg_cmd(self, gfa_basename: str, svg_basename: str, params: str = ['--labels']):
        return (['gfaviz-mrtomrod', '--no-gui', '--render'] + params + ['--output', svg_basename, gfa_basename])

    # def _gfa_to_svg_cmd(self, gfa_basename: str, svg_basename: str, params: str = '--labels'):
    #     return f'gfaviz-mrtomrod --no-gui --render {params} --output "{svg_basename}" "{gfa_basename}"'
    #
    # def _gfa_to_svg_wrap(self, cmd: str):
    #     return (f'podman run --rm '
    #             f'-v .:/data:Z '
    #             f'troder/gfaviz:latest '
    #             f'{cmd}')
