import os.path
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import gfapy

from .utils import AssemblyFailedException
from .Assembly import Assembly
from .Contig import Contig
from .ContigGroup import ContigGroup


class AssemblyImporter(ABC):
    assembler: str = None
    assembly_dir: str = None  # The directory where the assembly is located, relative to sample_dir

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
        if not os.path.isfile(fasta):
            raise AssemblyFailedException(f'{self.name}: FASTA file {fasta} does not exist')
        with open(fasta) as f:
            data = f.read().strip()
        if not data:
            raise AssemblyFailedException(f'{self.name}: FASTA file {fasta} is empty')

        data = [e.strip() for e in data.split('>')]

        contigs = {}
        for entry in data:
            if entry == '': continue  # start or end of FASTA file
            contig_header, sequence = entry.split('\n', 1)
            sequence = sequence.replace('\n', '')
            assert set(sequence) == {'A', 'T', 'C', 'G'}

            contig = Contig(
                importer=self,
                fasta_file=fasta,
                original_contig_header=contig_header,
                sequence=sequence,
                assembler=self.assembler
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
        assembly = Assembly(assembler=self.assembler, assembly_dir=self.assembly_dir, sample_dir=self._sample_dir)
        for group in groups:
            contig_group = ContigGroup()
            for segment in group:
                if segment in contigs:
                    contig_group.contigs.append(contigs[segment])
            if contig_group.contigs:
                assembly.contig_groups.append(contig_group)

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
        gfa = gfapy.Gfa.from_file(gfa)

        connections = {}  # {segment:  set(segment)}
        circular = set()  # {segment}

        def connect(segments1, segments2):
            for segment1 in segments1:
                for segment2 in segments2:
                    connections.setdefault(segment1, set()).add(segment2)
                    connections.setdefault(segment2, set()).add(segment1)
                    if segment1 == segment2:
                        circular.add(segment1)

        for line in gfa.lines:
            if line.record_type in '#HAS':
                pass
            elif line.record_type in 'L':  # Link
                connect([line.from_name], [line.to_name])
            elif line.record_type in 'P':  # Path
                connect([line.name], [l.name for l in line.segment_names])
            else:
                raise AssertionError(f'Unknown record type: {line.record_type}')
        return gfa, connections, circular
