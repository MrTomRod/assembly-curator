import os
import logging
import pandas as pd
from assembler_tools.utils import AssemblyFailedException
from assembler_tools.AssemblyImporter import AssemblyImporter
from assembler_tools.Assembly import Assembly
from assembler_tools.Contig import Contig


class LjaImporter(AssemblyImporter):
    assembler = 'lja'
    assembly_dir = 'lja'
    assembly = 'assembly.fasta'
    gfa = 'mdbg.gfa'
    plot = 'mdbg.gfa.svg'

    def load_assembly(self) -> Assembly:
        contigs = self.load_fasta(f"{self._assembly_dir_abs}/{self.assembly}")
        gfa, connections, circular = self.load_gfa(f"{self._assembly_dir_abs}/{self.gfa}")

        for contig in contigs:
            assert contig in connections, f'{contig=} not in {connections=}'

        self.declare_topology(contigs, circular)

        groups = self.create_groups(connections)
        assembly = self.create_assembly(groups, contigs)
        assembly.assembly_dir = self.assembly_dir
        assembly.plot = self.plot
        assembly.gfa = self.gfa

        return assembly
