import logging
import pandas as pd
from assembler_tools.utils import AssemblyFailedException
from assembler_tools.AssemblyImporter import AssemblyImporter
from assembler_tools.Assembly import Assembly
from assembler_tools.Contig import Contig


class FlyeImporter(AssemblyImporter):
    assembler = 'flye'
    assembly_dir = 'flye'
    assembly = 'assembly.fasta'
    gfa = 'assembly_graph.gfa'
    plot = 'assembly_graph.gfa.svg'

    def load_assembly(self) -> Assembly:
        contigs = self.load_fasta(f"{self._assembly_dir_abs}/{self.assembly}")
        gfa, connections, circular = self.load_gfa(f"{self._assembly_dir_abs}/{self.gfa}")

        for contig in contigs:
            assert contig in connections, f'{contig=} not in {connections=}'

        self.declare_topology(contigs, circular)
        self.load_assembly_info(contigs, f"{self._assembly_dir_abs}/assembly_info.txt")

        groups = self.create_groups(connections)
        assembly = self.create_assembly(groups, contigs)
        assembly.assembly_dir = self.assembly_dir
        assembly.plot = self.plot
        assembly.gfa = self.gfa

        return assembly

    def declare_topology(self, contigs: [Contig], circular):
        # Fix the topology declaration
        circular = {c.replace('edge', 'contig') for c in circular}
        super().declare_topology(contigs, circular)

    def load_assembly_info(self, contigs: [Contig], assembly_info):
        """
        Parse the Flye-specific assembly_info.txt file and:
            1) add coverage information to contigs
            2) sanity check: make sure gfa was parsed correctly
        """
        df = pd.read_csv(assembly_info, sep='\t', index_col=0)

        assert set(df.index) == set(contigs), \
            f'Mismatch in {assembly_info} and {contigs=}!\n{set(df.index)=}\n{set(contigs)=}'

        for contig_name, contig in contigs.items():
            assert contig_name in df.index, f'{contig_name=} not in {df.index}'
            info_length = df.loc[contig_name, 'length']
            if info_length != len(contig):
                logging.warning(f'Warning: {info_length=} != {len(contig)}')
            info_circular = df.loc[contig_name, 'circ.'] == 'Y'
            if (contig.topology == 'circular') != info_circular:
                logging.warning(f'Warning: {contig_name=} {contig.topology=} != {df.loc[contig_name, 'circ.']}')
            # Append coverage information
            contig.coverage = int(df.loc[contig_name, 'cov.'])
