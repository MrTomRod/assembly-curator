import json
import os.path

from .utils import human_bp
from .ContigGroup import ContigGroup


class Assembly:
    assembler: str
    assembly_dir: str
    _sample_dir: str
    contig_groups: [ContigGroup]
    plot: str = None
    gfa: str = None
    busco: {} = None

    def __init__(self, assembler: str, assembly_dir: str, sample_dir: str, contig_groups: [ContigGroup] = None):
        self.assembler = assembler
        self.assembly_dir = assembly_dir
        self._sample_dir = sample_dir
        self.contig_groups = contig_groups if contig_groups else []

    def __len__(self):
        return sum([len(contig_group) for contig_group in self.contig_groups])

    def len_human(self) -> str:
        return human_bp(len(self))

    def __str__(self):
        return f"<Assembly: {self.assembler} {self.len_human()} {len(self.contig_groups)} contig groups>"

    def has_contig(self, contig):
        for group in self.contig_groups:
            if contig in group.contigs:
                return True
        return False

    def sort(self):
        for group in self.contig_groups:
            group.sort()
        self.contig_groups = sorted(self.contig_groups, key=lambda group: len(group), reverse=True)

    def load_busco(self, busco_file: str):
        with open(busco_file) as busco_file:
            self.busco = json.load(busco_file)

    def pprint(self):
        print(f'#### {self} ####')
        if self.plot:
            print(f'Plot: file:///{os.path.abspath(f'{self.assembly_dir}/{self.plot}')}')
        for group in self.contig_groups:
            group.pprint()
        print(f'#### --- end --- ####\n')

    def gc(self) -> int:
        return sum([contig.gc_abs for contig in self.contig_groups])

    def gc_content(self) -> float:
        return self.gc() / len(self)

    def plot_svg(self):
        if self.plot is None:
            return None
        try:
            with open(f'{self._sample_dir}/{self.assembly_dir}/{self.plot}') as f:
                return f.read()
        except Exception as e:
            return f'Error: {e}'

    def to_json(self, sequence: bool = False):
        return {
            'assembler': self.assembler,
            'assembly_dir': self.assembly_dir,
            'sample_dir': self._sample_dir,
            'len': len(self),
            'plot': self.plot,
            'gfa': self.gfa,
            'busco': self.busco,
            'contig_groups': {group.id: group.to_json(sequence) for group in self.contig_groups}
        }
