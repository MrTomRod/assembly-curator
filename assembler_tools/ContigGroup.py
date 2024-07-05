from .utils import human_bp
from .Contig import Contig


class ContigGroup:
    contigs: [Contig]

    def __init__(self, contigs: [Contig] = None):
        self.contigs = contigs if contigs else []

    def __len__(self):
        return sum([len(contig) for contig in self.contigs])

    def len_human(self) -> str:
        return human_bp(len(self))

    def __str__(self):
        return f'<ContigGroup: {self.assembler}:{self.id} {self.len_human()}>'

    def pprint(self):
        if len(self.contigs) == 1:
            print(self.contigs[0])
        else:
            print(self)
            for contig in self.contigs:
                print(f"  - {contig}")

    @property
    def id(self):
        return self.contigs[0].original_id

    @property
    def id_nice(self):
        n_contigs = len(self.contigs)
        first_contig_id = self.contigs[0].original_id
        if n_contigs == 1:
            return first_contig_id
        return f"{first_contig_id}+{n_contigs - 1}"

    @property
    def assembler(self):
        return self.contigs[0].assembler

    def sort(self):
        self.contigs = sorted(self.contigs, key=lambda contig: len(contig), reverse=True)

    def gc(self) -> int:
        return sum([contig.gc for contig in self.contigs])

    def gc_content(self) -> float:
        return self.gc() / len(self)
