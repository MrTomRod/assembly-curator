from .utils import human_bp
from .Contig import Contig


class ContigGroup:
    contigs: [Contig]
    cluster_id: int = None
    cluster_color: (float, float, float) = None
    cluster_color_rgb: str = None

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
        additional_contigs = '' if len(self.contigs) == 1 else f'+{len(self.contigs) - 1}'
        return f'{self.assembler}#{self.contigs[0].original_id}{additional_contigs}'

    @property
    def assembler(self):
        return self.contigs[0].assembler

    def sort(self):
        self.contigs = sorted(self.contigs, key=lambda contig: len(contig), reverse=True)

    @property
    def gc_abs(self) -> int:
        return sum([contig.gc_abs for contig in self.contigs])

    @property
    def gc_rel(self) -> float:
        return self.gc_abs / len(self)

    def set_cluster_color(self, r: float, g: float, b: float):
        self.cluster_color = (r, g, b)
        self.cluster_color_rgb = f'rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})'

    def topology_or_n_contigs(self, short: bool = False) -> str:
        if len(self.contigs) == 1:
            topology = self.contigs[0].topology
            return topology[0] if short else topology
        return f'n={len(self.contigs)}' if short else f'{len(self.contigs)} contigs'

    def encode_sequences(self) -> [bytes]:
        # Format used by skani
        return [contig.sequence.encode('ascii') for contig in self.contigs]

    def to_json(self, sequence: bool = False):
        info = {}
        if self.cluster_id:
            info['cluster_id'] = self.cluster_id
        if self.cluster_color:
            info['cluster_color'] = self.cluster_color
        if self.cluster_color_rgb:
            info['cluster_color_rgb'] = self.cluster_color_rgb
        res = {
            'id': self.id,
            'len': len(self),
            'gc_abs': self.gc_abs,
            'gc_rel': self.gc_rel,
            'assembler': self.assembler,
            'contigs': {contig.id: contig.to_json(sequence, self.id, info) for contig in self.contigs},
            'topology_or_n_contigs': self.topology_or_n_contigs()
        }
        res.update(info)
        return res
