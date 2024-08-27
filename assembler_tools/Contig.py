import json
import os
from collections import Counter

from .utils import human_bp, invariant_atgc_count


class Contig:
    importer: object
    fasta_file: str
    original_id: str
    original_contig_header: str
    sequence: str
    assembler: str = None
    topology: str = None
    location: str = None
    coverage: int = None
    additional_info: list[str] = []
    _len = None

    def __init__(
            self,
            importer,
            fasta_file: str,
            original_contig_header: str,
            assembler: str,
            sequence: str = None,
    ):
        self.importer = importer
        self.fasta_file = fasta_file
        self.original_contig_header = original_contig_header
        self.original_id = original_contig_header.split(' ', 1)[0].rsplit('|', 1)[-1]
        self.sequence = sequence
        if sequence is not None:
            atgc_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0} | Counter(sequence)
            assert set(atgc_count) == {'A', 'T', 'C', 'G'}, \
                f'Error in {self}: Invalid characters in sequence: atgc_count={self.atgc_count}'
            self.atgc_count = invariant_atgc_count(atgc_count)
        self.assembler = assembler

    @property
    def id(self):
        return f'{self.assembler}@{self.original_id}'

    def header(self, contig_name: str, plasmid_name: str = None) -> str:
        self.sanity_check()

        header = f'>{contig_name} [length={len(self)}]'
        if self.topology == 'circular':
            header += f' [topology=circular] [completeness=complete]'
        elif self.topology == 'linear':
            header += f' [topology=linear]'
        if self.location == 'chromosome':
            header += f' [location=chromosome]'
        elif self.location == 'plasmid':
            assert plasmid_name is not None, f'Error in {self}: no plasmid_name!'
            header += f' [location=plasmid]'
            header += f' [plasmid-name={plasmid_name}]'
        if self.coverage:
            header += f' [coverage={self.coverage}x]'
        if self.assembler:
            header += f' [assembler={self.assembler}]'
        if self.original_id:
            header += f' [old-id={self.original_id}]'

        for info in self.additional_info:
            header += f' {info}'

        return header

    def __repr__(self) -> str:
        return f'{self.fasta_file}:{self.original_id}'

    def __str__(self):
        return f'<Contig: {self.assembler}:{self.original_id} {self.len_human()} {self.topology}>'

    def __len__(self) -> int:
        if self._len:
            return self._len
        return len(self.sequence)

    def len_human(self) -> str:
        return human_bp(len(self))

    def sanity_check(self):
        if self.topology:
            assert self.topology in ['circular', 'linear'], f'Error in {self}: Invalid topology: {self.topology}'
        if self.location:
            assert self.location in ['chromosome', 'plasmid'], f'Error in {self}: Invalid location: {self.location}'

    @property
    def gc_abs(self) -> int:
        return self.atgc_count['G'] + self.atgc_count['C']

    @property
    def gc_rel(self) -> float:
        return self.gc_abs / len(self)

    @property
    def topology_badge(self):
        template = '<div class="badge rounded-pill bg-{color} me-1">{topology}</div>'
        if self.topology == 'circular':
            return template.format(topology='c', color='success')
        elif self.topology == 'linear':
            return template.format(topology='l', color='warning')
        elif self.topology == 'unknown':
            return template.format(topology='u', color='info')
        else:
            return template.format(topology='?', color='danger')

    @property
    def has_coverage(self) -> bool:
        return self.coverage is not None

    @property
    def coverage_badge(self):
        template = '<div class="badge rounded-pill bg-{color} me-1">{coverage}</div>'
        if not self.has_coverage:
            return template.format(coverage='no coverage!', color='info')
        try:
            coverage = round(float(self.coverage))
        except (ValueError, TypeError):
            return template.format(coverage=self.coverage, color='secondary')
        if coverage >= 50:
            return template.format(coverage=f'{coverage}x', color='success')
        elif coverage >= 30:
            return template.format(coverage=f'{coverage}x', color='warning')
        else:
            return template.format(coverage=f'{coverage}x', color='danger')

    @property
    def atgc_badge(self):
        template = '<div class="badge rounded-pill bg-{color} me-1">atgc</div>'
        # if any in self.atgc_count below 2 %: danger
        if any([count / len(self) < 0.02 for count in self.atgc_count.values()]):
            return template.format(color='danger')
        # if any in self.atgc_count below 5 %: warning
        elif any([count / len(self) < 0.05 for count in self.atgc_count.values()]):
            return template.format(color='warning')
        else:
            return ''




    def to_json(self, sequence: bool = False, contig_group: str = None, additional_data: dict = {}) -> dict:
        res = {
            'id': self.id,
            'original_id': self.original_id,
            'assembler': self.assembler,
            'len': len(self),
            'atgc_count': self.atgc_count,
            'gc_abs': self.gc_abs,
            'gc_rel': self.gc_rel,
            'coverage': self.coverage,
            'topology': self.topology,
            'location': self.location,
            'additional_info': self.additional_info,
            'test-header': self.header('test_scf0', plasmid_name='test-plasmid')
        }
        if sequence: res['sequence'] = self.sequence
        if contig_group: res['contig_group'] = contig_group
        res.update(additional_data)
        return res

    @classmethod
    def from_json(cls, data: str | dict):
        if type(data) is str:
            with open(data) as f:
                data = json.load(f)
        contig = cls(
            importer=None,
            fasta_file=None,
            original_contig_header=data['original_id'],
            sequence=None,
            assembler=data['assembler']
        )
        contig.original_id = data['original_id']
        contig._len = data.get('len', None)
        contig.atgc_count = data.get('atgc_count', None)
        contig.topology = data['topology']
        contig.location = data['location']
        contig.coverage = data['coverage']
        contig.additional_info = data['additional_info']
        contig.sequence = data.get('sequence', None)
        contig.contig_group = data.get('contig_group', None)
        return contig
