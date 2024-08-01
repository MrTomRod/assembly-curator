from pywebio.output import put_text, put_markdown

from assembler_tools.ContigGroup import ContigGroup, Contig
from pywebio.input import input, FLOAT, SELECT, TEXT, select


def scf_names_dialact(sample_name, scf_id: int):
    return f"{sample_name.rsplit('-', 1)[0]}_scf{scf_id}"


def plasmid_names_dialact(sample_name, plasmid_id: int):
    return f"p{sample_name.rsplit('-', 2)[0]}_{plasmid_id}"


def determine_location(contig: Contig) -> None:
    if contig.location:
        return  # nothing to do

    # determine location based on length
    if contig.topology == 'circular':
        if 5_000 < len(contig) < 500_000:
            contig.location = 'plasmid'
            return
        elif len(contig) > 1_700_000:
            contig.location = 'chromosome'
            return

    contig.location = select(
        f'Please provide a location for {contig.id}',
        options=['chromosome', 'plasmid', 'unknown'], value='unknown'
    )


def determine_topology(contig: Contig) -> None:
    if contig.topology:
        return  # nothing to do

    # determine topology based on input
    contig.topology = select(
        f'Please provide a topology for {contig.id}',
        options=['circular', 'linear', 'unknown'], value='unknown'
    )


def create_headers(sample, contig_groups: [ContigGroup]) -> [str]:
    headers = {}
    contig_count, plasmid_count = 0, 0
    for cg_id, cg in enumerate(contig_groups):
        for contig in cg.contigs:
            contig_count += 1
            determine_location(contig)
            determine_topology(contig)
            if contig.location == 'plasmid':
                plasmid_count += 1
            headers[contig.id] = dict(
                header=create_header(
                    contig,
                    name=scf_names_dialact(sample, contig_count),
                    contig_group_id=cg_id,
                    plasmid_name=plasmid_names_dialact(sample, plasmid_count),
                    # name=f'{sample}_scf{contig_count}',
                    # plasmid_name=f'p{sample}_{contig_count}',
                ),
                contig=contig,
                contig_group=cg
            )

    markdown = f'## {sample}\n```text\n'
    for contig_id, data in headers.items():
        markdown += f"{data['header']}\n"
        markdown += f"{data['contig'].sequence[0:20]}...\n"
    markdown += '```'

    put_markdown(markdown)
    return headers


def create_header(contig: Contig, name: str, contig_group_id: int, plasmid_name: str = None) -> str:
    header = f'>{name} [length={len(contig)}]'
    if contig.topology == 'circular':
        header += f' [topology=circular] [completeness=complete]'
    elif contig.topology == 'linear':
        header += f' [topology=linear]'
    else:
        topology = input('Please provide a topology', type=SELECT, options=['circular', 'linear', 'unknown'])
        header += f' [topology={topology}]'
    if contig.location == 'chromosome':
        header += f' [location=chromosome]'
    elif contig.location == 'plasmid':
        if plasmid_name is None:
            plasmid_name = input(f'Please provide a plasmid name for {contig.id}', type=TEXT, value=name)
        header += f' [location=plasmid]'
        header += f' [plasmid-name={plasmid_name}]'
    if contig.coverage:
        header += f' [coverage={contig.coverage}x]'
    if contig.assembler:
        header += f' [assembler={contig.assembler}]'
    header += f' [contig-group={contig_group_id}]'
    if contig.original_id:
        header += f' [old-id={contig.original_id}]'

    for info in contig.additional_info:
        header += f' {info}'

    return header
