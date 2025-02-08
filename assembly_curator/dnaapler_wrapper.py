import os
import click
import logging
import subprocess
from tempfile import TemporaryDirectory
import pandas as pd
from Bio.Seq import reverse_complement

CPU_COUNT = os.cpu_count()

MESSAGE_DICT = {
    'Contig_already_reoriented': 'already-reoriented',
    'No_MMseqs2_hits': 'failed-no-hits',
}


def read_fasta(fasta_file) -> dict:
    # returns {header: {header: header, seq: sequence}}
    with open(fasta_file) as f:
        sequences = [s for s in f.read().lstrip().split('>') if s]
    sequences = [s.split('\n', 1) for s in sequences]
    sequences = [(h.split(' ', 1)[0], dict(header=h, seq=s)) for h, s in sequences]
    contig_ids = set([s[0] for s in sequences])
    assert len(contig_ids) == len(sequences), f"Non-unique contig_ids in {fasta_file}: {contig_ids=}"
    return dict(sequences)


def write_fasta(fasta: dict, output_file: str, replace_newlines: bool = True):
    with open(output_file, 'w') as f:
        for scf_id in fasta:
            header = fasta[scf_id]['header'].strip()
            sequence = fasta[scf_id]['seq'].strip()
            if replace_newlines:
                sequence = sequence.replace('\n', '')
            f.write(f">{header}\n{sequence}\n")


def read_output_table(output_table: str) -> pd.DataFrame:
    return pd.read_csv(output_table, sep='\t', index_col=0, dtype=str)


def run_dnaapler(
        input_fasta,
        workdir: str = None,
        n_cpu: int = None,
        arguments: [str] = ['-t', str(CPU_COUNT), '-f']
) -> (dict, pd.DataFrame):
    """
    Run dnaapler on the input fasta file and return the reoriented fasta and the reorientation summary table.§
    :param input_fasta: Path to the input fasta file
    :param workdir: Path to the output directory
    :param n_cpu: Number of CPUs to use
    :param arguments: Additional arguments to pass to dnaapler
    :return: (fasta as dict, output_table as pd.DataFrame)
    """
    if n_cpu is None:
        n_cpu = os.cpu_count()
    assert type(n_cpu) is int, f"n_cpu must be an integer, not {type(n_cpu)}"

    if workdir:
        os.makedirs(workdir, exist_ok=True)
    else:
        tempdir = TemporaryDirectory()
        workdir = tempdir.name

    input_fasta = os.path.abspath(input_fasta)
    input_fasta_temp = os.path.join(workdir, 'input.fasta')
    os.symlink(input_fasta, input_fasta_temp)

    arguments = [
        'dnaapler', 'all',
        '-i', 'input.fasta',
        '-o', 'out',
        *arguments
    ]

    logging.info(f"Running dnaapler in {workdir=}")
    logging.info(f"Running dnaapler with {arguments=}")
    logging.info(f"Running dnaapler with arguments={' '.join(arguments)}")

    res = subprocess.run(arguments, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error_msg = f"Error running dnaapler: {res.returncode}\n{res.stdout.decode()}\n{res.stderr.decode()}"
    assert res.returncode == 0, error_msg
    output_fasta = os.path.join(workdir, 'out', 'dnaapler_reoriented.fasta')
    assert os.path.isfile(output_fasta), error_msg
    output_table = os.path.join(workdir, 'out', 'dnaapler_all_reorientation_summary.tsv')
    assert os.path.isfile(output_table), error_msg

    fasta = read_fasta(output_fasta)
    output_table = read_output_table(output_table)

    logging.info(f"Finished dnaapler. Table:\n\n{output_table.to_string()}\n")

    if workdir is None:
        tempdir.cleanup()

    return fasta, output_table


def merge_fastas(
        fasta_orig: dict,
        fasta_dnaapler: dict,
        output_table: pd.DataFrame,
        default_topology: str = 'raise'
) -> dict:
    # make sure both dictionaries have the same keys
    assert list(fasta_orig.keys()) == list(fasta_dnaapler.keys()), (
        f"Keys do not match!\n"
        f"{list(fasta_orig.keys())} != {list(fasta_dnaapler.keys())}")

    assert default_topology in ['raise', 'linear', 'circular'], f"Unknown default_topology: {default_topology}"

    for scf_id in fasta_orig:
        row = output_table.loc[fasta_orig[scf_id]['header']]
        logging.info(f'{scf_id}: Dnaapler result: {row.Gene_Reoriented}')

        if not row.Start.isdigit():
            assert row.Gene_Reoriented in MESSAGE_DICT, f"Unknown dnaapler result: {row.Gene_Reoriented}\n{row=}"
            msg = MESSAGE_DICT[row.Gene_Reoriented]
            fasta_orig[scf_id]['header'] += f' [dnaapler={msg}]'
            continue

        compressed_header = fasta_orig[scf_id]['header'].replace(' ', '')

        def process_linear():
            fasta_orig[scf_id]['seq'] = reverse_complement(fasta_orig[scf_id]['seq'])
            fasta_orig[scf_id]['header'] += f' [dnaapler=reverse-complement] [dnaapler-gene={row.Gene_Reoriented}]'

        def process_circular():
            fasta_orig[scf_id]['seq'] = fasta_dnaapler[scf_id]['seq']
            fasta_orig[scf_id]['header'] += f' [dnaapler=rotated] [dnaapler-gene={row.Gene_Reoriented}]'

        if '[topology=linear]' in compressed_header:
            process_linear()
        elif '[topology=circular]' in compressed_header:
            process_circular()
        elif default_topology == 'linear':
            logging.warning(f'Topology not specified in header, falling back to default: {default_topology}')
            process_linear()
        elif default_topology == 'circular':
            logging.warning(f'Topology not specified in header, falling back to default: {default_topology}')
            process_circular()
        else:
            logging.error(f"Unknown topology in >{fasta_orig[scf_id]['header']}")
            raise ValueError(f"Unknown topology in {fasta_orig[scf_id]['header']}")
    return fasta_orig


def run_all(
        input_fasta,
        output_fasta,
        replace_newlines: bool = True,
        default_topology: str = 'raise',
        dnaapler_arguments: [str] = ['-t', str(CPU_COUNT), '-f'],
        n_cpu=None,
        workdir=None
):
    fasta_orig = read_fasta(input_fasta)
    fasta_dnaapler, output_table = run_dnaapler(input_fasta, workdir, n_cpu, dnaapler_arguments)
    fasta_result = merge_fastas(fasta_orig, fasta_dnaapler, output_table, default_topology)
    write_fasta(fasta_result, output_fasta, replace_newlines)   


@click.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.argument('output_fasta', type=click.Path())
@click.option('--replace-newlines/--no-replace-newlines', default=True)
@click.option('--default-topology', type=click.Choice(['raise', 'linear', 'circular']),
              default='raise', help='Default topology')
@click.option('--dnaapler-arguments', type=str, default=f'-t {CPU_COUNT} -f',
              help='Additional arguments to pass to dnaapler')
@click.option('--n-cpu', type=int, default=None, help='Number of CPUs to use')
@click.option('--workdir', type=click.Path(), default=None, help='Path to the output directory')
def cli(input_fasta, output_fasta, replace_newlines, default_topology, workdir, dnaapler_arguments, n_cpu):
    """
    Run dnaapler on the input FASTA file and save the reoriented FASTA to the output file.
    """
    dnaapler_arguments = dnaapler_arguments.split()

    logging.basicConfig(level=logging.INFO)
    run_all(input_fasta, output_fasta, replace_newlines, default_topology, dnaapler_arguments, n_cpu, workdir)


if __name__ == '__main__':
    cli()
