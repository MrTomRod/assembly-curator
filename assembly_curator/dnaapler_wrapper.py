import os
import shutil
from pathlib import Path

import click
from click.testing import CliRunner

# Create a dummy click context
def create_dummy_context():
    runner = CliRunner()
    with runner.isolated_filesystem():
        ctx = click.Context(click.Command('dummy'))
        ctx.params = {
            'input': 'dummy_input',
            'output': 'dummy_output',
            'threads': 4,
            'prefix': 'dummy_prefix',
            'evalue': 1e-10,
            'force': True,
            'autocomplete': 'none',
            'seed_value': 13,
            'ignore': None,
            'db': 'all',
            'custom_db': '',
        }
        return ctx

# Example usage
ctx = create_dummy_context()

from loguru import logger

from dnaapler import all
from dnaapler.utils.all import all_process_blast_output_and_reorient
from dnaapler.utils.bulk import bulk_process_blast_output_and_reorient, run_bulk_blast
from dnaapler.utils.cds_methods import (
    run_blast_based_method,
    run_largest,
    run_mystery,
    run_nearest,
)
from dnaapler.utils.constants import DNAAPLER_DB
from dnaapler.utils.external_tools import ExternalTool
from dnaapler.utils.util import (
    begin_dnaapler,
    check_duplicate_headers,
    end_dnaapler,
    get_version,
    print_citation,
    run_autocomplete,
)
from dnaapler.utils.validation import (
    check_evalue,
    instantiate_dirs,
    validate_choice_autocomplete,
    validate_choice_db,
    validate_choice_mode,
    validate_custom_db_fasta,
    validate_fasta,
    validate_fasta_all,
    validate_fasta_bulk,
    validate_ignore_file,
)


def all(
        input='data/15_N/lja/assembly.fasta',
        output='data/15_N/lja/dnaapler',
        threads=8,
        prefix='prefix',
        evalue=1e-10,
        force=True,
        autocomplete='none',
        seed_value=13,
        ignore='',  # Text file listing contigs (one per row) that are to be ignored
        db='all',
        custom_db='',
):
    """Reorients contigs to begin with any of dnaA, repA, terL or archaeal COG1474 Orc1/cdc6"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
        "--ignore": ignore,
        "--custom_db": custom_db,
        "--force": force,
        "--db": db,
    }

    # defines gene
    gene = "all"

    # other options

    if db == "dnaa":
        gene = "dnaA"
    elif db == "repa":
        gene = "repA"
    elif db == "terl":
        gene = "terL"
    elif db == "dnaa,repa":
        gene = "dnaA,repA"
    elif db == "dnaa,terl":
        gene = "dnaA,terL"
    elif db == "repa,terl":
        gene = "repA,terL"
    elif db == "cog1474":
        gene = "cog1474"

    # custom
    if custom_db != "":
        gene = "custom"

    # initial logging etc
    # start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta
    validate_fasta_all(input)
    check_duplicate_headers(input)

    # validate e value
    check_evalue(evalue)

    # create flag for ignore
    if ignore == "":
        ignore_flag = False
    else:
        ignore_flag = True
    # checks if the ignore file exists and contains text
    if ignore_flag == True:
        logger.info(f"You have specified contigs to ignore in {ignore}.")
        exists_contains_txt = validate_ignore_file(ignore)

    if gene == "custom":
        # validates custom fasta input for database
        validate_custom_db_fasta(Path(custom_db))

        # make db
        db_dir = os.path.join(output, "custom_db")
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        custom_db_fasta = os.path.join(db_dir, "custom_db.faa")
        shutil.copy2(custom_db, custom_db_fasta)

        logdir = Path(f"{output}/logs")

        # custom db
        # make custom db
        custom_database = os.path.join(db_dir, "custom_db")
        makeblastdb = ExternalTool(
            tool="makeblastdb",
            input=f"-in {custom_db_fasta}",
            output=f"-out {custom_database}",
            params="-dbtype prot ",
            logdir=logdir,
        )

        ExternalTool.run_tool(makeblastdb, ctx)
    else:
        custom_db = None

    # runs bulk BLAST
    run_bulk_blast(
        ctx, input, output, prefix, gene, evalue, threads, custom_db=custom_db
    )

    # rerorients blast
    blast_file = os.path.join(output, f"{prefix}_blast_output.txt")

    ### ignore
    # list is empty
    ignore_list = []
    if ignore_flag == True:
        if exists_contains_txt is False:
            logger.warning(f"{ignore} contains no text. No contigs will be ignored")
        else:
            # gets all contigs in the ignore
            # will split by space so short_contig only (to match BLAST)
            with open(ignore) as f:
                ignore_dict = {x.rstrip().split()[0] for x in f}
            ignore_list = list(ignore_dict)

    all_process_blast_output_and_reorient(
        input,
        blast_file,
        output,
        prefix,
        ignore_list,
        autocomplete,
        seed_value,
        custom_db=custom_db,
    )

    # end dnaapler
    # end_dnaapler(start_time)


all()
