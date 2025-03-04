"""
MutPredPy CLI

This module provides a command-line interface (CLI) for MutPredPy using Typer.
It allows users to:
- Prepare input data for MutPred2 (`prepare` command)
- Check job statuses (`check_status` command)
- Merge MutPred2 output files (`merge_results` command)
- Catalog MutPred2 results (`catalog_results` command)

Usage:
    python -m MutPredPy [COMMAND] [OPTIONS]

Author: Timothy Bergquist
"""

from typing import Optional
import typer

# Local module imports
from .prep import prep
from .status import status
from .merge import merge
from .catalog import Catalog

app = typer.Typer()


def time_minimum(x: int):
    """Ensure minimum time requirement is met."""
    if x < 2:
        raise typer.BadParameter("Minimum time is 2 hours")
    return x


## PREPARE
@app.command()
def prepare(  # pylint: disable=R0913,R0917
    input_file: str = typer.Option(help="Path to the input file"),
    working_dir: str = typer.Option(help="Path to output directory"),
    time: int = typer.Option(24, callback=time_minimum, help="Target time in hours"),
    canonical: bool = typer.Option(False, help="Only prepare canonical isoforms"),
    all_possible: bool = typer.Option(
        False, help="Prepare all possible amino acid substitutions"
    ),
    truncate: int = typer.Option(
        200000,
        help="Break up sequences into 100 AA windows on sequences \
            of length greater the input number.",
    ),
    verbose: bool = typer.Option(
        False,
        help="Print additional debugging and status information during the preparation",
    ),
    users: int = typer.Option(
        1,
        help="Designates how many users will be running the \
            MutPred2 scripts to better parallelize the jobs.",
    ),
    fasta: str = typer.Option(
        "",
        help="Pathway to a FASTA file from which to map the \
            proteins/transcripts from the input file",
    ),
    dry_run: bool = typer.Option(False, help="Perform a dry run without writing files"),
):
    """Prepare MutPred input files."""
    prep.Prepare(
        input_path=input_file,
        working_dir=working_dir,
        time=time,
        dry_run=dry_run,
        canonical=canonical,
        all_possible=all_possible,
        truncate=truncate,
        verbose=verbose,
        users=users,
        fasta_path=fasta,
    ).prepare_mutpred_input()


@app.command()
def check_status(  # pylint: disable=R0913,R0917
    job_dir: str = typer.Option(..., help="Path to the MutPred2 jobs directory"),
    summary: Optional[str] = typer.Option(
        None, help="Path to a pre-generated summary file"
    ),
    show_all: bool = typer.Option(False, help="Show all job statuses"),
    show_incomplete: bool = typer.Option(False, help="Show all incomplete jobs"),
    logs: str = typer.Option("", help="Path to job logs"),
    script: Optional[str] = typer.Option(
        None, help="Path to the script for error adjustment"
    ),
):
    """Check the status of a given job."""
    status.Status(
        job_dir=job_dir,
        summary=summary,
        show_all=show_all,
        show_incomplete=show_incomplete,
        logs=logs,
        script=script,
    ).mutpred_summary()


@app.command()
def merge_results(
    output: str = typer.Option(..., help="Path to the output file"),
    job_dir: str = typer.Option(..., help="Path to MutPred2 jobs directory"),
    mechanisms: bool = typer.Option(
        False, help="Include mechanisms and probabilities in output"
    ),
    dry_run: bool = typer.Option(
        False, help="Run merging process without saving output"
    ),
):
    """Merge MutPred output files into a single file."""
    merge.Merge(
        output=output,
        job_dir=job_dir,
        dry_run=dry_run,
        mechanisms=mechanisms,
    ).merge()


@app.command()
def catalog_results(
    job_dir: str = typer.Option(..., help="Path to MutPred2 jobs directory"),
    mechanisms: bool = typer.Option(
        False, help="Include mechanisms in the cataloging process"
    ),
    dry_run: bool = typer.Option(
        False, help="Run merging process without saving output"
    ),
):
    """Catalog all mutpred scores from a given job directory"""
    Catalog(job_dir=job_dir, mechanisms=mechanisms, dry_run=dry_run).catalog_jobs()


def main():
    """Run the Typer CLI application."""
    app()


if __name__ == "__main__":
    main()
