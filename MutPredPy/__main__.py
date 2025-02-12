import typer
from typing import Optional
from .prep import prep
from .status import status
from .merge import merge
#from .remainder import remainder

app = typer.Typer()

def time_minimum(x: int):
    """Ensure minimum time requirement is met."""
    if x < 2:
        raise typer.BadParameter("Minimum time is 2 hours")
    return x

## PREPARE
@app.command()
def prepare(
    input: str = typer.Option(help="Path to the input file"),
    working_dir: str = typer.Option(help="Path to output directory"),
    time: int = typer.Option(24, callback=time_minimum, help="Target time in hours"),
    dry_run: bool = typer.Option(False, help="Perform a dry run without writing files"),
    canonical: bool = typer.Option(False, help="Only prepare canonical isoforms"),
    all_possible: bool = typer.Option(False, help="Prepare all possible amino acid substitutions")
):
    """Prepare MutPred input files."""
    prep.Prepare(
        input=input,
        working_dir=working_dir,
        time=time,
        dry_run=dry_run,
        canonical=canonical,
        all_possible=all_possible,
    ).prepare_mutpred_input()

## STATUS
@app.command()
def status(
    job_dir: str = typer.Option(..., help="Path to the MutPred2 jobs directory"),
    summary: Optional[str] = typer.Option(None, help="Path to a pre-generated summary file"),
    all: bool = typer.Option(False, help="Show all job statuses"),
    show_incomplete: bool = typer.Option(False, help="Show all incomplete jobs"),
    logs: str = typer.Option("", help="Path to job logs"),
    script: Optional[str] = typer.Option(None, help="Path to the script for error adjustment"),
):
    """Check the status of a given job."""
    status.Status(
        job_dir=job_dir,
        summary=summary,
        all=all,
        show_incomplete=show_incomplete,
        logs=logs,
        script=script,
    ).mutpred_summary()

## DEBUG
@app.command()
def debug(
    job_dir: str = typer.Option(..., help="Path to the MutPred2 jobs directory"),
    summary: Optional[str] = typer.Option(None, help="Path to a pre-generated summary file"),
    all: bool = typer.Option(False, help="Show all job statuses"),
    show_incomplete: bool = typer.Option(False, help="Show all incomplete jobs"),
    logs: str = typer.Option("", help="Path to job logs"),
):
    """Debug a given MutPred job."""
    status.Status(
        job_dir=job_dir,
        summary=summary,
        all=all,
        show_incomplete=show_incomplete,
        logs=logs,
    ).mutpred_debugging()

## MERGE
@app.command()
def merge(
    output: str = typer.Option(..., help="Path to the output file"),
    job_dir: str = typer.Option(..., help="Path to MutPred2 jobs directory"),
    mechanisms: bool = typer.Option(False, help="Include mechanisms and probabilities in output"),
    dry_run: bool = typer.Option(False, help="Run merging process without saving output"),
):
    """Merge MutPred output files into a single file."""
    merge.Merge(
        output=output,
        job_dir=job_dir,
        dry_run=dry_run,
        mechanisms=mechanisms,
    ).merge()

#@app.command()
#def remainder(
#    working_dir: str = typer.Argument(..., help="Path to the working directory of MutPred2 project"),
#    time: int = typer.Option(24, callback=time_minimum, help="Target time in hours"),
#    dry_run: bool = typer.Option(False, help="Perform a dry run"),
#):
#    """Find the leftover unscored variants and create new jobs."""
#    remainder.Remaining(
#        working_dir=working_dir,
#        time=time,
#        dry_run=dry_run,
#    ).split_and_build_lsf()

#if __name__ == "__main__":
def main():
    app()
