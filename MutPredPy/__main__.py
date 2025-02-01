import pandas as pd
import argparse
import numpy as np
import os


from . import prep
from . import status
from . import merge
from . import remainder


def time_minimum(x):
    x = int(x)

    if x < 2:
        raise argparse.ArgumentTypeError("Minimum time is 2 hours")
    return x


def command_mutpred_prep(args):
    
    mut = prep.MutPredpy(
        input=args.input,
        working_dir=args.working_dir,
        time=args.time,
        dry_run=args.dry_run,
        canonical=args.canonical,
        database=args.database,
        all_possible=args.all_possible,
        fasta_path=args.fasta
    )
    
    mut.prepare_mutpred_input()


def command_status(args):
    
    summary = status.Status(
        job_dir=args.job_dir,
        summary=args.summary,
        all=args.all,
        show_incomplete=args.show_incomplete,
        logs=args.logs,
        script=args.script
    )

    summary.mutpred_summary()


def command_debug(args):
    mutpred_stat = status.Status(
        job_dir=args.job_dir,
        summary=args.summary,
        all=args.all,
        show_incomplete=args.show_incomplete,
        logs=args.logs
    )

    mutpred_stat.mutpred_debugging()



def command_merge(args):

    mutpred_merge = merge.Merge(
        output=args.output,
        job_dir=args.job_dir,
        dry_run=args.dry_run,
        mechanisms=args.mechanisms
    )
    mutpred_merge.merge()



def command_remainder(args):
    
    
    R = remainder.Remaining(
        working_dir=args.working_dir,
        time=args.time,
        dry_run=args.dry_run
    )
    
    R.split_and_build_lsf()



def build_parser():
    
    parser = argparse.ArgumentParser(description="MutPredPy functions")

    subparsers = parser.add_subparsers(
        title="commands",
        description="The following commands are available:",
        help='For additional help: "MutPredPy <COMMAND> -h"',
    )

    
    # ============== MutPredPy Prepare ==============
    parser_mutpredPrepare = subparsers.add_parser(
            "prepare", 
            help="Takes an input file and outputs faa files in preparation to run MutPred suite"
        )

    parser_mutpredPrepare.add_argument(
                '--time', type=time_minimum, nargs="?", default=24, 
                help="Target time in hours to run the jobs."
            )
    parser_mutpredPrepare.add_argument(
                '--input', type=str, nargs='?',
                help='Path to the input file with VEP or other variant annotations.'
            )
    parser_mutpredPrepare.add_argument(
                '--working_dir', type=str, nargs='?',
                help='Path to the directory where the outputs are written.'
            )
    parser_mutpredPrepare.add_argument(
                "--dry_run", action="store_true", 
                help='Calculate the resources that would be used with the current parameters, but write nothing to file.'
            )
    parser_mutpredPrepare.add_argument(
                "--canonical", action="store_true", 
                help='Only prepare canonical isoforms for mutpred scoring'
            )
    parser_mutpredPrepare.add_argument(
                "--all_possible", action="store_true", 
                help='Prepare to score all possible amino acid substitutions '
            )
    parser_mutpredPrepare.add_argument(
                "--database", type=str, nargs='?', default="None",
                help="Path to option config.yaml file for linking to mysql database. Include the name of the configuration in the config.yaml file after an @ symbol (ex. /path/to/file@Remote). If no config name included, program will default to 'Local'."
            )
    parser_mutpredPrepare.add_argument(
                "--fasta", type=str, nargs='?',
                help='Path to the combined fasta file from Ensembl.'
    )
    parser_mutpredPrepare.set_defaults(func=command_mutpred_prep)
    


    # ============== MutPredPy Status ==============
    parser_mutpredStatus = subparsers.add_parser(
            "status", 
            help="Check the status of a given job with the given parameters"
        )
    parser_mutpredStatus.add_argument(
                '--job_dir', type=str, nargs='?',
                help='Path to the directory containing the MutPred2 jobs.'
            )
    parser_mutpredStatus.add_argument(
                '--summary', type=str, nargs='?', required=False,
                help='(Optional) Path to a summary file that was previously generated by the status command. Will only use the summary file and will not compute the latest status of the jobs. Will ignore the job_dir variable.'
            )
    parser_mutpredStatus.add_argument(
                "--all", action="store_true", 
                help='Show the progress status of all files and jobs'
            )
    parser_mutpredStatus.add_argument('--show_incomplete', action="store_true", required=False,
                help='When listing the remaining jobs, show all incomplete jobs and not just jobs with zero output.'
            )
    parser_mutpredStatus.add_argument(
                '--logs', type=str, nargs='?', required=False, default="",
                help='Path to the directory containing the job logs.'
            )
    parser_mutpredStatus.add_argument(
                '--script', type=str, nargs='?', required=False, default=False,
                help='If path is included, MutPredPy will output new scripts with adjustments for memory and time errors. Path to one of the scripts that was run for this project.'
            )
    parser_mutpredStatus.set_defaults(func=command_status)



    # ============== MutPredPy Debug ==============
    parser_mutpredDebug = subparsers.add_parser(
            "debug", 
            help="Check the status of a given job with the given parameters"
        )

    parser_mutpredDebug.add_argument(
                '--job_dir', type=str, nargs='?',
                help='Path to the directory containing the MutPred2 jobs.'
            )
    parser_mutpredDebug.add_argument(
                '--summary', type=str, nargs='?', required=False,
                help='(Optional) Path to a summary file that was previously generated by the status command. Will only use the summary file and will not compute the latest status of the jobs. Will ignore the job_dir variable.'
            )
    parser_mutpredDebug.add_argument(
                "--all", action="store_true", 
                help='Show the progress status of all files and jobs'
            )
    parser_mutpredDebug.add_argument('--show_incomplete', action="store_true", required=False,
                help='When listing the remaining jobs, show all incomplete jobs and not just jobs with zero output.'
            )
    parser_mutpredDebug.add_argument(
                '--logs', type=str, nargs='?', required=False, default="",
                help='Path to the directory containing the job logs.'
            )

    parser_mutpredDebug.set_defaults(func=command_debug)



    # ============== MutPredPy Merge ==============
    parser_mutpredMerge = subparsers.add_parser(
            "merge", 
            help="Combine the output scored mutations from MutPred into a single file."
        )
    parser_mutpredMerge.add_argument(
                "--output", type=str, nargs='?',
                help='Path to the output file'
            )
    parser_mutpredMerge.add_argument(
                '--job_dir', type=str, nargs='?',
                help='Path to the directory containing the MutPred2 jobs.'
            )
    parser_mutpredMerge.add_argument(
                '--mechanisms', action="store_true", 
                help='Should the output merged results include the mechanisms and their posterior probabilities'
            )
    parser_mutpredMerge.add_argument(
                "--dry_run", action="store_true", 
                help='Run through the merging process without saving the output'
            )
    parser_mutpredMerge.set_defaults(func=command_merge)

    

    # ============== MutPredPy Remainder ==============
    parser_mutpredRemainder = subparsers.add_parser(
            "remainder", 
            help="Find the leftover unscored variants and create new jobs."
        )


    parser_mutpredRemainder.add_argument('--working_dir', type=str, nargs='?',
                        help='Path to the working directory of the ongoing mutpred2 project. A "jobs" folder and "scripts" folder should be in this working directory.')
    parser_mutpredRemainder.add_argument('--time', type=time_minimum, nargs="?", default=24,
                        help="Target time in hours to run the jobs")
    parser_mutpredRemainder.add_argument("--dry_run", action="store_true")

    parser_mutpredRemainder.set_defaults(func=command_remainder)

    return parser



def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError:
        parser.print_help()
        parser.exit()


if __name__ == "__main__":
    main()