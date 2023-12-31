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
        project=args.project,
        time=args.time,
        dry_run=args.dry_run,
        canonical=args.canonical,
        database=args.database,
        assembly=args.assembly,
    )
    
    mut.prepare_mutpred_input()


def command_status(args):
    summary = status.Status(
        input=args.input,
        project=args.project,
        all=args.all
    )

    summary.mutpred_summary()


def command_debug(args):
    mutpred_stat = status.Status(
        input=args.input,
        project=args.project,
        all=args.all
    )

    mutpred_stat.mutupred_debugging()


def command_merge(args):

    merge.merge(dry_run=args.dry_run)


def command_remainder(args):
    
    
    R = remainder.Remaining(
        input=args.input,
        project=args.project,
        exclude=args.exclude,
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
                help="Target time in hours to run the jobs"
            )
    parser_mutpredPrepare.add_argument(
                '--input', type=str, nargs='?',
                help='The name of the input filename that is located in the data folder'
            )
    parser_mutpredPrepare.add_argument(
                '--project', type=str, nargs='?',
                help='The name of the project for organization purposes'
            )
    parser_mutpredPrepare.add_argument(
                "--dry_run", action="store_true", 
                help='Calculate the resources that would be used with the current parameters'
            )
    parser_mutpredPrepare.add_argument(
                "--canonical", action="store_true", 
                help='Only prepare canonical isoforms for mutpred scoring'
            )
    parser_mutpredPrepare.add_argument(
                "--database", type=str, nargs='?', default="None",
                help="Path to option config.yaml file for linking to mysql database. Include the name of the configuration in the config.yaml file after an @ symbol (ex. /path/to/file@Remote). If no config name included, program will default to 'Local'."
            )
    parser_mutpredPrepare.add_argument(
                "--assembly", type=str, nargs='?', default="hg38",
                help="Identify which genome assembly the variants are built from (hg19, hg38)"
    )
    parser_mutpredPrepare.set_defaults(func=command_mutpred_prep)
    


    # ============== MutPredPy Status ==============
    parser_mutpredStatus = subparsers.add_parser(
            "status", 
            help="Check the status of a given job with the given parameters"
        )

    parser_mutpredStatus.add_argument(
                '--input', type=str, nargs='?',
                help='The name of the input filename that is located in the data folder'
            )
    parser_mutpredStatus.add_argument(
                '--project', type=str, nargs='?',
                help='The name of the project for organization purposes'
            )
    parser_mutpredStatus.add_argument(
                "--all", action="store_true", 
                help='Show the progress status of all files and jobs'
            )
    parser_mutpredStatus.set_defaults(func=command_status)



    # ============== MutPredPy Status ==============
    parser_mutpredDebug = subparsers.add_parser(
            "debug", 
            help="Check the status of a given job with the given parameters"
        )

    parser_mutpredDebug.add_argument(
                '--input', type=str, nargs='?',
                help='The name of the input filename that is located in the data folder'
            )
    parser_mutpredDebug.add_argument(
                '--project', type=str, nargs='?',
                help='The name of the project for organization purposes'
            )
    parser_mutpredDebug.add_argument(
                "--all", action="store_true", 
                help='Show the progress status of all files and jobs'
            )
    parser_mutpredDebug.set_defaults(func=command_debug)



    # ============== MutPredPy Merge ==============
    parser_mutpredMerge = subparsers.add_parser(
            "merge", 
            help="Combine the output scored mutations from MutPred into a single file."
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

    parser_mutpredRemainder.add_argument('--input', type=str, nargs='?',
                        help='The name of the input filename that is located in the data folder.')
    parser_mutpredRemainder.add_argument('--project', type=str, nargs='?',
                        help='The name of the project for organization purposes')
    parser_mutpredRemainder.add_argument('--exclude', type=str, nargs='?', required=False, default='[]',
                        help='LSF sequence of jobs to exlude from the remainder function')
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