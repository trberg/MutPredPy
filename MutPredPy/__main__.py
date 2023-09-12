import pandas as pd
import argparse
import numpy as np
import os

from . import ( prep,
                status,
                merge)


def time_minimum(x):
    x = int(x)

    if x < 5:
        raise argparse.ArgumentTypeError("Minimum time is 5 hours")
    return x


def command_mutpred_prep(args):
    
    mut = prep.MutPredpy(
        input=args.input,
        project=args.project,
        time=args.time,
        dry_run=args.dry_run,
        canonical=args.canonical
    )
    
    mut.prepare_mutpred_input()


def command_status(args):
    summary = status.Status(
        input=args.input,
        project=args.project
    )

    summary.mutpred_status()


def command_merge(args):
    
    merge.merge()



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
                '--input', type=str, nargs=1,
                help='The name of the input filename that is located in the data folder'
            )
    parser_mutpredPrepare.add_argument(
                '--project', type=str, nargs=1,
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
    parser_mutpredPrepare.set_defaults(func=command_mutpred_prep)
    


    # ============== MutPredPy Status ==============
    parser_mutpredStatus = subparsers.add_parser(
            "status", 
            help="Check the status of a given job with the given parameters"
        )

    parser_mutpredStatus.add_argument(
                '--input', type=str, nargs=1,
                help='The name of the input filename that is located in the data folder'
            )
    parser_mutpredStatus.add_argument(
                '--project', type=str, nargs=1,
                help='The name of the project for organization purposes'
            )
    parser_mutpredStatus.set_defaults(func=command_status)



    # ============== MutPredPy Merge ==============
    parser_mutpredMerge = subparsers.add_parser(
            "merge", 
            help="Combine the output scored mutations from MutPred into a single file."
        )
    
    parser_mutpredMerge.set_defaults(func=command_merge)

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