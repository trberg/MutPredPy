"""
Merge Module for MutPredPy

This module handles the merging of MutPred2 job outputs into a consolidated 
file, including mutation scores, mechanisms, and additional metadata.
"""

import os
import pandas as pd


class Merge:
    """
    Handles the merging of MutPred2 job outputs into a single consolidated file.

    This class reads individual job output files, processes mutation scores and
    associated data, and writes the merged results to an output file.

    Attributes:
        output (str): Path to the final merged output file.
        job_dir (str): Directory containing MutPred2 job outputs.
        mechanisms (bool): Whether to include mechanism-related data.
        dry_run (bool): If True, performs a dry run without writing output.
        scored (pd.DataFrame): DataFrame containing merged MutPred2 scores.
    """

    def __init__(self, output, job_dir, mechanisms, dry_run):

        self.output = output
        self.job_dir = job_dir
        self.mechanisms = mechanisms
        self.dry_run = dry_run
        self.scored = self.mutpred_scores()

    def write_output(self, data, output_file):
        """
        Writes the merged mutation data to a specified output file.

        Args:
            data (pd.DataFrame): DataFrame containing merged mutation scores and metadata.
            output_file (str): Path to the output file.

        Returns:
            None
        """
        data.to_csv(output_file, sep="\t", index=False)

    def mutpred_scores(self):
        """
        Collects and merges MutPred2 output scores from all job directories.

        Returns:
            pd.DataFrame: A DataFrame containing merged mutation scores, including
                        mechanisms and additional metadata if enabled.
        """
        if self.mechanisms:
            included_columns = [
                "ID",
                "Substitution",
                "MutPred2 score",
                "Molecular mechanisms with Pr >= 0.01 and P < 1.00",
                "Motif information",
                "Remarks",
            ]

        else:
            included_columns = ["ID", "Substitution", "MutPred2 score"]

        print("MERGE")
        scores = pd.concat(
            [
                pd.read_csv(f"{self.job_dir}/{m}/output.txt", usecols=included_columns)
                for m in os.listdir(self.job_dir)
                if os.path.exists(f"{self.job_dir}/{m}/output.txt")
            ]
        )

        scores = scores.drop_duplicates()

        return scores

    def merge(self):
        """
        Merges MutPred2 job outputs into a single file and writes the results.

        Returns:
            None
        """
        self.write_output(self.scored, self.output)
