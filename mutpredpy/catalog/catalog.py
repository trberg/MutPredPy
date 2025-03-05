"""
Catalog Module for MutPredPy

This module handles job cataloging, including checking directories, 
loading features, and collecting mechanisms from MutPred2 outputs.
"""

import os
import re
import importlib.resources as pkg_resources
import logging
import yaml
import numpy as np
import pandas as pd
import scipy.io as sio

from .catalog_job import CatalogJob
from ..utils import utils as u
from ..fasta import Protein


class Catalog:
    """Handles job cataloging for MutPred2 results."""

    def __init__(self, job_dir, mechanisms=False, dry_run=False):
        """
        Initializes the Catalog class.

        Args:
            job_dir (str): Path to job directory.
            mechanisms (bool): Whether to include mechanisms.
            dry_run (bool): Run without saving changes.
        """

        u.set_logger_level(log_level=20)

        self.job_dir = self.check_jobs_directory(job_dir)
        self.catalog_location = u.catalog_directory()

        self.dry_run = dry_run
        self.mechanisms = mechanisms

        if self.mechanisms:
            self.load_mechanism_properties()

    def load_mechanism_properties(self):
        """Loads mechanism-related properties to avoid reloading in child classes"""
        neutral_properties = self.load_neutral_distributions()

        self.neutral_property_cols = [
            neuts[0].replace("-", "_")
            for neuts in neutral_properties.get("neut_properties").flatten()
        ]
        self.property_names = [
            p.replace("_gain", "")
            for p in self.neutral_property_cols
            if p.endswith("_gain") or p == "Stability"
        ]
        neutral_property_gain_cols = np.array(
            [
                self.neutral_property_cols.index(col)
                for col in self.neutral_property_cols
                if col.endswith("_gain") or col == "Stability"
            ]
        )
        neutral_property_loss_cols = np.array(
            [
                self.neutral_property_cols.index(col)
                for col in self.neutral_property_cols
                if col.endswith("_loss") or col == "Stability"
            ]
        )

        # Load null distributions
        self.null_distributions_array = neutral_properties.get("neut_ntrans_feats")
        self.null_gain_distributions_array = self.null_distributions_array[
            :, neutral_property_gain_cols
        ]
        self.null_loss_distributions_array = self.null_distributions_array[
            :, neutral_property_loss_cols
        ]

        self.n = self.null_distributions_array.shape[0]
        self.null_distributions = pd.DataFrame(
            data=self.null_distributions_array, columns=self.neutral_property_cols
        )

    def load_neutral_distributions(self):
        """
        Loads the null distributions for the MutPred2 features to calculated p-values
        """
        dists = sio.loadmat(
            pkg_resources.files("mutpredpy").joinpath(
                "pkg_resources", "pu_features_null_distributions.mat"
            )
        )
        return dists

    def check_jobs_directory(self, job_dir):
        """
        Validates the given job directory and its contents.

        Args:
            job_dir (str): Path to job directory.

        Returns:
            str: Validated job directory.

        Raises:
            FileNotFoundError: If the job directory does not exist.
        """

        discrepancies = []
        job_pattern = re.compile(
            r"^(job_\d+|\d+)$"
        )  # Pattern for "job_<number>" or just a number.

        # Check if the job directory exists
        if not os.path.exists(job_dir):
            raise FileNotFoundError(f"Job directory {job_dir} not found")

        # Iterate through items in the job directory
        for job in os.listdir(job_dir):
            job_path = os.path.join(job_dir, job)

            # Check if it's a valid job folder
            if os.path.isdir(job_path) and job_pattern.match(job):

                # Check required files inside the job folder
                required_files = ["input.faa", "output.txt"]
                missing_files = [
                    file
                    for file in required_files
                    if not os.path.isfile(os.path.join(job_path, file))
                ]

                if missing_files:
                    discrepancies.append(
                        f"Missing files in job '{job}': {', '.join(missing_files)}\n"
                    )
            else:
                # If it's not a valid job folder, flag it
                if os.path.isdir(job_path):
                    discrepancies.append(f"Unexpected folder: {job}\n")
                else:
                    discrepancies.append(f"Unexpected file: {job}\n")

        if len(discrepancies) > 0:
            raise ValueError(
                f"Issues found in job directory {job_path}.\n{''.join(discrepancies)}"
            )

        return job_dir

    def write_mutation_to_catalog(self, row):
        """
        Writes mutation data to the catalog directory in YAML format.

        Args:
            row (pd.Series): A Pandas Series containing mutation details including
                            substitution, MutPred2 score, and mechanisms.

        Returns:
            None
        """
        mutation = row["Substitution"]
        pos = re.search(r"\d+", mutation).group()
        _, alt = mutation.split(pos)

        if self.mechanisms:
            mutpred_data = {
                "Substitution": row["Substitution"],
                "MutPred2 Score": row["MutPred2 score"],
                "Mechanisms": row["Mechanisms"],
            }
        else:
            mutpred_data = {
                "Substitution": row["Substitution"],
                "MutPred2 Score": row["MutPred2 score"],
            }

        mutation_path = os.path.join(
            self.catalog_location, "scores", row["sequence_hash"], pos, alt
        )
        if self.dry_run:
            logging.info("Would've created catalog directory %s", mutation_path)

        else:
            os.makedirs(mutation_path, exist_ok=True)

        mutpred_output_yaml_path = os.path.join(mutation_path, "mutpred2_output.yaml")
        if self.dry_run:
            logging.info("Would've created the file %s", mutpred_output_yaml_path)

        else:
            with open(mutpred_output_yaml_path, "wt", encoding="utf-8") as yaml_file:
                yaml.dump(
                    mutpred_data, yaml_file, default_flow_style=False, sort_keys=False
                )

    def print_progress(self, cur_job, number_of_jobs, end="\r"):
        """
        Displays a progress bar for tracking job completion.

        Args:
            cur_job (int): The current job index.
            number_of_jobs (int): The total number of jobs.
            end (str, optional): The end character for the print statement. Defaults to "\r".

        Returns:
            None
        """
        cur_percentage = cur_job / number_of_jobs
        cur_progress = int(50 * cur_percentage)
        remaining_progress = 50 - cur_progress
        print(
            f" [{'=' * cur_progress}{' ' * remaining_progress}] {round((cur_percentage)*100, 1)}%",
            end=end,
        )

    def catalog_jobs(self):
        """
        Processes all jobs in the directory, extracting relevant features,
        and stores the information in the catalog.

        Returns:
            pd.DataFrame: DataFrame containing cataloged job information.
        """
        catalog_index = []

        job_dirs = os.listdir(self.job_dir)
        job_dirs = ["263"]

        number_of_jobs = len(job_dirs)
        cur_job = 0
        self.print_progress(cur_job, number_of_jobs)

        for job in job_dirs:

            job_path = os.path.join(self.job_dir, job)
            catalog_job = CatalogJob(job_path, self)
            job_info = catalog_job.process_job(self)
            print(job_info)
            exit()
            catalog_index.append(job_info)

            # output_path = os.path.join(self.job_dir, job, "output.txt")
            # input_path = os.path.join(self.job_dir, job, "input.faa")

            # output = pd.read_csv(output_path)
            # input_faa = fasta.read_mutpred_input_fasta(input_path)
            # input_faa["Substitution"] = input_faa["Substitution"].apply(
            #     lambda x: x.split(",")
            # )
            # input_faa = input_faa.explode("Substitution")

            # sequenced_data = output.merge(
            #     input_faa, on=["ID", "Substitution"], how="left"
            # )

            # if (
            #     self.mechanisms
            #     and "Molecular mechanisms with Pr >= 0.01 and P < 1.00"
            #     in output.columns
            # ):
            #     keep_cols = [
            #         "ID",
            #         "Substitution",
            #         "MutPred2 score",
            #         "Molecular mechanisms with Pr >= 0.01 and P < 1.00",
            #         "sequence",
            #     ]
            # else:
            #     keep_cols = ["ID", "Substitution", "MutPred2 score", "sequence"]

            # sequenced_data = sequenced_data.filter(keep_cols)

            # job_path = os.path.join(self.job_dir, job)

            # cur_job = Catalog_Job(job_path=job_path)

            # job_information = cur_job.get_catalog_information(self)
            # print(job_information)

            sequenced_data["sequence_hash"] = sequenced_data["sequence"].apply(
                u.get_seq_hash
            )

            sequenced_data.apply(self.write_mutation_to_catalog, axis=1)

            cur_index = Protein.split_mutpred_output_ids(
                sequenced_data[["ID", "sequence_hash"]].drop_duplicates()
            )

            catalog_index.append(cur_index)

            cur_job += 1
            self.print_progress(cur_job, number_of_jobs)

        self.print_progress(cur_job, number_of_jobs, end="\n")

        catalog_index = pd.concat(catalog_index)

        return catalog_index.drop_duplicates()


if __name__ == "__main__":

    pass
