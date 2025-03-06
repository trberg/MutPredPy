"""
Catalog Module for MutPredPy

This module handles job cataloging, including checking directories, 
loading features, and collecting mechanisms from MutPred2 outputs.
"""

import os
import re
import importlib.resources as pkg_resources
import logging
import numpy as np
import pandas as pd
import scipy.io as sio

from .catalog_job import CatalogJob
from ..utils import utils as u

logger = logging.getLogger()


class Catalog:
    """Handles job cataloging for MutPred2 results."""

    def __init__(self, job_dir, mechanisms=False, features=False, dry_run=False):
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
        self.features = features

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

        # Feature columns
        self.feature_columns = self.load_feature_columns()

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

    def load_feature_columns(self):
        """
        Loads a list of feature columns names for feature files from mutpred2
        """
        if hasattr(self, "feature_columns"):
            return self.feature_columns

        feat_cols = pd.read_csv(
            pkg_resources.files("mutpredpy").joinpath(
                "pkg_resources", "mp2_features.csv"
            ),
            sep=",",
        )
        return feat_cols

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

        job_dirs = os.listdir(self.job_dir)
        # job_dirs = ["263"]

        number_of_jobs = len(job_dirs)
        cur_job = 0
        self.print_progress(cur_job, number_of_jobs)

        for job in job_dirs:

            job_path = os.path.join(self.job_dir, job)

            catalog_job = CatalogJob(job_path, self)

            catalog_job.process_job(self)

            cur_job += 1

            self.print_progress(cur_job, number_of_jobs)

        self.print_progress(cur_job, number_of_jobs, end="\n")


if __name__ == "__main__":

    pass
