"""
Catalog Module for MutPredPy

This module handles job cataloging, including checking directories,
loading features, and collecting mechanisms from MutPred2 outputs.
"""

import os
import re
import sys
import importlib.resources as pkg_resources
import logging
import numpy as np
import pandas as pd
from mpi4py import MPI

from .catalog_job import CatalogJob
from ..utils import utils as u

# Initialize Logger
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

        os.environ["KMP_WARNINGS"] = "0"  # Disables OpenMP warnings
        os.environ["OMP_DISPLAY_ENV"] = "FALSE"  # Hides OpenMP environment messages

        u.set_logger_level(log_level=20)

        self.__job_dir, self.__valid_jobs = self.check_jobs_directory(job_dir)
        self.catalog_location = u.catalog_directory()
        print(f"Writing to {self.catalog_location}")

        self.dry_run = dry_run
        self.mechanisms = mechanisms
        self.features = features

        if self.mechanisms:
            self.load_mechanism_properties()

    def load_mechanism_properties(self):
        """Loads mechanism-related properties to avoid reloading in child classes"""
        neutral_properties = self.load_neutral_distributions()

        self.neutral_property_cols = [
            col.replace("-", "_") for col in neutral_properties.columns
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
        null_distributions_array = neutral_properties.to_numpy()

        self.null_gain_distributions_array = null_distributions_array[
            :, neutral_property_gain_cols
        ]

        self.null_loss_distributions_array = null_distributions_array[
            :, neutral_property_loss_cols
        ]

        self.n = null_distributions_array.shape[0]

        self.null_distributions = neutral_properties

        # Feature columns
        self.feature_columns = self.load_feature_columns()

    def load_neutral_distributions(self):
        """
        Loads the null distributions for the MutPred2 features to calculated p-values
        """
        dists = pd.read_csv(
            pkg_resources.files("mutpredpy").joinpath(
                "pkg_resources", "neutral_property_distributions.csv"
            ),
            sep="\t",
            dtype=str,
        ).astype(np.float64)

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

    def get_job_dir(self):
        """
        Retrieves the absolute path of the MutPred2 job directory.

        Returns:
            str: The absolute path of the job directory.
        """
        return self.__job_dir

    def get_valid_jobs(self):
        """
        Retrieves a list of valid jobs.

        Returns:
            list: A list of valid job folders available in the job directory.
        """
        return self.__valid_jobs

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
        jobs = []
        job_pattern = re.compile(
            r"^(job(s)?_\d+|\d+)$"
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
                    jobs.append(job)
            else:
                # If it's not a valid job folder, flag it
                if os.path.isdir(job_path):
                    discrepancies.append(f"Unexpected folder: {job}\n")
                else:
                    discrepancies.append(f"Unexpected file: {job}\n")

        if len(discrepancies) > 0:
            logger.info("Issues found in job directory %s.", job_dir)
            for disc in discrepancies:
                logger.info(disc)

        if len(jobs) == 0:
            raise ValueError(f"No valid job folders found in {job_path}")

        return job_dir, jobs

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
            f"[{'=' * cur_progress}{' ' * remaining_progress}] {round((cur_percentage)*100, 1)}%",
        )
        sys.stdout.flush()

    def catalog_jobs(self):
        """
        Processes all jobs in the directory, extracting relevant features,
        and stores the information in the catalog.

        Returns:
            pd.DataFrame: DataFrame containing cataloged job information.
        """

        # Initialize MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()  # Process ID
        size = comm.Get_size()  # Total number of processes

        if rank == 0:

            job_dirs = self.get_valid_jobs()
            # job_dirs = ["7716"]
            job_chunks = [job_dirs[i::size] for i in range(size)]

            number_of_processes = len(job_chunks)
            total_jobs = sum([len(j) for j in job_chunks])

            logger.info(
                "%s processes processing %d jobs", number_of_processes, total_jobs
            )

        else:
            job_chunks = None

        job_chunk = comm.scatter(job_chunks, root=0)

        if job_chunk is None:
            job_chunk = []

        for job in job_chunk:
            job_path = os.path.join(self.get_job_dir(), job)

            catalog_job = CatalogJob(job_path, self)

            catalog_job.process_job()


if __name__ == "__main__":

    pass
