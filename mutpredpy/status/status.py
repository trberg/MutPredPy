"""
Status Module for MutPredPy

This module provides functionality to check the status of ongoing or completed 
MutPred2 jobs, including monitoring job progress, retrieving logs, and handling 
error reporting.
"""

import re
import os
import argparse
from datetime import datetime

import pandas as pd

from ..fasta import fasta
from ..computing import lsf


class Status:
    """
    Manages the tracking and reporting of MutPred2 job statuses.

    This class provides methods to retrieve job directories, logs, summaries,
    and error reports, helping users monitor and troubleshoot MutPred2 runs.

    Attributes:
        job_dir (str): Path to the MutPred2 jobs directory.
        show_all (bool): Whether to show all job statuses.
        show_incomplete (bool): Whether to show only incomplete jobs.
        log_dir (str): Path to job logs directory.
        write_new_scripts (bool): Whether to generate new job scripts for failed jobs.
        summary (str or bool): Path to a pre-generated summary file or False if not used.
    """

    def __init__(
        self, job_dir, show_all, show_incomplete, summary="", logs="", script=False
    ):

        self.__job_dir = job_dir
        self.__show_all = show_all  # pylint: disable=unused-private-member
        self.__show_incomplete = show_incomplete
        self.__log_dir = logs
        self.__write_new_scripts = bool(script)
        if self.__write_new_scripts:
            self.__template_script = script

        if summary == "":
            self.__summary = False
        else:
            self.__summary = summary

    def write_script_status(self):
        """
        Determines whether new job scripts should be generated for failed MutPred2 jobs.

        Returns:
            bool: True if new job scripts should be written, otherwise False.
        """
        return self.__write_new_scripts

    def get_script(self):
        """
        Retrieves the absolute path of the script used for job execution.

        Returns:
            str: Absolute path of the script.
        """
        return os.path.abspath(self.__template_script)

    def get_job_dir(self):
        """
        Retrieves the absolute path of the MutPred2 job directory.

        Returns:
            str: The absolute path of the job directory.
        """
        return os.path.abspath(self.__job_dir.rstrip("/"))

    def get_log_dir(self):
        """
        Retrieves the absolute path of the log directory.

        Returns:
            str: The absolute path of the log directory.
        """
        if self.__log_dir == "":
            return self.__log_dir
        return os.path.abspath(self.__log_dir.rstrip("/"))

    def get_summary(self):
        """
        Retrieves the absolute path of the job summary file if provided.

        Returns:
            str or bool: Absolute path of the summary file if available, otherwise False.
        """
        if not self.__summary:
            return False
        return os.path.abspath(self.__summary)

    def get_mutpred_output_file(self, jobid):
        """
        Retrieves the output file name for a given MutPred2 job.

        Args:
            jobid (str): The job identifier.

        Returns:
            str: The name of the output file.

        Raises:
            Exception: If no output file or multiple output files are found.
        """

        job_dir = f"{self.get_job_dir()}/{jobid}"

        o_reg = re.compile(r"output.txt")

        out_files = [o for o in os.listdir(job_dir) if o_reg.match(o)]

        if len(out_files) == 1:
            out_file = out_files[0]
        elif len(out_files) > 1:
            raise ValueError(f"More than 1 output file for job {jobid}")
        elif len(out_files) == 0:
            raise FileNotFoundError(f"No output file for job {jobid}")
        else:
            raise FileNotFoundError(
                f"Error with reading output file for job {jobid}. File may not exist"
            )

        return out_file

    def get_mutpred_output_file_path(self, index):
        """
        Retrieves the absolute path of the MutPred2 output file for a given job.

        Args:
            index (str): The job identifier.

        Returns:
            str: Absolute path to the output file.
        """
        directory = self.get_job_dir()
        file = self.get_mutpred_output_file(index)
        return f"{directory}/{file}"

    def get_mutpred_input_file(self, index):
        """
        Retrieves the input FASTA file name for a given MutPred2 job.

        Args:
            index (str): The job identifier.

        Returns:
            str: The name of the input FASTA file.
        """
        file = self.retrieve_faas(jobid=index)

        return file

    def get_mutpred_input_file_path(self, index):
        """
        Retrieves the absolute path of the MutPred2 input FASTA file for a given job.

        Args:
            project (str): The project name or identifier.
            index (str): The job identifier.

        Returns:
            str: Absolute path to the input FASTA file.
        """
        directory = self.get_job_dir()
        file = self.get_mutpred_input_file(index)

        return f"{directory}/{file}"

    def read_mutpred_output(self, file, jobid):
        """
        Reads and processes the MutPred2 output file for a given job.

        Args:
            file (str): Path to the output file.
            jobid (str): The job identifier.

        Returns:
            pd.DataFrame: A DataFrame containing mutation scores grouped by ID.
        """
        if os.path.isfile(file) and os.path.getsize(file) > 0:

            scores = pd.read_csv(file)

            scores = scores.drop_duplicates()

            scores = pd.DataFrame(
                scores.groupby("ID")["Substitution"].agg(list)
            ).reset_index()

            scores["num_mutations"] = scores["Substitution"].apply(len).astype(int)

            scores["index"] = jobid

            scores = scores[["ID", "num_mutations", "index"]]

        else:
            scores = pd.DataFrame(columns=["ID", "num_mutations", "index"])

        return scores

    def read_err_log(self, log):
        """
        Reads an error log file and extracts the relevant error message.

        Args:
            log (str): Name of the error log file.

        Returns:
            str: Extracted error message or an empty string if no error is found.
        """
        error = ""

        with open(f"{self.get_log_dir()}/{log}", "r", encoding="utf-8") as l:

            for line in l:
                if len(line) > 1 and ("TERM" in line):
                    return line.strip("\n")
                if len(line) > 1 and (
                    "MATLAB is exiting because of fatal error" in line
                ):
                    next_line = next(l)
                    if "Killed" in next_line or "Bus" in next_line:
                        return line.strip("\n") + ": " + "Memory Error"
                    if "Bus" in next_line:
                        return line.strip("\n") + ": " + "Memory Error"
                    return line.strip("\n") + ": " + next_line.strip("\n")
        return error

    def has_err_log(self, log):
        """
        Checks if an error log file exists and contains relevant errors.

        Args:
            log (str): Name of the error log file.

        Returns:
            pd.Series: A Series containing:
                - bool: True if an error log exists and contains errors, otherwise False.
                - str: Extracted error message or an empty string.
        """
        if self.get_log_dir() == "":
            error = ""
            return pd.Series([False, error])

        cur_log = f"{self.get_log_dir()}/{log}"

        if os.path.exists(cur_log) and os.path.getsize(cur_log) > 0:
            error = self.read_err_log(log)
            return pd.Series([True, error])
        error = ""
        return pd.Series([False, error])

    def get_job_index(self, log):
        """
        Extracts the job index from a log file.

        Args:
            log (str): Path to the job log file.

        Returns:
            pd.Series: A Series containing:
                - str: Job identifier.
                - str: Job index.
        """
        with open(log, "r", encoding="utf-8") as log_file:
            for line in log_file:

                if line.startswith("Subject"):
                    job, index = line.split(":")[1].strip(" Job ").strip("]").split("[")
                    return pd.Series([job, index])
                    # Subject: Job 137347991[10]: <phase3[1-4000]> in cluster <chimera> Exited
        return pd.Series([None, None])

    def retrieve_logs(self):
        """
        Retrieves and processes all error logs associated with MutPred2 jobs.

        Returns:
            pd.DataFrame: A DataFrame containing job indices, error statuses, and error messages.
        """
        print(self.get_log_dir())
        if not os.path.isdir(self.get_log_dir()):

            log_files = pd.DataFrame({"logs": []})
            log_files["job"] = ""
            log_files["index"] = ""
            log_files["hasError"] = ""
            log_files["Error"] = ""
            log_files = log_files[["index", "hasError", "Error"]]

        else:

            log_files = pd.DataFrame(
                {
                    "logs": [
                        l for l in os.listdir(f"{self.get_log_dir()}/") if "out" in l
                    ]
                }
            )

            log_files[["job", "index"]] = log_files["logs"].apply(
                lambda x: self.get_job_index(f"{self.get_log_dir()}/{x}")
            )

            latest_err_logs = pd.DataFrame(
                log_files.groupby(["index"])["job"].max()
            ).reset_index()

            log_files = log_files.merge(
                latest_err_logs, on=["index", "job"], how="inner"
            )
            # print (log_files)

            log_files[["hasError", "Error"]] = log_files.apply(
                lambda row: self.has_err_log(row["logs"]), axis=1
            )
            log_files = log_files[["index", "hasError", "Error"]]

            log_files["Error"].fillna("No Errors", inplace=True)

        return log_files

    def retrieve_faas(self, jobid):
        """
        Retrieves and parses the input FASTA file for a given MutPred2 job.

        Args:
            jobid (str): The job identifier.

        Returns:
            pd.DataFrame: A DataFrame containing parsed FASTA sequence data for the job.

        Raises:
            Exception: If no input FASTA file is found or multiple files exist.
        """
        faa_dir = f"{self.get_job_dir()}/{jobid}"

        f_reg = re.compile(r".*.faa")

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f)]

        if len(faa_files) == 1:
            faa_file = faa_files[0]
        elif len(faa_files) > 1:
            raise ValueError(f"More than 1 input fasta file in job {jobid}")
        elif len(faa_files) == 0:
            raise FileNotFoundError(f"No input fasta file found for job {jobid}")
        else:
            raise FileNotFoundError(f"Error finding input fasta file for job {jobid}")

        job_faa = fasta.read_mutpred_input_fasta(f"{faa_dir}/{faa_file}")

        job_faa["index"] = jobid

        return job_faa

    def retrieve_outputs(self, jobid):
        """
        Retrieves the output data for a given MutPred2 job.

        Args:
            jobid (str): The job identifier.

        Returns:
            pd.DataFrame: A DataFrame containing the parsed output data from the job.

        Raises:
            Exception: If no output file is found, multiple output files exist,
                    or there is an error in reading the output file.
        """
        job_dir = f"{self.get_job_dir()}/{jobid}"

        # o_reg = re.compile(f".*.missense_output_{jobid}.txt$")
        o_reg = re.compile(r"output.txt$")

        out_files = [o for o in os.listdir(job_dir) if o_reg.match(o)]
        # if jobid == 71

        if len(out_files) == 1:
            out_file = out_files[0]
        elif len(out_files) == 0:
            out_file = "no_out_file.txt"
        elif len(out_files) > 1:
            raise ValueError(f"More than 1 MutPred2 output file for job {jobid}")
        else:
            raise FileNotFoundError(
                f"Error in reading output file for job {jobid}. File may be \
                    missing or inaccessible."
            )

        return self.read_mutpred_output(f"{job_dir}/{out_file}", jobid)

    def mutpred_status(self):
        """
        Retrieves and processes the status of all MutPred2 jobs in the job directory.

        Returns:
            pd.DataFrame: A DataFrame containing job statuses, including the number of
                        mutations scored, completeness status, and job index.
        """
        all_statuses = []

        for job in os.listdir(f"{self.get_job_dir()}"):

            if os.path.isdir(f"{self.get_job_dir()}/{job}"):

                if not os.path.exists(f"{self.get_job_dir()}/{job}/completed.txt"):

                    faa = self.retrieve_faas(job)
                    # print (faa)

                    scores = self.retrieve_outputs(job)
                    # print (scores)

                    # exit()

                    current_status = faa.merge(
                        scores,
                        on=["ID", "index"],
                        how="outer",
                        suffixes=("_faa", "_scored"),
                    ).fillna(0)
                    current_status["complete"] = current_status.apply(
                        lambda x: x["num_mutations_faa"] == x["num_mutations_scored"],
                        axis=1,
                    )

                    if sum(current_status["num_mutations_faa"]) == sum(
                        current_status["num_mutations_scored"]
                    ):
                        self.mark_as_completed(job)

                    all_statuses.append(current_status)

                else:

                    faa = self.retrieve_faas(job)

                    faa.rename(
                        columns={"num_mutations": "num_mutations_faa"}, inplace=True
                    )
                    faa["num_mutations_scored"] = faa["num_mutations_faa"]
                    # faa["percent"] = 100.0
                    faa["complete"] = True

                    all_statuses.append(faa)

        all_statuses = pd.concat(all_statuses)

        return all_statuses

    def mutpred_logs(self, current_status):
        """
        Retrieves and processes error logs associated with MutPred2 jobs.

        Args:
            status (pd.DataFrame): DataFrame containing job status information.

        Returns:
            pd.DataFrame: A DataFrame with job indices, error statuses, error messages,
                        and completion percentages.
        """
        logs = self.retrieve_logs()

        current_status["ID"] = (
            current_status["ID"]
            .str.replace("NP_", "NP")
            .str.replace("NM_", "NM")
            .str.split("_")
            .str[0]
        )

        summary = (
            current_status.groupby(["index", "ID"])[
                ["num_mutations_faa", "num_mutations_scored"]
            ]
            .sum()
            .reset_index()
        )

        summary["index"] = summary["index"].apply(
            lambda x: str(x) if "job_" not in x else x.split("_")[-1]
        )

        summary["percent"] = round(
            (summary["num_mutations_scored"] / summary["num_mutations_faa"]) * 100, 2
        )
        summary["remaining_mutations"] = (
            summary["num_mutations_faa"] - summary["num_mutations_scored"]
        )

        summary = summary.merge(logs, on="index", how="left")

        return summary

    def mark_as_completed(self, job):
        """
        Marks a MutPred2 job as completed by creating a 'completed.txt' file in the job directory.

        Args:
            job (str): The job identifier.

        Returns:
            bool: True if the job was successfully marked as completed, otherwise False.
        """
        if self.get_job_dir() != "":
            with open(
                f"{self.get_job_dir()}/{job}/completed.txt", "w", encoding="utf-8"
            ) as complete:
                complete.write("completed")
                return True

        return False

    def create_scripts(self, job_status, error):
        """
        Generates new job scripts to rerun failed MutPred2 jobs based on error type.

        Args:
            job_status (pd.DataFrame): DataFrame containing job status information.
            error (str): The type of error encountered (e.g., "Memory", "Time").

        Returns:
            None
        """
        script = self.get_script()
        script_dir = "/".join(script.split("/")[:-1])

        if error == "Memory":
            output_script = f"{script_dir}/mutpred2_high_mem_remaining.lsf"
            with open(output_script, "w", encoding="utf-8") as out:
                with open(script, "r", encoding="utf-8") as s:
                    for line in s:
                        if "#BSUB -R rusage[mem=" in line:
                            mem_line = "#BSUB -R rusage[mem=6000]\n"
                            out.write(mem_line)

                        elif "#BSUB -J" in line:
                            job_line = f'{line.split("[")[0]}{self.unfinished_jobs(job_status)}\n'
                            out.write(job_line)

                        else:
                            out.write(line)

        elif error == "Time":
            output_script = f"{script_dir}/mutpred2_long_time_remaining.lsf"
            with open(output_script, "w", encoding="utf-8") as out:
                with open(script, "r", encoding="utf-8") as s:
                    for line in s:
                        if "#BSUB -q" in line:
                            out.write("#BSUB -q long\n")

                        elif "#BSUB -W" in line:
                            out.write("#BSUB -W 249:59\n")

                        elif "#BSUB -J" in line:
                            job_line = f'{line.split("[")[0]}{self.unfinished_jobs(job_status)}\n'
                            out.write(job_line)

                        else:
                            out.write(line)

    def unfinished_jobs(self, jobs):
        """
        Identifies and categorizes unfinished MutPred2 jobs.

        Args:
            jobs (pd.DataFrame): DataFrame containing job status information.

        Returns:
            list: A list of job indices that are incomplete and need to be rerun.
        """
        leftover_jobs = (
            jobs["index"]
            .drop_duplicates()
            .astype(str)
            .str.split("_")
            .str[-1]
            .astype(int)
        )

        job_arrays = lsf.build_job_array(leftover_jobs)

        return job_arrays

    def show_job_summary(self, job):
        """
        Displays a summary of a MutPred2 job's progress using a visual progress bar.

        Args:
            job (pd.Series): A Pandas Series containing job status details, including
                            completion percentage and error messages.

        Returns:
            None
        """
        # print (job)
        status_bar_size = 25
        percent_complete = int((job["percent"] / 100.0) * status_bar_size)
        percent_incomplete = status_bar_size - percent_complete

        status_bar = "#" * percent_complete + " " * percent_incomplete
        # print (job)

        job_status = f"""Job {job['index']} [{status_bar}] {int(job['percent'])}%\t{job['Error']}"""
        print(job_status)

    def mutpred_summary(self):
        """
        Generates and displays a summary of MutPred2 job statuses, including
        completed, partially completed, and unscored jobs.

        The summary includes:
        - Total scored and remaining mutations.
        - Breakdown of job completion statuses.
        - Error categorization and troubleshooting suggestions.

        If a summary file is available, it will be loaded; otherwise, a new summary
        will be generated and saved.

        Returns:
            None
        """
        if self.get_summary():
            summary = pd.read_csv(
                self.get_summary(),
                sep="\t",
                dtype={
                    "num_mutations_faa": int,
                    "num_mutations_scored": int,
                    "percent": float,
                },
            )
        else:
            # print ("start")
            cur_mutpred_status = self.mutpred_status()
            # print (status)

            # print ("before summary")
            summary = self.mutpred_logs(cur_mutpred_status).sort_values("index")
            # print (summary)

            # print ("after summary")
            summary["hasError"].fillna(False, inplace=True)
            summary.loc[summary["percent"] == 100, "Error"] = ""
            summary["Error"].fillna("Unknown", inplace=True)

            summary = summary.sort_values("index")

            cur_date = datetime.today().strftime("%m-%d-%Y_%H:%M:%S")
            summary.to_csv(f"mutpred2.summary.{cur_date}.csv", sep="\t", index=False)

        # print (summary)

        gene_summary = (
            summary.groupby("ID")[
                ["num_mutations_scored", "remaining_mutations", "num_mutations_faa"]
            ]
            .sum()
            .reset_index()
        )
        gene_summary["percent"] = round(
            (gene_summary["num_mutations_scored"] / gene_summary["num_mutations_faa"])
            * 100,
            2,
        )
        # print (gene_summary)

        job_summary = (
            summary.groupby("index")[
                ["num_mutations_scored", "remaining_mutations", "num_mutations_faa"]
            ]
            .sum()
            .reset_index()
        )
        job_summary["percent"] = round(
            (job_summary["num_mutations_scored"] / job_summary["num_mutations_faa"])
            * 100,
            2,
        )
        # print (job_summary)
        print(
            f"""
    ===== Mutation Summary =====
    Scored Mutations:    {sum(summary['num_mutations_scored'])}
    Remaining Mutations: {sum(summary['remaining_mutations'])}

    ===== Job Summary =====
    Completed Jobs:           {len(job_summary[job_summary['percent'] == 100])}
    Partially Completed Jobs: \
        {len(job_summary[(job_summary['percent'] < 100) & job_summary['percent'] > 0])}
    Non-starter Jobs:         {len(job_summary[job_summary['percent'] == 0])}

    ===== Protein Summary =====
    Fully Scored Genes:     {len(gene_summary[gene_summary['percent'] == 100])}
    Partially Scored Genes: \
        {len(gene_summary[(gene_summary['percent'] < 100) & gene_summary['percent'] > 0])}
    Unscored Genes:         {len(gene_summary[gene_summary['percent'] == 0])}
        """
        )

        summary.loc[summary["percent"] == 100, "Error"] = ""
        summary["Error"].fillna("Unknown", inplace=True)

        summary["Type"] = summary["Error"].map(
            {
                "TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.": "Memory",
                "MATLAB is exiting because of fatal error: Memory Error": "Memory",
                "TERM_RUNLIMIT: job killed after reaching LSF run time limit.": "Time",
                "Unknown": "Unknown",
            }
        )

        error_types = (
            summary[summary["percent"] == 0]
            .groupby("Type")["index"]
            .count()
            .reset_index()
        )

        if self.__show_incomplete:
            print(
                self.unfinished_jobs(
                    summary[(summary["percent"] > 0) & (summary["percent"] < 100)]
                )
            )
        else:

            # print (summary[summary["percent"] < 100])

            for error_type in error_types["Type"]:
                print("Error:", error_type)
                print(
                    self.unfinished_jobs(
                        summary[
                            (summary["percent"] < 100) & (summary["Type"] == error_type)
                        ]
                    )
                )
                print(" ")

        if self.write_script_status():
            for error in error_types["Type"]:
                self.create_scripts(
                    summary[(summary["percent"] < 100) & (summary["Type"] == error)],
                    error,
                )
        else:
            pass


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Check the status of a currently running or a previously ran MutPred2 job."
    )

    parser.add_argument(
        "--job_dir", type=str, nargs="?", help="Path to the job directory for MutPred2"
    )
    parser.add_argument(
        "--show_incomplete",
        action="store_true",
        help="When listing the remaining jobs, show all incomplete jobs and \
            not just jobs with zero output.",
    )

    args = parser.parse_args()

    project = args.project
