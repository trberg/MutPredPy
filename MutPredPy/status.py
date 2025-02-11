import pandas as pd
import argparse
from datetime import datetime
from pkg_resources import Requirement, resource_filename

import re
import os

from . import fasta
from . import lsf


class Status:
    def __init__(self, job_dir, all, show_incomplete, summary="", logs="", script=False):
        
        self.__job_dir = job_dir
        self.__show_all = all
        self.__show_incomplete = show_incomplete
        self.__log_dir = logs
        self.__write_new_scripts = bool(script)
        if self.__write_new_scripts:
            self.__template_script = script


        ## Keyword   Term_reason  Integer_logged
        #self.__log_messages = pd.read_csv(os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/LSF_errors.csv")))
        
        if summary == "":
            self.__summary = False
        else:
            self.__summary = summary
    

    def write_script_status(self):
        return self.__write_new_scripts


    def get_script(self):
        return os.path.abspath(self.__template_script)
        

    def get_job_dir(self):
        return os.path.abspath(self.__job_dir.rstrip("/"))
     
    
    def get_log_dir(self):
        if self.__log_dir == "":
            return self.__log_dir
        else:
            return os.path.abspath(self.__log_dir.rstrip("/"))
    

    def get_summary(self):
        if not self.__summary:
            return False
        else:
            return os.path.abspath(self.__summary)
    

    def get_mutpred_output_file(self, jobid):

        job_dir = f"{self.get_job_dir()}/{jobid}"
        
        #o_reg = re.compile(f".*.missense_output_{jobid}.txt")
        o_reg = re.compile(f"output.txt")

        out_files = [o for o in os.listdir(job_dir) if o_reg.match(o)]

        if len(out_files) == 1:
            out_file = out_files[0]
        elif len(out_files) > 1:
            raise Exception(f"More than 1 output file for job {jobid}")
        elif len(out_files) == 0:
            raise Exception(f"No output file for job {jobid}")
        else:
            raise Exception(f"Error with reading output file for job {jobid}")

        return out_file
    

    def get_mutpred_output_file_path(self, index):

        directory = self.get_job_dir()
        file = self.get_mutpred_output_file(index)
        return f"{directory}/{file}"


    def get_mutpred_input_file(self, index):

        file = self.retrieve_faas(jobid=index)
        
        return file
    

    def get_mutpred_input_file_path(self, project, index):

        directory = self.get_job_dir()
        file = self.get_mutpred_input_file(index)

        return f"{directory}/{file}"
    

    def read_mutpred_output(self, file, jobid):

        index = file.split("/")[-1].split(".")[-2].split("_")[-1]
        
        if os.path.isfile(file) and os.path.getsize(file) > 0:
            
            scores = pd.read_csv(file)

            scores = scores.drop_duplicates()
            
            scores = pd.DataFrame(scores.groupby("ID")["Substitution"].agg(list)).reset_index()
            
            scores["num_mutations"] = scores["Substitution"].apply(len).astype(int)
            
            scores["index"] = jobid

            scores = scores[["ID","num_mutations","index"]]
            
        else:
            scores = pd.DataFrame(columns=["ID","num_mutations","index"])

        
        return scores
    
    
    def read_err_log(self, log):
        
        error = ""
        
        with open(f"{self.get_log_dir()}/{log}") as l:
            error_code = 0
            for line in l:
                if len(line) > 1 and ("TERM" in line):
                    return line.strip("\n")
                elif len(line) > 1 and ("MATLAB is exiting because of fatal error" in line):
                    next_line = next(l)
                    if "Killed" in next_line or "Bus" in next_line:
                        return line.strip("\n") + ": " + "Memory Error"
                    elif "Bus" in next_line:
                        return line.strip("\n") + ": " + "Memory Error"
                    else:
                        return line.strip("\n") + ": " + next_line.strip("\n") 
        return error

    
    
    def has_err_log(self, log):

        if self.get_log_dir() == "":
            error = ""
            return pd.Series([False, error])
        
        cur_log = f"{self.get_log_dir()}/{log}"
        
        if os.path.exists(cur_log) and os.path.getsize(cur_log) > 0:
            error = self.read_err_log(log)
            return pd.Series([True, error])
        
        else:
            error = ""
            return pd.Series([False, error])


    def get_job_index(self, log):

        with open(log) as L:
            for line in L:
                
                if line.startswith("Subject"):
                    job,index = line.split(":")[1].strip(" Job ").strip("]").split("[")
                    return pd.Series([job, index])
                    #Subject: Job 137347991[10]: <phase3[1-4000]> in cluster <chimera> Exited
            else:
                return ["",""]


    def retrieve_logs(self):
        
        print (self.get_log_dir())
        if not os.path.isdir(self.get_log_dir()):
            
            log_files = pd.DataFrame({"logs":[]})
            log_files["job"] = ""
            log_files["index"] = ""
            log_files["hasError"] = ""
            log_files["Error"] = ""
            log_files = log_files[["index", "hasError", "Error"]]

        else:
            
            log_files = pd.DataFrame({"logs":[l for l in os.listdir(f"{self.get_log_dir()}/") if "out" in l]})

            log_files[["job","index"]] = log_files["logs"].apply(lambda x: self.get_job_index(f"{self.get_log_dir()}/{x}"))
        
            latest_err_logs = pd.DataFrame(log_files.groupby(["index"])["job"].max()).reset_index()
            
            log_files = log_files.merge(latest_err_logs, on=["index","job"], how="inner")
            #print (log_files)
            
            log_files[["hasError","Error"]] = log_files.apply(lambda row: self.has_err_log(row["logs"]), axis=1)
            log_files = log_files[["index", "hasError", "Error"]]
            
            log_files["Error"].fillna("No Errors", inplace=True)
        
        return log_files
    

    def retrieve_faas(self, jobid):

        faa_dir = f"{self.get_job_dir()}/{jobid}"
        
        f_reg = re.compile(f".*.faa")

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f)]

        if len(faa_files) == 1:
            faa_file = faa_files[0]
        elif len(faa_files) > 1:
            raise Exception(f"More than 1 input fasta file in job {jobid}")
        elif len(faa_files) == 0:
            raise Exception(f"No input fasta file found for job {jobid}")
        else:
            raise Exception(f"Error finding input fasta file for job {jobid}")

        job_faa = fasta.read_mutpred_input_fasta(f"{faa_dir}/{faa_file}")

        job_faa["index"] = jobid

        return job_faa
    

    def retrieve_outputs(self, jobid):

        job_dir = f"{self.get_job_dir()}/{jobid}"
        
        #o_reg = re.compile(f".*.missense_output_{jobid}.txt$")
        o_reg = re.compile(f"output.txt$")

        out_files = [o for o in os.listdir(job_dir) if o_reg.match(o)]
        #if jobid == 71

        if len(out_files) == 1:
            out_file = out_files[0]
        elif len(out_files) == 0:
            out_file = "no_out_file.txt"
        elif len(out_files) > 1:
            raise Exception(f"More than 1 MutPred2 output file for job {jobid}")
        else:
            raise Exception(f"Error in reading output file for job {jobid}")

        return self.read_mutpred_output(f"{job_dir}/{out_file}", jobid)

    
    def mutpred_status(self):

        all_statuses = []
        
        for job in os.listdir(f"{self.get_job_dir()}"):

            if os.path.isdir(f"{self.get_job_dir()}/{job}"):

                if not os.path.exists(f"{self.get_job_dir()}/{job}/completed.txt"):
                    
                    faa = self.retrieve_faas(job)
                    #print (faa)

                    scores = self.retrieve_outputs(job)
                    #print (scores)
                    
                    #exit()

                    status = faa.merge(scores, on=["ID", "index"], how="outer", suffixes=("_faa","_scored")).fillna(0)
                    status["complete"] = status.apply(lambda x: x["num_mutations_faa"]==x["num_mutations_scored"], axis=1)

                    if sum(status["num_mutations_faa"]) == sum(status["num_mutations_scored"]):
                        self.mark_as_completed(job)
                        
                    all_statuses.append(status)
                    
                else:

                    faa = self.retrieve_faas(job)
                    
                    faa.rename(columns={"num_mutations":"num_mutations_faa"}, inplace=True)
                    faa["num_mutations_scored"] = faa["num_mutations_faa"]
                    #faa["percent"] = 100.0
                    faa["complete"] = True

                    all_statuses.append(faa)


        all_statuses = pd.concat(all_statuses)

        return all_statuses

    
    def mutpred_logs(self, status):

        logs = self.retrieve_logs()
    
        status["ID"] = status["ID"].str.replace("NP_","NP").str.replace("NM_","NM").str.split("_").str[0]

        summary = status.groupby(["index","ID"])[["num_mutations_faa","num_mutations_scored"]].sum().reset_index()

        
        summary["index"] = summary["index"].apply(lambda x: str(x) if "job_" not in x else x.split("_")[-1])

        summary["percent"] = round((summary["num_mutations_scored"]/summary["num_mutations_faa"])*100, 2)
        summary["remaining_mutations"] = summary["num_mutations_faa"] - summary["num_mutations_scored"]

        summary = summary.merge(logs, on="index", how="left")
        
        return summary


    def mark_as_completed(self, job):

        if self.get_job_dir() != "":
            with open(f"{self.get_job_dir()}/{job}/completed.txt", "w") as complete:
                complete.write("completed")
                return True    
            
        return False
            

    def create_scripts(self, job_status, error):

        script = self.get_script()
        script_dir = "/".join(script.split("/")[:-1])

        if error == "Memory":
            output_script = f"{script_dir}/mutpred2_high_mem_remaining.lsf"
            with open(output_script, "w") as out:
                with open(script) as s:
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
            with open(output_script, "w") as out:
                with open(script) as s:
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
        
        leftover_jobs = jobs["index"].drop_duplicates().astype(str).str.split("_").str[-1].astype(int)
        
        job_arrays = lsf.build_job_array(leftover_jobs)
        
        return job_arrays
    

    def show_job_summary(self, job):
        #print (job)
        status_bar_size = 25
        percent_complete = int((job['percent']/100.0)*status_bar_size)
        percent_incomplete = status_bar_size - percent_complete
        
        status_bar = "#"*percent_complete + " "*percent_incomplete
        #print (job)
        
        job_status = f"""Job {job['index']} [{status_bar}] {int(job['percent'])}%\t{job['Error']}"""
        print (job_status)


    def mutpred_summary(self):
        
        if self.get_summary():
            summary = pd.read_csv(self.get_summary(), sep="\t", dtype={
                "num_mutations_faa":int,
                "num_mutations_scored":int,
                "percent":float
            })
        else:
            #print ("start")
            status = self.mutpred_status()
            #print (status)
            
            #print ("before summary")
            summary = self.mutpred_logs(status).sort_values("index")
            #print (summary)
            
            #print ("after summary")
            summary["hasError"].fillna(False, inplace=True)
            summary.loc[summary["percent"] == 100, "Error"] = ""
            summary["Error"].fillna("Unknown", inplace=True)

            summary = summary.sort_values("index")

            cur_date = datetime.today().strftime('%m-%d-%Y_%H:%M:%S')
            summary.to_csv(f"mutpred2.summary.{cur_date}.csv", sep="\t", index=False)

        #print (summary)

        gene_summary = summary.groupby("ID")[["num_mutations_scored","remaining_mutations","num_mutations_faa"]].sum().reset_index()
        gene_summary["percent"] = round((gene_summary["num_mutations_scored"]/gene_summary["num_mutations_faa"])*100, 2)
        #print (gene_summary)

        job_summary = summary.groupby("index")[["num_mutations_scored","remaining_mutations","num_mutations_faa"]].sum().reset_index()
        job_summary["percent"] = round((job_summary["num_mutations_scored"]/job_summary["num_mutations_faa"])*100, 2)
        #print (job_summary)
        print(f"""
    ===== Mutation Summary =====
    Scored Mutations:    {sum(summary['num_mutations_scored'])}
    Remaining Mutations: {sum(summary['remaining_mutations'])}

    ===== Job Summary =====
    Completed Jobs:           {len(job_summary[job_summary['percent'] == 100])}
    Partially Completed Jobs: {len(job_summary[(job_summary['percent'] < 100) & job_summary['percent'] > 0])}
    Non-starter Jobs:         {len(job_summary[job_summary['percent'] == 0])}

    ===== Protein Summary =====
    Fully Scored Genes:     {len(gene_summary[gene_summary['percent'] == 100])}
    Partially Scored Genes: {len(gene_summary[(gene_summary['percent'] < 100) & gene_summary['percent'] > 0])}
    Unscored Genes:         {len(gene_summary[gene_summary['percent'] == 0])}
        """)
        
        summary.loc[summary["percent"] == 100, "Error"] = ""
        summary["Error"].fillna("Unknown", inplace=True)

        summary["Type"] = summary["Error"].map({
            "TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.": "Memory",
            "MATLAB is exiting because of fatal error: Memory Error":"Memory",
            "TERM_RUNLIMIT: job killed after reaching LSF run time limit.": "Time",
            "Unknown":"Unknown"
        })
        

        error_types = summary[summary["percent"] == 0].groupby("Type")["index"].count().reset_index()

        if self.__show_incomplete:
            print (self.unfinished_jobs(summary[(summary["percent"] > 0) & (summary["percent"] < 100)]))
        else:
            
            #print (summary[summary["percent"] < 100])

            for error_type in error_types["Type"]:
                print ("Error:", error_type)
                print (self.unfinished_jobs(summary[(summary["percent"] < 100) & (summary["Type"] == error_type)]))
                print (" ")
        
        
        if self.write_script_status():
            for error in error_types["Type"]:
                self.create_scripts(summary[(summary["percent"] < 100) & (summary["Type"] == error)], error)
        else:
            pass



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Check the status of a currently running or a previously ran MutPred2 job.')
    
    parser.add_argument('--job_dir', type=str, nargs='?', 
                        help='Path to the job directory for MutPred2')
    parser.add_argument('--show_incomplete', action="store_true", 
                        help='When listing the remaining jobs, show all incomplete jobs and not just jobs with zero output.')
    
    args = parser.parse_args()

    project = args.project
    
    status = Status(project=project)
    
    status.mutpred_status()