import pandas as pd
import argparse

import re
import os

from . import fasta


class Status:
    def __init__(self, input, project):
        
        self.__intermediate_dir = "intermediates"
        self.__project = project
        self.__input = self.input_path(input)
        self.__base = input.split("/")[-1].split(".")[0]

    
    def get_intermediate_dir(self):
        return self.__intermediate_dir
    
    def get_project(self):
        return self.__project
    
    def get_base(self):
        return self.__base


    def input_path(self, input):
        
        if os.path.exists(input):
            path = input
        
        else:
            print (f"File input {input} not found...trying inputs/{self.__project}/{input}")
            path = f"inputs/{self.__project}/{input}"
        
        
        if os.path.exists(path):
            return path
        
        else:
            print (f"Input file {path} not found.")
            exit()



    def read_mutpred_output(self, file):

        index = file.split("/")[-1].split(".")[-2].split("_")[-1]

        if os.path.isfile(file) and os.path.getsize(file) > 0:
            scores = pd.read_csv(file)
            scores = scores.drop_duplicates()

            scores = pd.DataFrame(scores.groupby("ID")["Substitution"].agg(list)).reset_index()
            scores["num_mutations"] = scores["Substitution"].apply(len).astype(int)
            scores["index"] = index
        else:
            scores = pd.DataFrame(columns=["ID","Substitution","num_mutations","index"])

        return scores
    

    
    def read_err_log(self, log):
        
        error = ""
        
        with open(f"logs/{self.get_project()}/{log}") as l:
            for line in l:
                if len(line) > 1:
                    error += line
        
        return error

    
    
    
    def has_err_log(self, log, index, job):
        cur_log = f"logs/{self.get_project()}/{log}"
        #print (cur_log, os.path.exists(cur_log), os.path.getsize(cur_log))
        
        if os.path.exists(cur_log) and os.path.getsize(cur_log) > 0:
            error = self.read_err_log(log).replace("\n",", ")
            return pd.Series([True, error])
        
        else:
            error = ""
            return pd.Series([False, error])



    def retrieve_logs(self):
        l_reg = re.compile(f"err_{self.get_base()}.\d+.faa_file_\d+$")

        log_files = pd.DataFrame({"logs":[l for l in os.listdir(f"logs/{self.get_project()}/") if l_reg.match(l)]})
        log_files["type"] = log_files["logs"].str.split("_").str[0]
        log_files["job"] = log_files["logs"].str.split(".").str[1]
        log_files["index"] = log_files["logs"].str.split("_").str[-1].astype(int)
        
        latest_err_logs = pd.DataFrame(log_files.groupby(["index","type"])["job"].max()).reset_index()
        
        log_files = log_files.merge(latest_err_logs, on=["index","type","job"], how="inner")

        print (log_files)
        print (latest_err_logs)
        exit()

        latest_job = max(log_files["job"])

        log_files = log_files[log_files["job"]==latest_job]

        log_files[["hasError","Error"]] = log_files.apply(lambda row: self.has_err_log(row["logs"], row["index"], latest_job), axis=1)
        log_files = log_files[["index", "hasError", "Error"]]
        return log_files
    

    def retrieve_faas(self):
        faa_dir = f"{self.get_intermediate_dir()}/faa/{self.get_project()}"
        f_reg = re.compile(f"{self.get_base()}.missense_\d+.faa")

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f)]

        return pd.concat([fasta.read_mutpred_input_fasta(f"{faa_dir}/{file}") for file in faa_files])
    

    def retrieve_outputs(self):
        o_reg = re.compile(f"{self.get_base()}.missense_output_\d+.txt$")
        out_dir = f"{self.get_intermediate_dir()}/scores"
        
        out_files = [o for o in os.listdir(out_dir) if o_reg.match(o)]

        return pd.concat([self.read_mutpred_output(f"{out_dir}/{s}") for s in out_files])

    
    def mutpred_status(self):

        base = self.get_base()

        faa = self.retrieve_faas()
        
        scores = self.retrieve_outputs()

        logs = self.retrieve_logs()

        status = faa.merge(scores, on=["ID","index"], how="outer", suffixes=("_faa","_scored")).fillna(0)
        status["complete"] = status.apply(lambda x: x["num_mutations_faa"]==x["num_mutations_scored"], axis=1)

        summary = pd.DataFrame(status.groupby("index")["num_mutations_faa","num_mutations_scored"].agg(sum)).reset_index()
        summary["index"] = summary["index"].astype(int)

        summary["percent"] = round((summary["num_mutations_scored"]/summary["num_mutations_faa"])*100, 2)
        summary["remaining_mutations"] = summary["num_mutations_faa"] - summary["num_mutations_scored"]

        summary = summary.merge(logs, on="index", how="left")
        summary["hasError"].fillna(False, inplace=True)
        summary["Error"].fillna("", inplace=True)

        summary = summary.sort_values("index")

        for index, row in summary.iterrows():
            if row['percent'] < 100 and not row['hasError']:
                print (f"{base}.missense_{int(row['index'])}\t{row['percent']}%")
            elif row['percent'] < 100 and row['hasError']:
                print (f"{base}.missense_{int(row['index'])}\t{row['percent']}%\t{row['Error']}")
            else:
                print (f"{base}.missense_{int(row['index'])}\t{row['percent']}%\tComplete!")

        print (f"> Remaining Mutations {sum(summary['remaining_mutations'])}")




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Check the status of a currently running or a previously ran MutPred2 job.')
    

    parser.add_argument('--input', type=str, nargs='?',
                        help='The name of the input filename that is located in the data folder.')
    parser.add_argument('--project', type=str, nargs='?',
                        help='The name of the project for organization purposes')

    
    args = parser.parse_args()

    input   = args.input
    project = args.project
    
    status = Status(
        input=input,
        project=project
    )
    
    status.mutpred_status()