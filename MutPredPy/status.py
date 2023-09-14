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

            scores = pd.DataFrame(scores.groupby("ID")["Substitution"].agg(list)).reset_index()
            scores["num_mutations"] = scores["Substitution"].apply(len).astype(int)
            scores["index"] = index
        else:
            scores = pd.DataFrame(columns=["ID","Substitution","num_mutations","index"])

        return scores


    
    def mutpred_status(self):

        base = self.get_base()
        project = self.get_project()
        intermediate = self.get_intermediate_dir()

        
        faa_dir = f"{intermediate}/faa/{project}"
        f_reg = re.compile(f"{base}.missense_\d+.faa")

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f)]

        faa = pd.concat([fasta.read_mutpred_input_fasta(f"{faa_dir}/{file}") for file in faa_files]) 
        


        o_reg = re.compile(f"{base}.missense_output_\d+.txt$")
        out_dir = f"{intermediate}/scores"
        
        out_files = [o for o in os.listdir(out_dir) if o_reg.match(o)]

        scores = pd.concat([self.read_mutpred_output(f"{out_dir}/{s}") for s in out_files])

        #print (scores)

        status = faa.merge(scores, on=["ID","index"], how="outer", suffixes=("_faa","_scored")).fillna(0)
        status["complete"] = status.apply(lambda x: x["num_mutations_faa"]==x["num_mutations_scored"], axis=1)

        summary = pd.DataFrame(status.groupby("index")["num_mutations_faa","num_mutations_scored"].agg(sum)).reset_index()
        summary["index"] = summary["index"].astype(int)

        summary["percent"] = round((summary["num_mutations_scored"]/summary["num_mutations_faa"])*100, 2)
        summary["remaining_mutations"] = summary["num_mutations_faa"] - summary["num_mutations_scored"]

        summary = summary.sort_values("index")

        for index, row in summary.iterrows():
            if row['percent'] < 100:
                print (f"{base}.missense_{int(row['index'])}\t{row['percent']}%")

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