import pandas as pd
import argparse
import re
import os

#from . import fasta
import fasta


class Remaining:
    def __init__(self, input, project):
        
        self.__intermediate_dir = "intermediates"
        self.__project = project
        self.__input = self.input_path(input)
        self.__base = input.split("/")[-1].split(".")[0]

        self.__fasta_location = "resources/Homo_sapiens.GRCh38.combined.pep.all.fa"

    
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
    


    def retrieve_faas(self):
        faa_dir = f"{self.get_intermediate_dir()}/faa/{self.get_project()}"
        f_reg = re.compile(f"{self.get_base()}.missense_\d+.faa")

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f)]

        return pd.concat([fasta.read_mutpred_input_fasta(f"{faa_dir}/{file}") for file in faa_files])



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
    


    def retrieve_outputs(self):
        o_reg = re.compile(f"{self.get_base()}.missense_output_\d+.txt$")
        out_dir = f"{self.get_intermediate_dir()}/scores"
        
        out_files = [o for o in os.listdir(out_dir) if o_reg.match(o)]

        return pd.concat([self.read_mutpred_output(f"{out_dir}/{s}") for s in out_files])
    


    def remainder(self):

        inputs  = self.retrieve_faas()
        inputs["mutations"] = inputs["mutations"].str.split(",").apply(set)

        outputs = self.retrieve_outputs()

        both = inputs.merge(outputs, on="ID", suffixes=["_faas", "_scored"], how="left")
        both["Substitution"] = both["Substitution"].fillna("").apply(set)
        #print (both)

        both["mutations"] = both["mutations"] - both["Substitution"]
        both["num_mutations_faas"] = both["mutations"].apply(len)
        both = both[both["num_mutations_faas"]>0]
        both[["Ensembl_proteinid_v","gene_symbol"]] = both["ID"].str.split("|", expand=True)
        both = both[["ID","Ensembl_proteinid_v","gene_symbol","mutations","num_mutations_faas","index_faas","Substitution","num_mutations_scored"]]

        #print (both)

        fasta_file = fasta.collect_fasta(self.__fasta_location)
        fasta_file = fasta_file[["Ensembl_proteinid_v","sequence","Memory Estimate (MB)","Time per Mutation (hrs)"]]

        both = both.merge(fasta_file, on="Ensembl_proteinid_v", how="left")
        print (both)

        print (f"Memory: {max(both['Memory Estimate (MB)'])}")
        print (f"Time: {sum(both['Time per Mutation (hrs)'])}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Check the status of a currently running or a previously ran MutPred2 job.')
    

    parser.add_argument('--input', type=str, nargs='?',
                        help='The name of the input filename that is located in the data folder.')
    parser.add_argument('--project', type=str, nargs='?',
                        help='The name of the project for organization purposes')

    
    args = parser.parse_args()

    input   = args.input
    project = args.project
    
    R = Remaining(
        input=input,
        project=project
    )
    
    R.remainder()

    