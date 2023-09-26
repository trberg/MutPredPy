import pandas as pd
import argparse
import re
import os

from . import fasta
from . import fasta
from . import prep
from . import lsf

#import fasta
#import prep
#import lsf


class Remaining:
    def __init__(self, input, project, exclude, time, dry_run):
        
        self.__intermediate_dir = "intermediates"
        self.__project = project
        self.__input = self.input_path(input)
        self.__base = input.split("/")[-1].split(".")[0]
        self.__exclude_indices = self.generate_exclude_indices(exclude)

        self.time = time
        self.dry_run = dry_run

        self.__fasta_location = "resources/Homo_sapiens.GRCh38.combined.pep.all.fa"

    
    def get_intermediate_dir(self):
        return self.__intermediate_dir
    
    def get_project(self):
        return self.__project
    
    def get_base(self):
        return self.__base
    
    def get_exclude_indices(self):
        return self.__exclude_indices
    
    def generate_exclude_indices(self, exclude):
        
        exclude = re.sub(r'\[|\]', '', exclude)
        exclude = exclude.split(",")

        indices = []
        for e in exclude:
            if '-' in e:
                indices += list(range(int(e.split('-')[0]),int(e.split('-')[1])+1))
            elif e == '':
                pass
            else:
                indices.append(int(e))
        indices = [str(i) for i in indices]
        #print (indices)
        return indices


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
        f_exc = re.compile(f"{self.get_base()}.missense_({'|'.join(self.get_exclude_indices())}).faa")
        

        faa_files = [f for f in os.listdir(faa_dir) if f_reg.match(f) and not f_exc.match(f)]
        #faa_files_exclude = [f for f in os.listdir(faa_dir) if f_exc.match(f)]

        try:
            return pd.concat([fasta.read_mutpred_input_fasta(f"{faa_dir}/{file}") for file in faa_files])
        except ValueError:
            return pd.DataFrame()



    def read_mutpred_output(self, file):

        index = file.split("/")[-1].split(".")[-2].split("_")[-1]

        if os.path.isfile(file) and os.path.getsize(file) > 0:
            scores = pd.read_csv(file)

            scores = pd.DataFrame(scores.groupby("ID")["Substitution"].agg(set)).reset_index()
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
        #print (inputs)
        
        gene_of_interest = "|TTN"

        inputs = inputs.groupby("ID").agg({'mutations':lambda x: set.union(*x)}).reset_index()
        inputs = inputs[inputs["ID"].str.contains(gene_of_interest)]
        #print (inputs)
        
        
        outputs = self.retrieve_outputs()
        
        outputs = outputs.groupby("ID").agg({'Substitution':lambda x: set.union(*x)}).reset_index()
        outputs = outputs[outputs["ID"].str.contains(gene_of_interest)]
        #print (outputs)

        both = inputs.merge(outputs, on="ID", suffixes=["_faas", "_scored"], how="left")
        both["Substitution"] = both["Substitution"].fillna("").apply(set)
        
        both["pre_filter_count"] = both["Substitution"].apply(len)
        print ("Pre Filter:",sum(both["pre_filter_count"]))
        

        both["mutations"] = both["mutations"] - both["Substitution"]
        both = both.groupby("ID").agg({'mutations':lambda x: set.union(*x)}).reset_index()

        both["num_mutations_faas"] = both["mutations"].apply(len)
        both = both[both["num_mutations_faas"]>0]
        
        both[["Ensembl_proteinid","gene_symbol"]] = both["ID"].str.split("|", expand=True)
        both = both[["ID","Ensembl_proteinid","gene_symbol","mutations","num_mutations_faas"]]
        both = both.rename(columns={"num_mutations_faas":"num_mutations"})


        fasta_file = fasta.collect_fasta(self.__fasta_location)
        fasta_file = fasta_file[["Ensembl_proteinid_v","sequence","Memory Estimate (MB)","Time per Mutation (hrs)"]]
        fasta_file = fasta_file.rename(columns={"Ensembl_proteinid_v":"Ensembl_proteinid"})

        both = both.merge(fasta_file, on="Ensembl_proteinid", how="left")
        both["mutations"] = both["mutations"].apply(lambda x: " ".join(sorted(x)))
        both = both.rename(columns={"mutations":"mutation"})
        
        both = both.drop_duplicates()
        
        both["Time Estimate (hrs)"] = both.apply(lambda row: row["num_mutations"]*row["Time per Mutation (hrs)"], axis=1)

        #print (set(both["Ensembl_proteinid"]))
        TTN = both[both["gene_symbol"]=="TTN"]#[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]]
        print (TTN)
        print (sum(TTN["num_mutations"]))
        #exit()
        non_TTN = both[both["gene_symbol"]!="TTN"]#[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]]
        print (non_TTN[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]].sort_values("Time Estimate (hrs)"))
        both = TTN

        print (both[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]].sort_values("Time Estimate (hrs)"))
        #exit()
        print (f"Memory: {max(both['Memory Estimate (MB)'])}")
        print (f"Time: {sum(both['Time Estimate (hrs)'])}")
        print (f"Mutations: {sum(both['num_mutations'])}")
        
        return both


    def split_and_build_lsf(self):
        
        remaining = self.remainder()
        print (remaining[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]].sort_values("Time Estimate (hrs)"))
        #exit()

        mut = prep.MutPredpy(
            input=self.__input,
            project=self.__project,
            time=self.time,
            dry_run=self.dry_run,
            canonical=False
        )
        file_number = max([int(f.split(".")[-2].split("_")[-1]) for f in os.listdir(mut.set_faa_output())]) + 1
        
        tech_requirements = mut.split_data(remaining, file_number=file_number)

        lsf.usage_report(tech_requirements)
        lsf.build_lsf_config_file(tech_requirements, mut.get_intermediate_dir(), mut.get_project(), mut.get_base(), user=4, dry_run=mut.dry_run)



if __name__ == "__main__":

    def time_minimum(x):
        x = int(x)

        if x < 5:
            raise argparse.ArgumentTypeError("Minimum time is 5 hours")
        return x

    parser = argparse.ArgumentParser(description='Check the status of a currently running or a previously ran MutPred2 job.')
    

    parser.add_argument('--input', type=str, nargs='?',
                        help='The name of the input filename that is located in the data folder.')
    parser.add_argument('--project', type=str, nargs='?',
                        help='The name of the project for organization purposes')
    parser.add_argument('--exclude', type=str, nargs='?', required=False, default='[]',
                        help='LSF sequence of jobs to exlude from the remainder function')
    parser.add_argument('--time', type=time_minimum, nargs="?", default=24,
                        help="Target time in hours to run the jobs")
    parser.add_argument("--dry_run", action="store_true")

    
    args = parser.parse_args()
    
    R = Remaining(
        input=args.input,
        project=args.project,
        exclude=args.exclude,
        time=args.time,
        dry_run=args.dry_run
    )
    
    remaining_data = R.split_and_build_lsf()


    