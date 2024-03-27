from pkg_resources import Requirement, resource_filename

import pandas as pd
import argparse
import json
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
    def __init__(self, working_dir, time, dry_run):
        
        self.__working_dir = working_dir
        self.time = time
        self.dry_run = dry_run


        #self.__fasta_location = "resources/Homo_sapiens.GRCh38.combined.pep.all.fa"
        self.__fasta_location = os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/Homo_sapiens.GRCh38.combined.pep.all.fa"))

    
    def get_working_dir(self):
        return os.path.abspath(f"{self.__working_dir}")
    

    def get_job_dir(self):
        return os.path.abspath(f"{self.__working_dir}/jobs")
    
    
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


    def read_mutpred_output(self, file):

        cur_job = file.split("/")[-2].split("_")[-1]

        if os.path.isfile(file) and os.path.getsize(file) > 0:
            scores = pd.read_csv(file)

            scores = pd.DataFrame(scores.groupby("ID")["Substitution"].agg(set)).reset_index()
            scores["num_mutations"] = scores["Substitution"].apply(len).astype(int)
            scores["job"] = cur_job
        else:
            scores = pd.DataFrame(columns=["ID","Substitution","num_mutations","job"])

        return scores
    


    def retrieve_outputs(self):
        o_reg = re.compile(f"{self.get_base()}.missense_output_\d+.txt$")
        out_dir = f"{self.get_intermediate_dir()}/scores"
        
        out_files = [o for o in os.listdir(out_dir) if o_reg.match(o)]

        return pd.concat([self.read_mutpred_output(f"{out_dir}/{s}") for s in out_files])
    

    def collect_input_file_name(self, job):
        r = re.compile(r'.*.faa')
        input_files = list(filter(r.match,os.listdir(f"{self.get_job_dir()}/{job}")))
        if len(input_files) > 1:
            return None
        elif len(input_files) == 0:
            return None
        else:
            return input_files[0]


    def collect_output_file_name(self, job):
        r = re.compile(r'.*.txt$')
        output_files = list(filter(r.match,os.listdir(f"{self.get_job_dir()}/{job}")))
        if len(output_files) > 1:
            return None
        elif len(output_files) == 0:
            return None
        else:
            return output_files[0]


    def retrieve_jobs(self):

        inputs = []
        outputs = []

        for job in os.listdir(self.get_job_dir()):
            input_file = self.collect_input_file_name(job)
            if input_file != None:
                inputs.append(fasta.read_mutpred_input_fasta(f"{self.get_job_dir()}/{job}/{input_file}"))

            output_file = self.collect_output_file_name(job)
            if output_file != None:
                outputs.append(self.read_mutpred_output(f"{self.get_job_dir()}/{job}/{output_file}"))
        

        inputs = pd.concat(inputs)
        inputs["mutations"] = inputs["mutations"].str.split(",").apply(set)

        outputs = pd.concat(outputs)

        return inputs, outputs



    def remainder(self):

        inputs, outputs  = self.retrieve_jobs()
        
        #exit()

        inputs = inputs.groupby("ID").agg({'mutations':lambda x: set.union(*x)}).reset_index()
        outputs = outputs.groupby("ID").agg({'Substitution':lambda x: set.union(*x)}).reset_index()
        
        both = inputs.merge(outputs, on="ID", suffixes=["_faas", "_scored"], how="left")
        both["Substitution"] = both["Substitution"].fillna("").apply(set)
        
        both["pre_filter_count"] = both["Substitution"].apply(len)
        print ("Pre Filter:",sum(both["pre_filter_count"]))
        

        both["mutations"] = both["mutations"] - both["Substitution"]

        both = both.groupby("ID").agg({'mutations':lambda x: set.union(*x)}).reset_index()
        #print (both)
        both["num_mutations_faas"] = both["mutations"].apply(len)
        both = both[both["num_mutations_faas"] > 0]

        if len(both) > 0:
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

            
            print (f"Memory: {max(both['Memory Estimate (MB)'])}")
            print (f"Time: {sum(both['Time Estimate (hrs)'])}")
            print (f"Mutations: {sum(both['num_mutations'])}")
        else:
            print ("No remaining mutations to score")
            exit()
        
        return both


    def split_and_build_lsf(self):
        
        remaining = self.remainder()
        #print (remaining[["ID","Ensembl_proteinid","mutation","num_mutations","Time Estimate (hrs)"]].sort_values("Time Estimate (hrs)"))
        
        ##def __init__(self, input, working_dir, time, dry_run, canonical, database, fasta_location):
        mut = prep.MutPredpy(
            input="",
            working_dir=self.get_job_dir(),
            time=self.time,
            dry_run=self.dry_run,
            canonical=False,
            database="None"
        )
        file_number = max([int(f.split("_")[-1]) for f in os.listdir(self.get_job_dir())]) + 1
        
        tech_requirements = mut.split_data(remaining, file_number=file_number)

        lsf.usage_report(tech_requirements)
        ## build_lsf_config_file(tech_requirements, working_dir, base, user, dry_run)
        lsf.build_lsf_config_file(tech_requirements, self.get_job_dir(), "regeneron", user=1, dry_run=self.dry_run)



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


    