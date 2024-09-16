from . import utils as u
from . import fasta
import pandas as pd
import os


class Merge:
    def __init__(self, input, output, job_dir, database, assembly, dry_run):
        
        if "@" in database:
            self.config_file, self.config_name = database.split("@")
        else:
            self.config_file = database
            self.config_name = "Local"
        
        self.input = input
        self.output = output
        self.job_dir = job_dir
        self.dry_run = dry_run
        self.assembly = assembly
        self.scored = self.mutpred_scores()


    def write_output(self, data, output_file):
        data.to_csv(output_file, sep="\t", index=False)


    def mutpred_scores(self):


        scores = pd.concat([pd.read_csv(f"{m}/output.txt") for m in os.listdir(self.job_dir) if os.path.exists()])
        #engine = u.get_sql_engine(config_file=self.config_file, config_name=self.config_name)
        print (scores)
        exit()
        #scores = pd.read_sql("SELECT * FROM mutations", con=engine)
        #if self.assembly == "hg19":
        #    scores = pd.read_csv("../NYCKidSeqVUS/hg19.missense_output_2.txt", sep=",")
        #elif self.assembly == "hg38":
        #    scores = pd.read_csv("../NYCKidSeqVUS/hg38.missense_output_1.txt", sep=",")
        #scores["Ensembl_proteinid_v"] = scores["ID"].str.split("|").str[0]
        #scores["mutation"] = scores["Substitution"]

        return scores


    def merge(self, input, output):

        input_data = pd.read_csv(input, sep="\t", skiprows=[i for i,line in enumerate(open(input)) if line.startswith("##")])
        input_data[["protein_id","mutations"]] = input_data["Extra"].apply(lambda x: u.collect_value(x, "HGVSp")).str.split(":", expand=True)
        
        if self.assembly == "hg38":
            fasta_data = fasta.collect_fasta("resources/Homo_sapiens.GRCh38.combined.pep.all.fa")[["Ensembl_proteinid_v","sequence"]]
        elif self.assembly == "hg19":
            fasta_data = fasta.collect_fasta("resources/Homo_sapiens.GRCh37.combined.pep.all.fa")[["Ensembl_proteinid_v","sequence"]]
        else:
            print ("Assembly not recognized, defaulting to hg38.")
        
        
        fasta_data = fasta_data.drop_duplicates(subset="Ensembl_proteinid_v", keep="first")
        
        scores = self.scored
        
        input_data = input_data.merge(fasta_data, left_on="protein_id", right_on="Ensembl_proteinid_v", how="left")
        input_data["seq_hash"] = input_data["sequence"].apply(lambda x: u.get_seq_hash(x) if not pd.isna(x) else x)
        input_data["mutation"] = input_data["mutations"].apply(lambda x: fasta.mutation_mapping(x.split(".")[1]) if not pd.isna(x) else x)
        #input_data = input_data.merge(scores, on=["seq_hash","mutation"], how="left")
        input_data = input_data.merge(scores, on=["Ensembl_proteinid_v","mutation"], how="left")
        input_data = input_data[["#Uploaded_variation","Location","Allele","Gene","Feature","protein_id","mutations","MutPred2 score","Molecular mechanisms with Pr >= 0.01 and P < 1.00","Motif information","Remarks"]]

        print (input_data)

        self.write_output(input_data, output)
