import pandas as pd
import argparse
import numpy as np

import os
import re
import hashlib


from . import fasta
from . import lsf
from . import sql_connection


class MutPredpy:
    def __init__(self, input, project, time, dry_run, canonical, database, assembly):
        
        self.__intermediate_dir = "intermediates"
        self.__project = project
        self.__input = self.input_path(input)
        self.__base = input.split("/")[-1].split(".")[0]
        self.__time = time

        self.__fasta_file_number_start = self.set_fasta_file_number_start()

        self.dry_run   = dry_run
        self.canonical = canonical

        if "@" in database:
            self.database = database.split("@")[0]
            self.db_config_name = database.split("@")[1]
        else:
            self.database = database
            self.db_config_name = "Local"

        if assembly == "hg38":
            self.__fasta_location = "resources/Homo_sapiens.GRCh38.combined.pep.all.fa"
        elif assembly == "hg19":
            self.__fasta_location = "resources/Homo_sapiens.GRCh37.combined.pep.all.fa"
        else:
            print ("assembly not recognized, using hg38 assembly")
            self.__fasta_location = "resources/Homo_sapiens.GRCh38.combined.pep.all.fa"

        self.fasta = fasta.collect_fasta(self.__fasta_location)

        self.header_data, self.original_data = self.read_input_data()
        self.file_format = self.input_format()
        self.annotation = self.get_annotation()


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



    def get_input_data(self):
        return self.original_data
    
    def get_intermediate_dir(self):
        return self.__intermediate_dir
    
    def get_project(self):
        return self.__project
    
    def get_base(self):
        return self.__base
    
    def get_time(self):
        return self.__time
    
    def get_database_status(self):
        if self.database == "None":
            return False
        else:
            return True
        
    def get_database_connection(self):
        
        con = sql_connection.SQL_Connection(config_name=self.db_config_name, config_file=self.database)
        engine = con.get_engine()
        return engine
    
    
    def get_annotation(self):
        if self.file_format == "VCF-info":
            if "SnpEff" in self.header_data:
                self.__annotation = "SnpEff"
            elif "ENSEMBL VARIANT EFFECT PREDICTOR" in self.header_data:
                self.__annotation = "VEP"
            else:
                self.__annotation = False

        elif self.file_format == "dbNSFP":
            self.__annotation = "dbNSFP"

        else:
            self.__annotation = False
        
        return self.__annotation
    


    def input_format(self):
        cols = set(self.original_data.columns)
        if len(set(["#chr","pos(1-based)","ref","alt","aaref","aaalt","HGVSp_VEP","HGVSp_ANNOVAR","VEP_canonical"]).intersection(cols)) == 9:
            file_format = "dbNSFP"
        elif len(set(["Location","Allele","Consequence","Protein_position","Amino_acids","Extra"]).intersection(cols)) == 6:
            file_format = "VEP"
        elif len(set(["#CHROM","POS","REF","ALT","INFO"]).intersection(cols)) == 5:
            file_format = "VCF-info"
        elif len(set(["#CHROM","POS","REF","ALT"]).intersection(cols)) == 4:
            file_format = "VCF"
        else:
            file_format = "Custom"

        print (f"{file_format} file format")
        return file_format
    

    def set_fasta_file_number_start(self):

        project_output = f"{self.__intermediate_dir}/{self.__project}"
        
        if not os.path.exists(f"{project_output}"):
            return 1
        elif os.path.exists(f"{project_output}") and len(os.listdir(f"{project_output}"))>0:
            return max([int(f) for f in os.listdir(f"{project_output}")]) + 1
        else:
            return 1


    def set_intermediate_directory(self):
        
        if not os.path.exists(f"{self.__intermediate_dir}"):
            if self.dry_run:
                print (f"(Dry Run) {self.__intermediate_dir} created")
            else:
                os.mkdir(f"{self.__intermediate_dir}")

        return self.__intermediate_dir
    

    def set_project_folder(self):

        project_output = f"{self.__intermediate_dir}/{self.__project}"
        
        if not os.path.exists(f"{project_output}"):
            if self.dry_run:
                print (f"(Dry Run) {project_output} created")
            else:
                os.mkdir(f"{project_output}")

        return project_output
    

    def set_job_folder(self, number):

        job_folder = f"{self.__intermediate_dir}/{self.__project}/{number}"
        
        if not os.path.exists(f"{job_folder}"):
            if self.dry_run:
                print (f"(Dry Run) {job_folder} created")
            else:
                os.mkdir(f"{job_folder}")
        
        return job_folder



    def set_mutpred_output(self, number):

        output_num = f"{self.__intermediate_dir}/{self.__project}/{number}"
        if not os.path.exists(f"{output_num}"):
            if self.dry_run:
                print (f"(Dry Run) {output_num} created")
            else:
                os.mkdir(f"{output_num}")
        
        return output_num
    

    def set_log_output(self):
        logs = f"./logs"
        if not os.path.exists(f"{logs}"):
            if self.dry_run:
                print (f"(Dry Run) {logs} created")
            else:
                os.mkdir(f"{logs}")

        log_output = f"{logs}/{self.__project}"
        if not os.path.exists(f"{log_output}"):
            if self.dry_run:
                print (f"(Dry Run) {log_output} created")
            else:
                os.mkdir(f"{log_output}")



    def read_input_data(self):
        
        header_rows = []
        header_info = ""
        for i, line in enumerate(open(self.__input)):
            if line.startswith('##'):
                header_rows.append(i)
                header_info += line

        input_data = pd.read_csv(self.__input, sep="\t", skiprows = header_rows)

        return header_info, input_data



    def collect_value(self, x, val):
        cur_dict = {k.split("=")[0]:k.split("=")[1] for k in x.split(";")}
        
        if val in cur_dict.keys():
            protein = cur_dict[val]#.split(":")[0]
        else:
            protein = "."
        return protein



    def hgvsp(self, data):

        input_cols = set(data.columns)
        possibilites = set(["hgvsp", "HGVSp", "HGVSP"])
        inter = input_cols.intersection(possibilites)

        if "hgvsp" in input_cols:
            return data, True
        
        elif len(inter) == 1:
            data = data.rename(columns={inter[0]:"hgvsp"})
            return data, True
        
        elif "Extra" in input_cols and self.file_format=="VEP":
            
            try:
                data["hgvsp"] = data["Extra"].apply(lambda x: self.collect_value(x, val="HGVSp"))
                return data, True
            except KeyError:
                return data, False
            
        elif self.file_format == "dbNSFP":
            
            try:
                data["Ensembl_proteinid"] = data["Ensembl_proteinid"].str.split(";")
            except KeyError:
                print ("No Ensembl_proteinid column found")
                return data, False
            
            try:
                data["aapos"] = data["aapos"].astype(str)
                data["aapos"] = data["aapos"].str.split(";")
            except KeyError:
                print ("No aapos column found")
                return data, False
            
            try:
                data["HGVSp_VEP"] = data["HGVSp_VEP"].str.split(";")
            except KeyError:
                print ("No HGVSp_VEP column found")
                return data, False
            
            try:
                data["HGVSp_ANNOVAR"] = data["HGVSp_ANNOVAR"].str.split(";")
            except KeyError:
                print ("No HGVSp_ANNOVAR column found")
                return data, False
            
            try:
                data["VEP_canonical"] = data["VEP_canonical"].str.split(";")
            except KeyError:
                print ("No VEP_canonical column found")
                return data, False
            
            data = data[["#chr","pos(1-based)","ref","alt","aaref","aaalt","aapos","Ensembl_proteinid","HGVSp_VEP","HGVSp_ANNOVAR","VEP_canonical"]]
            #print (data)
            data = data.explode(["aapos","Ensembl_proteinid","HGVSp_VEP","HGVSp_ANNOVAR","VEP_canonical"])

            data["hgvsp"] = data.apply(lambda row: fasta.collect_dbNSFP_hgvsp(row["Ensembl_proteinid"], fasta.collect_dbNSFP_mutations(row)), axis=1)

            #print (data)
            #exit()
            
            return data, True
        
        elif self.file_format == "VCF-info" and self.annotation == "SnpEff":
            
            columns = []
            for line in self.header_data.split("\n"):
                if "ID=ANN" in line:
                    columns = line.split("'")[1].split(" | ")

            if len(columns) > 0:
                data["INFO"] = data["INFO"].apply(lambda x: self.collect_value(x, val="ANN")).str.split(",")
                data = data.explode("INFO")
                data[columns] = data["INFO"].str.split("|", expand=True)
                data = data[["#CHROM","POS","ID","REF","ALT","Annotation","Gene_Name","Gene_ID","Feature_Type","Feature_ID",'Transcript_BioType','HGVS.p']]
                
                transcript_to_protein_map = fasta.collect_fasta(self.__fasta_location, primary="Ensembl_transcriptid")[["Ensembl_transcriptid", "Ensembl_proteinid_v"]]
                
                data = data.merge(transcript_to_protein_map, left_on="Feature_ID", right_on="Ensembl_transcriptid", how="left")

                data["hgvsp"] = data["Ensembl_proteinid_v"] + ":" + data["HGVS.p"]
                
                return data, True
            
        else:
            return data, False
        

    
    def filter_canonical(self, data):

        if self.canonical and self.file_format == "VEP":
            data["CANONICAL"] = data["Extra"].apply(lambda x: self.collect_value(x, val="CANONICAL"))
            data = data[data["CANONICAL"]=="YES"]

            return data
        
        return data


    def cur_max_sequence_version(self, data):

        data["version"] = data["Ensembl_transcriptid"].str.split(".").str[-1]

        left_data = data.copy()

        data = data[["Ensembl_proteinid","version"]]
        data = pd.DataFrame(data.groupby("Ensembl_proteinid")["version"].max()).reset_index()
        data.rename(columns={"version": "latest_version"}, inplace=True)

        data = left_data.merge(data, on="Ensembl_proteinid", how="left")

        data = data.drop_duplicates(subset=["Ensembl_proteinid","Ensembl_transcriptid"])

        return data
    

    def max_sequence_version(self, data, fasta):

        seq_latest_versions = pd.DataFrame(fasta.groupby("Ensembl_proteinid")["version"].max()).reset_index()
        seq_latest_versions.rename(columns={"version": "latest_version"}, inplace=True)

        data = data.merge(seq_latest_versions, on="Ensembl_proteinid", how="left")

        return data
    
    
    def alt_sequences(self, data):
        
        mutations = data["mutation"]

        fasta_seqs = fasta.collect_fasta(pep_file=self.__fasta_location, drop_dups=False)

        alts = fasta_seqs[(fasta_seqs["Ensembl_geneid"].isin(data["Ensembl_geneid"])) | (fasta_seqs["Ensembl_transcriptid"].isin(data["Ensembl_transcriptid"]))].drop_duplicates()
        alts = alts[["Ensembl_proteinid","Ensembl_geneid","Ensembl_transcriptid","sequence","version"]]
        
        alts = alts.merge(mutations, how="cross")
        print  (f"Alts: {len(alts)}")

        if len(alts) > 0:
            
            alts["alignment_score"] = alts.apply(fasta.alignment_score, axis=1)
            alts = alts[alts["alignment_score"]==1]
            alts["num_mutations"] = alts["mutation"].str.split(" ").apply(len)
            alts["length"] = alts["sequence"].apply(len)
            alts = self.max_sequence_version(alts, fasta_seqs)
        
        
        return alts
        


    def add_sequences(self, data):

        if self.file_format == "dbNSFP":
            #print (len(set(data["Ensembl_proteinid"])))
            FF = self.fasta
            
            aligning = data.merge(FF, on="Ensembl_proteinid", how="left")
            
            aligning["alignment_score"] = aligning.apply(fasta.alignment_score, axis=1)

            aligning = self.max_sequence_version(fasta=self.fasta, data=aligning)

            alignment    = aligning[aligning["alignment_score"] == 1]

            non_aligning = aligning[(aligning["alignment_score"]<1) & (~aligning["Ensembl_proteinid"].isin(alignment["Ensembl_proteinid"]))]

            new_aligned = self.alt_sequences(non_aligning)

            data = pd.concat([alignment,new_aligned])
            
            aligned = data[(data["alignment_score"] == 1)].drop("latest_version", axis=1)

            aligned = self.cur_max_sequence_version(data=aligned)

            data = aligned[aligned["version"]==aligned["latest_version"]]
            #print (data)
        else:

            self.fasta["Ensembl_proteinid"] = self.fasta["Ensembl_proteinid_v"]
            self.fasta.drop(["Ensembl_proteinid_v", "version"], inplace=True, axis=1)
            data = data.merge(self.fasta, on="Ensembl_proteinid", how="left")
        
        #exit()
        return data



    def collect_mutations(self, data):

        if self.file_format == "dbNSFP":
            data["mutation"] = data.apply(lambda x: fasta.collect_dbNSFP_mutations(x), axis=1)

        else:
            #print (data)
            data[["Ensembl_proteinid","mutation"]] = data["hgvsp"].str.split(":",expand=True)
            data["mutation"] = data["mutation"].str.split(".").str[1].apply(lambda x: fasta.mutation_mapping(x))
        
        return data
    


    def groupby_protein_id(self, data):

        data = data[["Ensembl_proteinid","mutation"]].drop_duplicates()
        data = data.groupby(["Ensembl_proteinid"])["mutation"].apply(' '.join).reset_index()
        data["num_mutations"] = data["mutation"].str.split(" ").apply(lambda x: len(x))
        
        return data
    


    def count_mutations(self, data):
        num_mutations = sum(data["mutation"].str.split(" ").apply(len))
        return num_mutations


    def filtered_scored(self, data):

        if self.get_database_status():
            scores = pd.read_sql("SELECT seq_hash,mutation,'True' as scored FROM mutations", con=self.get_database_connection())
            scores["scored"] = scores["scored"].astype(bool)

            print ("Filtering scored variants.")
            print (f"Pre-filter: {len(data)}")

            cur_cols = data.columns
            
            data["seq_hash"] = data["sequence"].apply(lambda x: hashlib.md5(x.encode(encoding='utf-8')).hexdigest())
            data["mutation"] = data["mutation"].str.split(" ")
            data = data.explode("mutation")
            
            data = data.merge(scores, on=["seq_hash","mutation"], how="left")
            data["scored"].fillna(False, inplace=True)
            
            data = data[~data["scored"]]
            data = data.drop_duplicates()
            columns = [col for col in data.columns if col != "mutation"]
            data = pd.DataFrame(data.groupby(columns)["mutation"].apply(lambda x: ' '.join(list(x)))).reset_index()
            data = data[cur_cols]
            print (f"Post-filter: {len(data)}")

        else:
            scores = pd.read_csv("scores/MutPred2.tsv", sep="\t")[["hgvsp","MutPred2 score"]].drop_duplicates()

            cur_cols = data.columns

            data["mutation"] = data["mutation"].str.split(" ")
            
            data = data.explode("mutation")
            print ("FIX NEEDED")
            exit()

            data = data.merge(scores, on="hgvsp", how="left")

            data = data[data["MutPred2 score"].isna()]

            data = data.drop("MutPred2 score", axis=1)

            exclude = ['ENSP00000332973.4', 'ENSP00000355192.3', 'ENSP00000262186.5', 'ENSP00000305692.3', 'ENSP00000218516.4', 'ENSP00000352608.2', 'ENSP00000304408.4', 'ENSP00000039007.4', 'ENSP00000303992.5', 'ENSP00000347507.3', 'ENSP00000303208.5', 'ENSP00000261584.4', 'ENSP00000265849.7', 'ENSP00000333984.5', 'ENSP00000347942.3', 'ENSP00000300036.5', 'ENSP00000327145.8', 'ENSP00000355533.2', 'ENSP00000325527.5', 'ENSP00000287878.3', 'ENSP00000495254.2', 'ENSP00000369497.3', 'ENSP00000369129.3', 'ENSP00000234420.5', 'ENSP00000258439.3', 'ENSP00000257555.5', 'ENSP00000261448.5', 'ENSP00000295754.5', 'ENSP00000350283.3', 'ENSP00000298552.3', 'ENSP00000398266.2', 'ENSP00000280904.6', 'ENSP00000233242.1', 'ENSP00000290378.4', 'ENSP00000361107.2', 'ENSP00000155840.2', 'ENSP00000324856.6', 'ENSP00000415516.5', 'ENSP00000219476.3', 'ENSP00000233146.2', 'ENSP00000257430.4', 'ENSP00000373574.4', 'ENSP00000267163.4', 'ENSP00000261590.8', 'ENSP00000341838.5', 'ENSP00000454071.1', 'ENSP00000357283.4', 'ENSP00000342800.5', 'ENSP00000301761.3', 'ENSP00000341551.3', 'ENSP00000256474.3', 'ENSP00000442795.1', 'ENSP00000224784.6', 'ENSP00000407590.2', 'ENSP00000242839.5', 'ENSP00000231790.3', 'ENSP00000344666.5', 'ENSP00000394933.3', 'ENSP00000385107.4', 'ENSP00000364649.3', 'ENSP00000356953.3', 'ENSP00000362299.4', 'ENSP00000361021.3', 'ENSP00000364699.3', 'ENSP00000364133.4', 'ENSP00000269305.4', 'ENSP00000499593.1', 'ENSP00000228841.8', 'ENSP00000262340.5', 'ENSP00000351490.4', 'ENSP00000417404.1', 'ENSP00000292327.4']
            data["protein_id"] = data["hgvsp"].str.split(":").str[0]
            #print (data)
            print (f"Pre-filter: {len(data)}")
            data = data[~data["protein_id"].isin(exclude)]
            print (f"Post-filter: {len(data)}")
            data.drop("protein_id", inplace=True, axis=1)

        return data



    def sequence_quality_check(self, data):
        

        data["sequence"] = data["sequence"].apply(lambda x: fasta.clean_FASTA_sequence(x))
        data[["status","Sequence Errors", "Mutation Errors"]] = data.apply(lambda x: pd.Series(fasta.check_sequences(x)), axis=1)
        
        return data
    


    def calc_time_estimate(self, data):

        estimated_time = str(np.sum(data["Time Estimate (hrs)"].values))

        total_hours = int(estimated_time.split(".")[0]) + 1

        return total_hours



    def write_sequence_to_file(self, number, file_number, header, sequence):

        output_folder = self.set_job_folder(number=file_number)

        output_file = f"{output_folder}/{self.__base}.missense_{file_number}.faa"
        #print (f"Writing > {output_file}")
        if number == 0:
            #if os.path.exists(output_file):
            output = open(output_file, "w")
            output.write(header)
            output.write(sequence)
        else:
            output = open(output_file, "a")
            output.write(header)
            output.write(sequence)
        output.close()



    def clean_splits(self, row):
        row["num_mutations"] = len(row["mutation"].split(" "))
        row["Time Estimate (hrs)"] = row["Time per Mutation (hrs)"]*row["num_mutations"]
        
        return row



    def split_mutations(self, x, threshold, num_overflow_mutations):
        mutations = x["mutation"].split(" ")

        overflow_mutations = mutations[0:num_overflow_mutations]
        leftover_mutations = mutations[num_overflow_mutations:len(mutations)]

        size_of_splits = int(len(leftover_mutations)/((x["Time Estimate (hrs)"] - (num_overflow_mutations * x["Time per Mutation (hrs)"]))/threshold))
        
        if len(overflow_mutations) > 0:
            return [" ".join(overflow_mutations)] + [" ".join(leftover_mutations[i:i+size_of_splits]) for i in range(0,len(mutations[num_overflow_mutations:len(mutations)]),size_of_splits)]
        else:
            return [" ".join(leftover_mutations[i:i+size_of_splits]) for i in range(0,len(mutations[num_overflow_mutations:len(mutations)]),size_of_splits)]



    def split_sequence(self, data, cur_number, variable, threshold):
        
        overflow = max(0, cur_number - threshold)

        num_overflow_mutations = int(overflow / data["Time per Mutation (hrs)"])

        data = data.to_frame().T

        data["mutation"] = data.apply(lambda x: self.split_mutations(x, threshold, num_overflow_mutations), axis=1)
        
        data = data.explode("mutation")
        data = data.apply(lambda row: self.clean_splits(row), axis=1)

        data["Time Estimate (hrs)"] = data.apply(lambda row: len(row["mutation"].split(" "))*row["Time per Mutation (hrs)"], axis=1)

        return data
         


    def split_data(self, data, file_number=1):
        
        number = 0
        job_information = []
        memory = []
        threshold = self.get_time()
        variable  = "Time Estimate (hrs)"
        
        print (data)

        for index,row in data.iterrows():

            if number + row[variable] > threshold + 2:
                
                splits = self.split_sequence(row, number, variable, threshold)
                num_of_splits = len(splits)
                
                cur_split = 0
                for index,split in splits.iterrows():
                    cur_split += 1
                    header = f">{split['Ensembl_proteinid']}|{split['gene_symbol']} {split['mutation']}\n"
                    sequence = f"{split['sequence']}\n"
                    
                    if self.dry_run:
                        pass
                        #self.set_mutpred_output(file_number)
                    else:
                        self.write_sequence_to_file(number, file_number, header, sequence)
                        self.set_mutpred_output(file_number)

                    memory.append(split["Memory Estimate (MB)"])
                    number = split["Time Estimate (hrs)"]

                    if cur_split != num_of_splits:
                        job_information.append({
                            "File": file_number,
                            "Time Estimate": number,
                            "Memory Minimum": max(memory)
                        })
                        memory = []
                        file_number += 1
            else:
                header = f">{row['Ensembl_proteinid']}|{row['gene_symbol']} {row['mutation']}\n"
                sequence = f"{row['sequence']}\n"
                #print (row)
                memory.append(row["Memory Estimate (MB)"])
                ## write to file 
                if self.dry_run:
                    pass
                    #self.set_mutpred_output(file_number)

                else:
                    self.write_sequence_to_file(number, file_number, header, sequence)
                    self.set_mutpred_output(file_number)
                
                number += row[variable]
                print (number, threshold)
            
                if number >= threshold:

                    job_information.append({
                        "File": file_number,
                        "Time Estimate": number,
                        "Memory Minimum": max(memory)
                    })

                    number = 0
                    file_number += 1
                    memory = []
                else:
                    pass
        else:
            job_information.append({
                "File": file_number,
                "Time Estimate": number,
                "Memory Minimum": max(memory)
            })

            if self.dry_run:
                    pass
                    #self.set_mutpred_output(file_number)

            else:
                self.write_sequence_to_file(number, file_number, header, sequence)
                self.set_mutpred_output(file_number)
                
        job_information = pd.DataFrame(job_information)

        job_information["Normal Memory"] = job_information["Memory Minimum"].apply(lambda x: x <= 5000)
        job_information["Middle Memory"] = job_information["Memory Minimum"].apply(lambda x: x > 5000 and x <= 10000)
        job_information["High Memory"]   = job_information["Memory Minimum"].apply(lambda x: x > 10000)

        return job_information

    

    def prepare_mutpred_input(self):

        input_data = self.get_input_data()
        
        self.variant_data, hgvsp_status = self.hgvsp(input_data)
        
        if hgvsp_status:

            self.variant_data = self.filter_canonical(self.variant_data)

            self.variant_data = fasta.filter_non_missense(self.variant_data, self.file_format, self.annotation)
            
            self.variant_data = self.collect_mutations(self.variant_data)

            print (f"Pre filtered variants: {self.count_mutations(self.variant_data)}")
            
            self.variant_data = self.groupby_protein_id(self.variant_data)

            self.variant_data = self.add_sequences(self.variant_data)

            self.variant_data = self.filtered_scored(self.variant_data)
            print (self.variant_data)
            
            self.variant_data["Time Estimate (hrs)"] = self.variant_data.apply(lambda row: row["num_mutations"]*row["Time per Mutation (hrs)"], axis=1)

            self.variant_data = self.sequence_quality_check(self.variant_data)

            ## Check if any sequences are invalid
            not_passed = self.variant_data[self.variant_data["status"]==False]
            if len(not_passed) > 0:
                print (f"{len(not_passed)} proteins failed quality check")
                #print (not_passed)
                #exit()
                self.variant_data = self.variant_data[self.variant_data["status"]]

            ## Check if any of the proteins are unmapped
            unmapped = self.variant_data[pd.isna(self.variant_data["sequence"])]
            if len(unmapped) > 0:
                print ("======= UNMAPPED PROTEINS ======")
                print (unmapped)
                exit()
            
            #print (self.variant_data.sort_values("Time Estimate (hrs)"))
            #print (self.variant_data)

            #self.summary(self.variant_data)
            #print (self.variant_data)

            print (f"Post filtered variants: {self.count_mutations(self.variant_data)}")

            self.set_intermediate_directory()
            
            self.set_project_folder()

            self.set_log_output()

            tech_requirements = self.split_data(self.variant_data, file_number=self.__fasta_file_number_start)
            
            user = 1
            if user == 1:
                per_user_tech_requirements = lsf.split_for_multiple_users(tech_requirements, users=1)

                for user in range(len(per_user_tech_requirements)):
                    lsf.usage_report(per_user_tech_requirements[user])

                    lsf.build_lsf_config_file(per_user_tech_requirements[user], self.get_intermediate_dir(), self.get_project(), self.get_base(), user, self.dry_run)
            else:
                lsf.usage_report(tech_requirements)

                lsf.build_lsf_config_file(tech_requirements, self.get_intermediate_dir(), self.get_project(), self.get_base(), user, self.dry_run)

        else:
            print (f"HGVSp key not found in input columns -> [{','.join(self.variant_data.columns)}]")
            exit()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Take a list of chromosomal variants, Ensembl protein mutations, or other formats and output faa protein sequence files ready for MutPred2 intake.')

    def time_minimum(x):
        x = int(x)

        if x < 3:
            raise argparse.ArgumentTypeError("Minimum time is 5 hours")
        return x
    

    parser.add_argument('--time', type=time_minimum, nargs="?", default=24,
                        help="Target time in hours to run the jobs")
    parser.add_argument('--input', type=str, nargs=1,
                        help='The name of the input filename that is located in the data folder.')
    parser.add_argument('--project', type=str, nargs=1,
                        help='The name of the project for organization purposes')
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--canonical", action="store_true")

    
    args = parser.parse_args()

    input   = args.input[0]
    project = args.project[0]
    time    = args.time
    dry_run = args.dry_run
    canonical = args.canonical

    
    mut = MutPredpy(
        input=input,
        project=project,
        time=time,
        dry_run=dry_run,
        canonical=canonical
    )
    
    mut.prepare_mutpred_input()