import pandas as pd
import numpy as np

import os
import re
import logging
logger = logging.getLogger()

from ..fasta import fasta
from ..computing import lsf
from ..utils import utils as u

from .input_processing import process_input
from .jobs import split_data

class Prepare:
    def __init__(self, input, working_dir, time, dry_run, canonical, all_possible, verbose, users, fasta):
        
        u.set_logger_level(log_level = 20 if verbose else 25)

        ## Track dry run file and directory creation for logging purposes
        self.logging_status = {
            "working_dir":      False,
            "jobs_directory":   False,
            "job_folder":       dict(),
            "log_directory":    False,
            "script_directory": False,
            "input_faa_files":  dict()
        }

        self.__working_dir  = working_dir
        self.__input        = input
        self.__time         = time
        self.__fasta        = fasta
        
        self.dry_run        = dry_run
        self.all_possible   = all_possible
        self.users          = users
        

        self.__fasta_file_number_start = self.set_fasta_file_number_start()
        
        self.canonical = canonical

        if self.__input != "":
            self.__header_data, self.__input_data = self.read_input_data(self.get_input_path())
        

    def get_input_path(self):
        
        if os.path.exists(self.__input):
            path = self.__input

        elif os.path.exists(os.path.abspath(self.__input)):
            path = os.path.abspath(self.__input)
        
        else:
            logger.info(f"File input {self.__input} not found...trying {self.get_working_dir()}/{self.__input}")
            path = f"{self.get_working_dir()}/{self.__input}"
        
        if os.path.exists(path):
            return path
        else:
            raise FileNotFoundError(f"Input file {path} not found.")
    

    def get_input(self):
        return self.__input_data
    

    def get_input_headers(self):
        return self.__header_data
    
    def read_input_data(self, input_path):
        
        header_rows = []
        header_info = ""
        for i, line in enumerate(open(input_path)):
            if line.startswith('##'):
                header_rows.append(i)
                header_info += line
            else:
                break
        
        input_data = pd.read_csv(input_path, sep="\t", skiprows = header_rows)

        return header_info, input_data
    

    def get_working_dir(self):
        
        working_dir, self.logging_status["working_dir"] = u.create_directory(self.__working_dir, dry_run=self.dry_run, logged_status=self.logging_status.get("working_dir"))
        
        return os.path.abspath(working_dir).rstrip("/")
    

    def get_base(self):
        return ".".join(str(self.get_input_path()).split("/")[-1].split(".")[:-1])
    
    def get_time(self):
        return self.__time
    

    def set_fasta_file_number_start(self):

        jobs_dir = f"{self.get_jobs_directory()}"

        def job_num(job):
            if isinstance(job, int):
                return int(job)
            elif isinstance(job, str) and re.search(r"job_\d+",job):
                return int(job.split("_")[-1])
            elif isinstance(job, str) and re.search(r'\d+',job):
                return int(re.search(r'\d+',job).group())
            else:
                return 0

        if not os.path.exists(f"{jobs_dir}"):
            return 1
        elif os.path.exists(f"{jobs_dir}") and len(os.listdir(f"{jobs_dir}"))>0:
            return max([job_num(f) for f in os.listdir(f"{jobs_dir}")]) + 1
        else:
            return 1
    

    def get_jobs_directory(self):

        jobs_directory, self.logging_status["jobs_directory"] = u.create_directory(f"{self.get_working_dir()}/jobs", dry_run=self.dry_run, logged_status=self.logging_status.get("jobs_directory"))
        return jobs_directory
        

    def get_job_folder(self, number):

        job_folder, self.logging_status["job_folder"][number] = u.create_directory(f"{self.get_jobs_directory()}/{number}", dry_run=self.dry_run, logged_status=self.logging_status.get("job_folder").get(number))
        return job_folder
    

    def get_log_folder(self):

        logs, self.logging_status["log_folder"] = u.create_directory(f"{self.get_working_dir()}/logs", dry_run=self.dry_run, logged_status=self.logging_status.get("log_folder"))
        return logs
    
        

    def mapping_unversioned_sequence(self, FF, data, col_mapping):
        
        mutation_col = col_mapping['mutation_column']
        id_cols = col_mapping['id_column']

        # Split version from IDs in FF
        for col in id_cols:
            FF[f'{col}_unversioned'] = FF[col].str.split('.').str[0]
            FF[f'{col}_version'] = FF[col].str.split('.').str[1].astype(float)

        # Split version from IDs in data
        for col in id_cols:
            data[f'{col}_unversioned'] = data[col].str.split('.').str[0]
            #data[f'{col}_version'] = data[col].str.split('.').str[1].astype(float)

        
        # Merge substitutions into FF using unversioned IDs
        FF = pd.merge(
            FF,
            data[[f'{col}_unversioned' for col in id_cols] + [mutation_col]],
            on=[f'{col}_unversioned' for col in id_cols],
            how='right'
        )
        

        # Function to check reference amino acids match
        def check_reference_match(row):

            if self.all_possible:
                substitutions = row[mutation_col]
            else:
                substitutions = row[mutation_col].split()
            sequence = row['sequence']

            if isinstance(sequence, float) and np.isnan(sequence):
                return False
            
            elif substitutions == None and self.all_possible:
                return True

            for sub in substitutions:
                ref_aa, pos, alt_aa = sub[0], int(sub[1:-1]), sub[-1]
                if pos <= len(sequence) and sequence[pos - 1] != ref_aa:
                    return False
            return True


        # Filter sequences that match reference amino acids
        matching_FF = FF[FF.apply(check_reference_match, axis=1)]
        
        # Find highest version per unversioned ID in matching sequences
        highest_version_FF = (
            matching_FF.groupby([f'{col}_unversioned' for col in id_cols] + ["sequence", 'Memory Estimate (MB)', 'Time per Mutation (hrs)']).agg({
                f"{col}_version":"max" for col in id_cols
            })
            .reset_index(drop=False)
        )

        # Merge data with the highest version entries from matching sequences in FF
        merged_df = pd.merge(
            data,
            highest_version_FF,
            on=[f'{col}_unversioned' for col in id_cols],
            how='left'
        )
        

        # Filter best matching sequences with reference match
        best_matches = merged_df[merged_df['sequence'].notna() & merged_df.apply(check_reference_match, axis=1)]
        

        # Versionize the non-versioned input IDs
        for col in id_cols:
            best_matches.loc[:, col] = best_matches.apply(lambda row: f"{row[col]}.{int(row[f'{col}_version'])}", axis=1)


        # Select relevant columns for output
        output = best_matches[
            id_cols + [mutation_col, 'sequence', 'Memory Estimate (MB)', 'Time per Mutation (hrs)']
        ]

        return output


    def add_sequences(self, data, validation_results):
        
        def get_protein_ids(fasta):
            
            ## First try to collect only versioned protein and transcript IDs
            protein_ids = [p for p in validation_results.keys() if p in fasta.columns and validation_results[p]["versioned"] and validation_results[p]["found"]]
            if protein_ids:
                return protein_ids, True
            
            ## If no versioned protein or transcript IDs are found, collect all unversioned IDs
            protein_ids = [p for p in validation_results.keys() if p in fasta.columns and validation_results[p]["found"]]
            if protein_ids:
                return protein_ids, False
            
            raise Exception(f"Input columns {', '.join([p for p in validation_results.keys() if validation_results[p]['found']])} not found in FASTA file")

        
        col_mapping = {"mutation_column": "Substitution"}

        ## Collect FASTA files
        if not self.__fasta:
            logger.info(f"--fasta flag not used, defaulting to Ensembl FASTA files (versions 90-110)")
            Ensembl_FF_grch37 = fasta.collect_ensembl_fasta(assembly="GRCh37")
            Ensembl_FF_grch38 = fasta.collect_ensembl_fasta(assembly="GRCh38")
            
            col_mapping["id_column"], versioned = get_protein_ids(Ensembl_FF_grch38)
            
            if not self.all_possible:
                data = self.groupby_id(data, col_mapping=col_mapping)
                pre_mapped_numbers = len(data)
            else:
                pre_mapped_numbers = len(data)

            if versioned:
            
                data_grch37 = data.merge(Ensembl_FF_grch37, on=col_mapping.get("id_columns"), how="left")
                unmapped_grch37 = len(data_grch37[data_grch37["sequence"].isna()])

                data_grch38 = data.merge(Ensembl_FF_grch38, on=col_mapping.get("id_columns"), how="left")
                unmapped_grch38 = len(data_grch38[data_grch38["sequence"].isna()])

            else:
                logger.info(f"Unversioned IDs used in input file. Finding the best matches in the FASTA file.")

                data_grch37 = self.mapping_unversioned_sequence(Ensembl_FF_grch37, data, col_mapping)
                unmapped_grch37 = len(data_grch37[data_grch37["sequence"].isna()])

                data_grch38 = self.mapping_unversioned_sequence(Ensembl_FF_grch38, data, col_mapping)
                unmapped_grch38 = len(data_grch38[data_grch38["sequence"].isna()])

            # Return the sequence mapped dataframe with the least number of unmapped IDs
            data = data_grch38 if unmapped_grch38 <= unmapped_grch37 else data_grch37
        else:
            logger.info(f"Using FASTA file {self.__fasta}")
            fasta_seqs = fasta.read_fasta(self.__fasta)
            col_mapping["id_column"], versioned = get_protein_ids(fasta_seqs)

            if not self.all_possible:
                data = self.groupby_id(data, col_mapping=col_mapping)
                pre_mapped_numbers = len(data)
            else:
                pre_mapped_numbers = len(data)

            if versioned:
                data = data.merge(fasta_seqs, on=col_mapping["id_column"], how="left")
            else:
                logger.info(f"Unversioned IDs used in input file. Finding the best matches in the FASTA file.")
                data = self.mapping_unversioned_sequence(fasta_seqs, data, col_mapping)

        
        data = data.filter(items=col_mapping['id_column'] + [col_mapping['mutation_column'], "sequence", "Memory Estimate (MB)", "Time per Mutation (hrs)"])
        logger.info(f"Mapped to protein/transcript sequences using columns {', '.join(col_mapping['id_column'])}")

        ## Check if any of the proteins are unmapped
        data = self.check_unmapped(data, col_mapping, pre_mapped_numbers)

        ## Check that everything is deduplicated
        data = data.drop_duplicates()

        ## Check that the sequences conform to MutPred2 standards
        logger.info(f"Running sequence quality checks")
        data = fasta.sequence_quality_check(data)

        ## Collect all possible amino acid substitutions from mapped sequence if --all-possible flag is used
        if self.all_possible:
            logger.info(f"--all-possible selected. Collecting all possible amino acid substitutions.")
            data = u.AminoAcidMap.get_all_possible_mutations(data, col_mapping)

        ## Check that the sequences and mutations are concordant with MutPred2 expectations
        logger.info(f"Running data quality checks of mutations")
        data = fasta.data_quality_check(self, data, col_mapping)

        return data, col_mapping


    def check_unmapped(self, data, col_mapping, pre_mapped_numbers):
        
        if pre_mapped_numbers > len(data):
            log_error_message = f"{pre_mapped_numbers - len(data)} proteins dropped during sequence mapping." 
            logger.warning(f"{log_error_message}")

        unmapped = data[pd.isna(data["sequence"])]
        if len(unmapped) > 0:

            ## Log missing proteins/transcripts
            log_error_message = f"{len(unmapped)} proteins have been dropped. No matching protein sequences found." 
            logger.warning(f"{log_error_message} (additional info in {self.get_log_folder()}/errors.log)")

            if self.dry_run:
                logger.dry_run(f"Would've writen logs to {self.get_log_folder()}/errors.log)")
            else:
                u.missing_sequence_error_logging(unmapped, self.get_log_folder(), col_mapping, "No matching protein sequences found.")
            
            ## Drop unmapped proteins and mutations
            data = data[~pd.isna(data["sequence"])]
        
        if data.empty:
            raise Exception("Could not map proteins to sequences. Try a different FASTA file?")

        return data


    def groupby_id(self, data, col_mapping):

        id_col = col_mapping["id_column"]
        muts   = col_mapping["mutation_column"]

        cols = id_col + [muts]

        data = data[cols].drop_duplicates()
        data = data.groupby(id_col)[muts].apply(' '.join).reset_index()
        data[f"num_{muts}"] = data[muts].str.split(" ").apply(lambda x: len(x))
        
        return data


    def filtered_scored(self, data):

        """
        TODO: Create function to filter already scored missense variants. This will need to be built in conjunction with catalog.
        """
        
        return data
    


    def calc_time_estimate(self, data):

        estimated_time = str(np.sum(data["Time Estimate (hrs)"].values))

        total_hours = int(estimated_time.split(".")[0]) + 1

        return total_hours



    def write_sequence_to_file(self, number, file_number, header, sequence):

        job_folder = self.get_job_folder(number=file_number)

        faa_input_file = f"{job_folder}/input.faa"
        
        if self.dry_run:
            if not self.logging_status.get("input_faa_files").get(file_number):
                logger.dry_run(f"Would've written job {file_number} input.faa")
                self.logging_status["input_faa_files"][file_number] = True

        else:
            if number == 0:
                output = open(faa_input_file, "w")
                output.write(header)
                output.write(sequence)
            else:
                output = open(faa_input_file, "a")
                output.write(header)
                output.write(sequence)

            output.close()
    

    def prepare_mutpred_input(self):
        
        ## Read in input data
        input_data = self.get_input()
        #print (input_data)
        #input_data = input_data.head(2000)

        ## Process input data for use with the MutPred Suite
        logger.info(f"Processing input data {self.get_input_path()}")
        self.variant_data, validation_results = process_input(input_data, canonical=self.canonical, all_possible=self.all_possible)

        ## Add protein/transcript sequences
        logger.info(f"Mapping protein/transcript sequences")
        self.variant_data, col_mapping = self.add_sequences(self.variant_data, validation_results)
        
        #from random import sample, randint

        #self.variant_data["Substitution"] = self.variant_data["Substitution"].apply(lambda x: sample(x.split(" "), randint(1, 15)))
        #self.variant_data["Substitution"] = self.variant_data["Substitution"].apply(lambda x: x.split(" "))
        #self.variant_data = self.variant_data.explode("Substitution")
        #self.variant_data["Ensembl_protein_id"] = self.variant_data["ENSP"].str.split(".").str[0]
        #self.variant_data["Ensembl_transcript_id"] = self.variant_data["ENST"].str.split(".").str[0]
        #self.variant_data["Ensembl_gene_id"] = self.variant_data["ENSG"].str.split(".").str[0]
        #print (self.variant_data)
        
        #self.variant_data[["Ensembl_protein_id","Ensembl_transcript_id","Ensembl_gene_id","Substitution"]].sample(n=15000).to_csv("/Users/bergqt01/Research/MutPredPy/MutPredPy/test/test_data/test.ENSP_ENST_ENSG.txt", sep="\t", index=False)
        

        # TODO: filter out already scored mutations
        #self.variant_data = self.filtered_scored(self.variant_data)

        ## Estimate time to complete MutPred2 computations on a per protein level
        logger.info(f"Estimating time to completion")
        self.variant_data["Time Estimate (hrs)"] = self.variant_data.apply(lambda row: len(row[col_mapping["mutation_column"]].split(" "))*row["Time per Mutation (hrs)"], axis=1)

        ## Calculate the computational resources and number of parallel jobs necessary to complete the MutPred runs in the designated time.
        logger.info(f"Calculating resource requirements and parallelizing jobs")
        tech_requirements = split_data(self, self.variant_data, col_mapping, file_number=self.__fasta_file_number_start)
        
        ## Divide MutPred2 jobs evenly amoung the designated users     
        per_user_tech_requirements = lsf.split_for_multiple_users(tech_requirements, users=self.users)
        
        for user in range(len(per_user_tech_requirements)):
            
            ## Generate the usage report if --verbose is active
            lsf.usage_report(per_user_tech_requirements[user])

            ## Build and output the 
            lsf.build_lsf_config_file(self, per_user_tech_requirements[user], user + 1)
       

if __name__ == "__main__":
    pass