"""
Catalog Module for MutPredPy

This module handles job cataloging, including checking directories, 
loading features, and collecting mechanisms from MutPred2 outputs.
"""

import os
import re
import importlib.resources as pkg_resources
import yaml
import numpy as np
import pandas as pd
import scipy.io as sio
import logging
from numba import njit, prange

from ..utils import utils as u
from ..fasta import fasta, Protein


class Catalog:
    """Handles job cataloging for MutPred2 results."""

    def __init__(self, job_dir, mechanisms, dry_run):
        """
        Initializes the Catalog class.

        Args:
            job_dir (str): Path to job directory.
            mechanisms (bool): Whether to include mechanisms.
            dry_run (bool): Run without saving changes.
        """
        u.set_logger_level(log_level = 20)
        self.job_dir = self.check_jobs_directory(job_dir)
        self.catalog_location = u.catalog_directory()

        self.dry_run = dry_run
        self.mechanisms = mechanisms

        if self.mechanisms:
            neutral_properties = sio.loadmat(
                pkg_resources.files("mutpredpy").joinpath("pkg_resources/pu_features_null_distributions.mat")
            )
            self.neutral_property_cols = [
                neuts[0].replace("-","_") for neuts in neutral_properties.get("neut_properties").flatten()
            ]
            self.property_names = [
                p.replace("_gain","") for p in self.neutral_property_cols if p.endswith("_gain") or p == "Stability"
            ]
            neutral_property_gain_cols = np.array([self.neutral_property_cols.index(col) for col in self.neutral_property_cols if col.endswith("_gain") or col == "Stability"])
            neutral_property_loss_cols = np.array([self.neutral_property_cols.index(col) for col in self.neutral_property_cols if col.endswith("_loss") or col == "Stability"])

            # Load null distributions
            self.null_distributions_array = neutral_properties.get("neut_ntrans_feats")
            self.null_gain_distributions_array = self.null_distributions_array[:,neutral_property_gain_cols]
            self.null_loss_distributions_array = self.null_distributions_array[:,neutral_property_loss_cols]
            
            self.n = self.null_distributions_array.shape[0]
            self.null_distributions = pd.DataFrame(data=self.null_distributions_array, columns=self.neutral_property_cols)

    
    def check_jobs_directory(self, job_dir):
        """
        Validates the given job directory and its contents.

        Args:
            job_dir (str): Path to job directory.

        Returns:
            str: Validated job directory.

        Raises:
            FileNotFoundError: If the job directory does not exist.
        """
        
        discrepancies = []
        job_pattern = re.compile(r"^(job_\d+|\d+)$")  # Pattern for "job_<number>" or just a number.

        # Check if the job directory exists
        if not os.path.exists(job_dir):
            raise FileNotFoundError(f"Job directory {job_dir} not found")
        
        # Iterate through items in the job directory
        for job in os.listdir(job_dir):
            job_path = os.path.join(job_dir, job)

            # Check if it's a valid job folder
            if os.path.isdir(job_path) and job_pattern.match(job):

                # Check required files inside the job folder
                required_files = ["input.faa", "output.txt"]
                missing_files = [file for file in required_files if not os.path.isfile(os.path.join(job_path, file))]

                if missing_files:
                    discrepancies.append(f"Missing files in job '{job}': {', '.join(missing_files)}\n")
            else:
                # If it's not a valid job folder, flag it
                if os.path.isdir(job_path):
                    discrepancies.append(f"Unexpected folder: {job}\n")
                else:
                    discrepancies.append(f"Unexpected file: {job}\n")

        if len(discrepancies) > 0:
            raise Exception(f"Issues found in job directory {job_path}.\n{''.join(discrepancies)}")
        
        return job_dir


    def mechanisms_to_dict(self, row, value_name):
        
        try:
            cur_row = row.drop("Substitution")
        except KeyError:
            cur_row = row
        
        properties = {i:{value_name: r} for i,r in cur_row.items()}

        return properties
    

    def mechanisms_np_to_dict(self, row, value_name):
        """
        Convert a NumPy array of property values into a dictionary.

        Args:
            row (np.ndarray): Array of property values.
            value_name (str): Name for the value field in the dictionary.

        Returns:
            dict: Dictionary with property names as keys and value_name as sub-key.
        """
        # Ensure the row is a NumPy array for fast indexing
        if not isinstance(row, np.ndarray):
            row = np.array(row)

        # Create a dictionary where each property maps to its corresponding value
        properties = {col: {value_name: row[i]} for i, col in enumerate(self.neutral_property_cols)}

        return properties
    

    def combine_mechs(self, row):
        """Combines two dictionaries of mechanisms."""
        cur_dict = {}
        for k in row["Mechanisms_pr"].keys():
            cur_dict[k] = dict()
            cur_dict[k].update(row["Mechanisms_pr"][k])
            cur_dict[k].update(row["Mechanisms_p"][k]) 

        return cur_dict
    

    def get_property_scores_and_pvalues(self, job, mutations):
        ## Collect predicted property scores
        prop_score_files = [prop_file for prop_file in os.listdir(job) if re.search(r"output.txt.prop_scores_pu_\d+.mat", prop_file)]
        labeled_property_scores = pd.concat(
            [pd.DataFrame(data=sio.loadmat(os.path.join(job, f"output.txt.prop_scores_pu_{p+1}.mat")).get("prop_scores_pu"), 
                          columns=self.property_names) for p in range(len(prop_score_files))]).reset_index(drop=True)
        labeled_property_scores["Substitution"] = mutations
        labeled_property_scores["Mechanisms"]    = labeled_property_scores.apply(lambda row: self.mechanisms_to_dict(row, "Posterior Probability"), axis=1)
        labeled_property_scores = labeled_property_scores.filter(["Substitution","Mechanisms"])

        ## Collect predicted property p-values
        prop_pvalues_files = [prop_file for prop_file in os.listdir(job) if re.search(r"output.txt.prop_pvals_pu_\d+.mat", prop_file)]
        labeled_property_pvalues = pd.concat(
            [pd.DataFrame(data=sio.loadmat(os.path.join(job, f"output.txt.prop_pvals_pu_{p+1}.mat")).get("prop_pvals_pu"), 
                          columns=self.property_names) for p in range(len(prop_pvalues_files))]).reset_index(drop=True)
        labeled_property_pvalues["Substitution"] = mutations
        labeled_property_pvalues["Mechanisms"]    = labeled_property_pvalues.apply(lambda row: self.mechanisms_to_dict(row, "P-value"), axis=1)
        labeled_property_pvalues = labeled_property_pvalues.filter(["Substitution","Mechanisms"])
        

        ## Combine scores and p-values
        labeled_properies = labeled_property_scores.merge(labeled_property_pvalues, on="Substitution", how="inner", suffixes=["_pr","_p"])
        

        ## Format combination as a combined dictionary
        labeled_properies["Mechanisms"] = labeled_properies.apply(self.combine_mechs, axis=1)
        labeled_properies = labeled_properies.filter(["Substitution", "Mechanisms"])
        
        return labeled_properies
    


    def get_properties_from_propX(self, job, mutations):
        
        output_mechanisms = pd.DataFrame()
        output_mechanisms["Substitution"] = mutations
        
        propx_pu_files = [propX for propX in os.listdir(job) if re.search(r"output.txt.propX_pu_\d+.mat", propX)]
        
        predicted_properties = pd.concat([pd.DataFrame(data=sio.loadmat(os.path.join(job, f"output.txt.propX_pu_{p+1}.mat")).get("propX_pu"), columns=self.neutral_property_cols) for p in range(len(propx_pu_files))]).reset_index(drop=True)

        
        loss_columns = [col for col in predicted_properties.columns if re.search(r'_loss$|Stability', col)]
        loss_props = predicted_properties.filter(loss_columns)


        gain_columns = [col for col in predicted_properties.columns if re.search(r'_gain$|Stability', col)]
        gain_props = predicted_properties.filter(gain_columns)


        def p_value(x, column):
            pval = np.sum(self.null_distributions[column] >= x) / self.n
            return pval


        property_scores = pd.DataFrame()
        property_pvalues = pd.DataFrame()
        prop_columns = [col.replace("_loss","") for col in loss_columns]
        
        ## 
        for i, row in loss_props.iterrows():
            for col in prop_columns:
    
                if col == "Stability":
                    prop_loss_score = row[col]
                    prop_gain_score = gain_props.loc[i,col]
                    property_pvalues.loc[i, col] = p_value(prop_gain_score, col) if prop_gain_score > prop_loss_score else p_value(prop_loss_score, col)
                
                else:
                    prop_loss_score = row[f"{col}_loss"]
                    prop_gain_score = gain_props.loc[i,f"{col}_gain"]
                    try:
                        property_pvalues.loc[i, col] = p_value(prop_gain_score, f"{col}_gain") if prop_gain_score > prop_loss_score else p_value(prop_loss_score, f"{col}_loss")
                    except ValueError:
                        print (prop_loss_score)
                        print (prop_gain_score)
                        exit()
                
                property_scores.loc[i, col] = prop_gain_score if prop_gain_score > prop_loss_score else prop_loss_score

        
        output_mechanisms["Mechanisms_pr"] = property_scores.apply(lambda row: self.mechanisms_to_dict(row, "Posterior Probability"), axis=1)
        output_mechanisms["Mechanisms_p"] = property_pvalues.apply(lambda row: self.mechanisms_to_dict(row, "P-value"), axis=1)

        output_mechanisms["Mechanisms"] = output_mechanisms.apply(lambda row: self.combine_mechs(row), axis=1)
        output_mechanisms = output_mechanisms.filter(["Substitution", "Mechanisms"])

        return output_mechanisms
    

    # JIT-compiled p-value computation (Numba for ultra-fast execution)
    @njit(parallel=True)
    def compute_p_values(self, loss_props, gain_props, null_gain_distributions, null_loss_distributions, N):
        """
        Compute p-values using parallelized Numba function.

        Args:
            loss_props (np.ndarray): Array of loss property values.
            gain_props (np.ndarray): Array of gain property values.
            null_gain_distributions (np.ndarray): Null distribution for gain properties.
            null_loss_distributions (np.ndarray): Null distribution for loss properties.
            N (int): Number of samples in null distribution.

        Returns:
            np.ndarray: Computed p-values.
        """

        rows = int(loss_props.shape[0])
        cols = int(loss_props.shape[1])

        p_values = np.empty(loss_props.shape)
        for i in prange(rows):
            for j in prange(cols):
                null_gain_col = null_gain_distributions[:, j]  # Get column as array
                null_loss_col = null_loss_distributions[:, j]  # Get column as array
                if gain_props[i, j] > loss_props[i, j]:
                    p_values[i, j] = np.sum(null_gain_col >= gain_props[i, j]) / N
                else:
                    p_values[i, j] = np.sum(null_loss_col >= loss_props[i, j]) / N
        return p_values


    def vectorized_get_properties_from_propx(self, job, mutations):
        """
        Extracts and processes property scores and p-values for given mutations.

        Args:
            job (str): Path to job directory.
            mutations (list): List of mutation identifiers.

        Returns:
            pd.DataFrame: Dataframe containing processed mutation data.
        """
        output_mechanisms = pd.DataFrame()
        output_mechanisms["Substitution"] = mutations
        
        # Load all propX files
        propx_pu_files = [
            os.path.join(job, f) for f in os.listdir(job) if re.search(r"output.txt.propX_pu_\d+.mat", f)
        ]

        # Load and concatenate all matrices
        prop_data = [sio.loadmat(f)["propX_pu"] for f in propx_pu_files]
        predicted_properties = pd.DataFrame(
            np.vstack(prop_data), columns=self.neutral_property_cols
        )

        # Identify loss/gain columns efficiently
        loss_mask = predicted_properties.columns.str.contains(r'_loss$|Stability')
        loss_columns = predicted_properties.columns[loss_mask]
        gain_columns = loss_columns.str.replace("_loss", "_gain")

        # Convert to NumPy for fast operations
        loss_props = predicted_properties[loss_columns].to_numpy()
        gain_props = predicted_properties[gain_columns].to_numpy()

        # Compute max score across loss and gain
        max_scores = np.maximum(loss_props, gain_props)
        
        p_values = self.compute_p_values(loss_props, gain_props, self.null_gain_distributions_array, self.null_loss_distributions_array, self.n)

        # Convert NumPy back to Pandas
        property_scores = pd.DataFrame(max_scores, columns=[col.replace("_loss", "") for col in loss_columns])
        property_pvalues = pd.DataFrame(p_values, columns=[col.replace("_loss", "") for col in loss_columns])

        
        output_mechanisms["Mechanisms_pr"] = property_scores.apply(lambda row: self.mechanisms_to_dict(row, "Posterior Probability"), axis=1)
        output_mechanisms["Mechanisms_p"] = property_pvalues.apply(lambda row: self.mechanisms_to_dict(row, "P-value"), axis=1)

        output_mechanisms["Mechanisms"] = output_mechanisms.apply(self.combine_mechs, axis=1)
        output_mechanisms = output_mechanisms.filter(["Substitution", "Mechanisms"])

        return output_mechanisms
    

    def collect_mechanisms(self, job):
        
        mutations = os.path.join(job, "output.txt.substitutions.mat")
        muts_df  = sio.loadmat(f"{mutations}").get("substitutions")
        mutations = np.array([muts[0] for mutation_list in muts_df.flatten() for muts in mutation_list.flatten()])

        propX_pu  = os.path.join(job, "output.txt.propX_pu_1.mat")
        property_scores  = os.path.join(job, "output.txt.prop_scores_pu_1.mat")
        property_pvalues = os.path.join(job, "output.txt.prop_pvals_pu_1.mat")

        if self.mechanisms:

            if os.path.exists(property_scores) and os.path.exists(property_pvalues):
                
                properties = self.get_property_scores_and_pvalues(job, mutations)
                #propX = self.vectorized_get_properties_from_propX(job, mutations)
                
                #propX = self.get_properties_from_propX(job, mutations)
                #print ("Got propX properties scores")

                #def evaluate(row):
                #    for k,v in row["Mechanisms_x"].items():
                #        if not row["Mechanisms_x"][k] == row["Mechanisms_y"][k]:
                #            print (k, row["Mechanisms_x"][k], row["Mechanisms_y"][k])
                #            exit()
                    

                #testing = properties.merge(propX, on="Substitution", how="inner")
                #testing.apply(lambda row: evaluate(row), axis=1)
                #print (testing)
                #exit()

                return properties


            elif os.path.exists(propX_pu):

                properties = self.vectorized_get_properties_from_propx(job, mutations)

                return properties
                

        else:
            return pd.DataFrame()


    def write_mutation_to_catalog(self, row):
        
        mutation = row["Substitution"]
        pos = re.search(r"\d+",mutation).group()
        ref,alt = mutation.split(pos)

        if self.mechanisms:
            mutpred_data = {
                "Substitution": row["Substitution"],
                "MutPred2 Score": row["MutPred2 score"],
                "Mechanisms": row["Mechanisms"]
            }
        else:
            mutpred_data = {
                "Substitution": row["Substitution"],
                "MutPred2 Score": row["MutPred2 score"]
            }

        mutation_path = os.path.join(self.catalog_location, "scores", row["sequence_hash"], pos, alt)
        if self.dry_run:
            logging.info(f"Would've created catalog directory {mutation_path}")
            pass
        else:
            os.makedirs(mutation_path, exist_ok=True)

        mutpred_output_yaml_path = os.path.join(mutation_path, "mutpred2_output.yaml")
        if self.dry_run:
            logging.info(f"Would've created the file {mutpred_output_yaml_path}")
            pass
        else:
            with open(mutpred_output_yaml_path, "wt") as yaml_file:
                yaml.dump(mutpred_data, yaml_file, default_flow_style=False, sort_keys=False)
        

    def print_progress(self, cur_job, number_of_jobs, end="\r"):
        cur_percentage = cur_job/number_of_jobs
        cur_progress = int(50*cur_percentage)
        remaining_progress = 50 - cur_progress
        print (f" [{'=' * cur_progress}{' ' * remaining_progress}] {round((cur_percentage)*100, 1)}%", end=end)


    def catalog_jobs(self):
        
        catalog_index = []

        job_dirs = os.listdir(self.job_dir)
        #job_dirs = ["1","2","69","263"]
        number_of_jobs = len(job_dirs)
        cur_job = 0
        self.print_progress(cur_job, number_of_jobs)
        
        for job in job_dirs:
            
            
            output_path = os.path.join(self.job_dir, job, "output.txt")
            input_path = os.path.join(self.job_dir, job, "input.faa")

            output = pd.read_csv(output_path)
            input_faa = fasta.read_mutpred_input_fasta(input_path)
            input_faa["Substitution"] = input_faa["Substitution"].apply(lambda x: x.split(","))
            input_faa = input_faa.explode("Substitution")

            sequenced_data = output.merge(input_faa, on=["ID","Substitution"], how="left")


            if self.mechanisms and "Molecular mechanisms with Pr >= 0.01 and P < 1.00" in output.columns:
                keep_cols = ["ID","Substitution","MutPred2 score","Molecular mechanisms with Pr >= 0.01 and P < 1.00","sequence"]
            else:
                keep_cols = ["ID","Substitution","MutPred2 score","sequence"]


            sequenced_data = sequenced_data.filter(keep_cols)

            mechanisms = self.collect_mechanisms(os.path.join(self.job_dir, job))
            
            if not mechanisms.empty:

                sequenced_data = sequenced_data.merge(mechanisms, on="Substitution", how="left")
            
            else:

                sequenced_data["Mechanisms"] = dict()


            sequenced_data["sequence_hash"] = sequenced_data["sequence"].apply(u.get_seq_hash)

            sequenced_data.apply(lambda row: self.write_mutation_to_catalog(row), axis=1)

            cur_index = Protein.split_mutpred_output_ids(sequenced_data[["ID","sequence_hash"]].drop_duplicates())

            catalog_index.append(cur_index)

            cur_job += 1
            self.print_progress(cur_job, number_of_jobs)

        else:
            self.print_progress(cur_job, number_of_jobs, end="\n")
        
        catalog_index = pd.concat(catalog_index)
        
        return catalog_index.drop_duplicates()
    


if __name__=="__main__":

    pass

