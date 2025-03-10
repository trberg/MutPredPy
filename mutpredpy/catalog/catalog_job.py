"""
Catalog Job Module for MutPredPy

This module handles job-specific cataloging operations for MutPred2, 
including extracting mutations, mechanisms, motifs, features, remarks 
and sequence positions from job output files.
"""

import os
import re
import logging
import yaml
import scipy.io as sio
import numpy as np
import pandas as pd

from mutpredpy.utils import utils as u
from mutpredpy.fasta import fasta
from .mechanisms import Mechanisms

logger = logging.getLogger()


class CatalogJob:
    """
    Handles job-specific cataloging operations for MutPred2.

    This class processes MutPred2 job outputs, extracting relevant data
    including mutations, mechanisms, motifs, and sequence positions.

    Attributes:
        job_path (str): The path to the job directory.
        catalog (Catalog): Instantiated Catalog object.
    """

    remarks_list = [
        "Sequence length less than 30 residues",
        "Substitution in incorrect format",
        "Position lies outside the sequence",
        "Reference residue does not match the one at the position in the sequence",
        "Reference residue and substituted residue are the same",
        "Non-standard amino acids found within 10 residues",
        "No PSSM",
        "No conservation scores",
        "Predicted conservation scores",
    ]
    updated_property_list = {
        "Intracellular": "Cytoplasmic",
        "Extracellular": "Non-cytoplasmic",
        "VSL2B_disorder": "Intrinsic_disorder",
        "Surface_accessibility": "Relative_solvent_accessibility",
        "GPI_anchor_amidation": "GPI-anchor_amidation",
        "B_factor": "B-factor",
        "N_terminal_signal": "N-terminal_signal",
        "C_terminal_signal": "C-terminal_signal",
    }
    To_skip = "Non_transmembrane"
    Altered = [
        "N_terminal_signal",
        "Signal_helix",
        "C_terminal_signal",
        "Signal_cleavage",
        "Cytoplasmic_loop",
        "Transmembrane_region",
        "Non_cytoplasmic_loop",
        "Non_transmembrane",
        "Coiled_coil",
        "Calmodulin_binding",
        "DNA_binding",
        "RNA_binding",
        "PPI_residue",
        "PPI_hotspot",
        "MoRF",
        "Stability",
    ]
    Region = Altered + [
        "Helix",
        "Strand",
        "Loop",
        "Intrinsic_disorder",
        "B_factor",
        "Relative_solvent_accessibility",
    ]

    def __init__(self, job_path, catalog):

        self.__job_path = job_path
        self.__catalog = catalog
        self._mutations = self.get_mutations()
        self._num_mutations = len(self._mutations)
        self._sequences = self.get_sequences()

        self.__num_files, self.__file_sizes = self.get_num_and_size_files()

        self._positions = self.get_positions()

    def get_job_path(self):
        """
        Retrieves the job directory path.

        Returns:
            str: The absolute path of the job directory.
        """
        return self.__job_path

    def get_num_and_size_files(self):
        """
        Retrieves the maximum number of files for any output type and calculates file sizes.

        Returns:
            tuple: (num_files, file_sizes), where:
                - num_files (int): The max number of files found per output type.
                - file_sizes (list): A list containing the sizes of the most frequent file type.
        """

        if hasattr(self, "__num_files") and hasattr(self, "__file_sizes"):
            return self.__num_files, self.__file_sizes

        compiled_patterns = [
            re.compile(pattern)
            for pattern in [
                r"output\.txt\.(?P<key>feats)_\d+\.mat",
                r"output\.txt\.(?P<key>models)_\d+\.mat",
                r"output\.txt\.(?P<key>motif_info)_\d+\.mat",
                r"output\.txt\.(?P<key>MutPred2Score)_\d+\.mat",
                r"output\.txt\.(?P<key>notes)_\d+\.mat",
                r"output\.txt\.(?P<key>positions)_\d+\.mat",
                r"output\.txt\.(?P<key>positions_pu)_\d+\.mat",
                r"output\.txt\.(?P<key>propX)_\d+\.mat",
            ]
        ]

        # Initialize a counter for each pattern
        file_counts = [0] * len(compiled_patterns)
        file_groups = [[] for _ in range(len(compiled_patterns))]
        key_map = {}  # Store mapping of files to their extracted keys

        job_path = self.get_job_path()

        # Scan files once, updating the respective counters
        for job_file in os.listdir(job_path):
            for i, pattern in enumerate(compiled_patterns):
                match = pattern.search(job_file)
                if match:
                    file_counts[i] += 1
                    full_path = os.path.join(job_path, job_file)
                    file_groups[i].append(full_path)
                    key_map[full_path] = match.group("key")

        # Determine the file type with the maximum occurrences
        max_num_files = max(file_counts)
        max_index = file_counts.index(
            max_num_files
        )  # Get the index of the most frequent file type

        # Extract the corresponding file group (files with max occurrences)
        max_files = file_groups[max_index][
            :max_num_files
        ]  # Ensure we take only the expected number

        # Load file sizes using the dynamically extracted key
        file_sizes = []
        for file in max_files:
            key = key_map.get(file)  # Get the key from the stored mapping
            mat_data = sio.loadmat(file)
            if key in mat_data:  # Ensure the key exists in the loaded file
                file_sizes.append(mat_data[key].shape[0])

        # Store the max count
        self.__num_files = max(file_counts)
        self.__file_sizes = file_sizes

        return self.__num_files, self.__file_sizes

    def get_file_sizes(self):
        """
        Retrieves the sizes of the relevant output files in the job directory.

        Returns:
            dict: A dictionary containing file names as keys and their sizes as values.
        """
        return self.__file_sizes

    def get_num_files(self):
        """
        Retrieves the total number of output files per type related to the job.

        Returns:
            int: The number of relevant output files in the job directory.
        """
        return self.__num_files

    def get_ids(self):
        """
        Extracts and returns the list of protein or transcript IDs from the job's input FASTA file.

        Returns:
            list: A list of IDs corresponding to the sequences in the input FASTA file.
        """
        input_file = fasta.read_mutpred_input_fasta(
            os.path.join(self.get_job_path(), "input.faa")
        )
        input_file = input_file.filter(["ID", "Substitution"])
        input_file["Substitution"] = input_file["Substitution"].str.split(",")
        input_file = input_file.explode("Substitution")
        input_file = input_file["ID"].to_list()
        return input_file

    def get_mutations(self):
        """
        Extracts mutation data from the job's MutPred2 output file.

        Returns:
            np.ndarray: An array of mutation identifiers.
        """

        if hasattr(self, "_mutations"):
            return self._mutations

        mutations = os.path.join(self.get_job_path(), "output.txt.substitutions.mat")
        if os.path.exists(mutations):
            muts_df = sio.loadmat(f"{mutations}").get("substitutions")
            mutations = [
                str(muts[0])
                for mutation_list in muts_df.flatten()
                for muts in mutation_list.flatten()
            ]
        else:
            original_info = fasta.read_mutpred_input_fasta(
                os.path.join(self.get_job_path(), "input.faa")
            )
            original_info = original_info["Substitution"].str.split(",")
            mutations = original_info.explode().to_list()

        self._mutations = mutations

        return self._mutations

    def get_mutpred2_scores(self):
        """
        Extracts MutPred2 scores from job output files.

        Returns:
            np.ndarray: An array containing all extracted MutPred2 scores.
        """
        if os.path.exists(
            os.path.join(self.get_job_path(), "output.txt.MutPred2Score_1.mat")
        ):
            mutpred2_score_files = [
                sio.loadmat(
                    os.path.join(
                        self.get_job_path(), f"output.txt.MutPred2Score_{i+1}.mat"
                    )
                ).get("S")
                for i in range(self.get_num_files())
            ]

            all_mutpred2_scores = [
                float(mutpred2_score)
                for mutpred2_score_file in mutpred2_score_files
                for mutpred2_score_array in mutpred2_score_file
                for mutpred2_score in mutpred2_score_array
            ]
        else:
            mutpred2_score_file = pd.read_csv(
                os.path.join(self.get_job_path(), "output.txt"),
                sep=",",
                usecols=["MutPred2 score"],
            )
            all_mutpred2_scores = mutpred2_score_file["MutPred2 score"].to_list()

        return all_mutpred2_scores

    def get_features(self):
        """
        Placeholder method for extracting feature-related data from MutPred2 output files.

        Returns:
            None
        """
        if not os.path.exists(
            os.path.join(self.get_job_path(), "output.txt.feats_1.mat")
        ):
            return [pd.DataFrame() for i in range(self._num_mutations)]

        features = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.feats_{i+1}.mat")
            ).get("feats")
            for i in range(self.get_num_files())
        ]

        all_features = np.array(
            [feat_arry for feat_file in features for feat_arry in feat_file]
        )

        feature_column_info = self.__catalog.load_feature_columns()

        def add_feature_weights(feature_weights, feature_columns):
            feature_columns["Feature weight"] = feature_weights
            return feature_columns

        output_features = [
            add_feature_weights(feature_array, feature_column_info.copy())
            for feature_array in all_features
        ]

        return output_features

    def get_sequences(self):
        """
        Extracts protein sequences from the MutPred2 output file.

        Returns:
            list: A list of sequence strings extracted from the job's output.
        """

        if hasattr(self, "_sequences"):
            return self._sequences

        sequence_file = os.path.join(self.get_job_path(), "output.txt.sequences.mat")
        if os.path.exists(sequence_file):
            sequences = sio.loadmat(os.path.join(sequence_file)).get("sequences")
            sequence_strings = [
                str(item)
                for sublist in sequences
                for arr in sublist
                for item in np.atleast_1d(arr)
            ]
        else:
            original_info = fasta.read_mutpred_input_fasta(
                os.path.join(self.get_job_path(), "input.faa")
            )
            original_info["Substitution"] = original_info["Substitution"].str.split(",")
            original_info = original_info.explode("Substitution")
            sequence_strings = original_info["sequence"].to_list()

        self._sequences = sequence_strings

        return self._sequences

    def get_sequence_hash(self):
        """
        Extracts protein sequences from the MutPred2 output file.

        Returns:
            list: A list of sequence strings extracted from the job's output
                and zipped out so that it's the full job length.
        """
        sequences = self.get_sequences()

        all_sequences = np.array(
            [
                u.get_seq_hash(sequences[i])
                for i, file_size in enumerate(self.get_file_sizes())
                for _ in range(file_size)
            ]
        )

        return all_sequences

    def get_positions(self):
        """
        Extracts position data for mutations from MutPred2 output files.

        Returns:
            np.ndarray: A nested array containing formatted position information
                        for each sequence.
        """
        if hasattr(self, "_positions"):
            return self._positions

        position_file = os.path.join(
            self.get_job_path(), "output.txt.positions_pu_1.mat"
        )
        if not os.path.exists(position_file):
            return []

        sequences = self.get_sequences()

        position_files = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.positions_pu_{i+1}.mat")
            )
            for i in range(self.get_num_files())
        ]

        positions = [
            np.array([f"{sequences[i][pos - 1]}{pos}" for pos in dim] + ["-"])
            for i, pos_files in enumerate(position_files)
            for dim in pos_files["positions_pu"]
        ]
        return positions

    def get_motifs(self):
        """
        Extracts motif-related information from MutPred2 output files.

        Returns:
            list: A list of motif strings extracted from the job's output.
        """
        motif_file = os.path.join(self.get_job_path(), "output.txt.motif_info_1.mat")
        if os.path.exists(motif_file):
            motifs = [
                sio.loadmat(
                    os.path.join(
                        self.get_job_path(), f"output.txt.motif_info_{i+1}.mat"
                    )
                ).get("motif_info")
                for i in range(self.get_num_files())
            ]

            motif_strings = [
                f"Altered Motifs ({str(motif_string[0])})"
                for motif_file in motifs
                for motif_arry in motif_file
                for motif_string in motif_arry
            ]
        else:
            motif_strings = ["-" for i in range(self._num_mutations)]

        return motif_strings

    def get_notes(self):
        """
        Extracts and maps notes related to mutations from MutPred2 output files.

        Returns:
            list: A list of formatted remark strings corresponding to mutation-related notes.
        """
        note_file = os.path.join(self.get_job_path(), "output.txt.notes_1.mat")
        if os.path.exists(note_file):
            note_indices = np.vstack(
                [
                    sio.loadmat(
                        os.path.join(self.get_job_path(), f"output.txt.notes_{i+1}.mat")
                    ).get("notes")
                    for i in range(self.get_num_files())
                ]
            )

            def map_indices_to_strings(index_lists):
                remarks = [self.remarks_list[i] for i in index_lists if i == 1]
                if len(remarks) == 0:
                    return "-"
                return ". ".join(remarks)

            remarks = [map_indices_to_strings(note) for note in note_indices]
        else:
            remarks = ["-" for sub in range(self._num_mutations)]

        return remarks

    def write_mutation_to_catalog(
        self,
        mutation,
        mutpred2_score,
        sequence_hash,
        mut_mechanisms,
        mut_motifs,
        mut_remarks,
        mut_features,
    ):
        """
        Writes mutation data to the catalog directory in YAML format.

        Args:
            row (pd.Series): A Pandas Series containing mutation details including
                            substitution, MutPred2 score, and mechanisms.

        Returns:
            None
        """

        pos = re.search(r"\d+", mutation).group()
        _, alt = mutation.split(pos)

        mutpred_data = {
            "Substitution": mutation,
            "MutPred2 Score": mutpred2_score,
            "Motifs": mut_motifs,
            "Notes": mut_remarks,
        }

        mutation_path = os.path.join(
            self.__catalog.catalog_location, "scores", sequence_hash, pos, alt
        )

        if self.__catalog.dry_run:
            logger.dry_run("Would've created catalog directory %s", mutation_path)

        else:
            os.makedirs(mutation_path, exist_ok=True)

        output_yaml_path = os.path.join(mutation_path, "output.yaml")
        mechanisms_path = os.path.join(mutation_path, "mechanisms.csv")
        features_path = os.path.join(mutation_path, "features.csv")

        if self.__catalog.dry_run:
            logger.dry_run("Would've created the file %s", output_yaml_path)
            logger.dry_run("Would've created the file %s", mechanisms_path)
            logger.dry_run("Would've created the file %s", features_path)

        else:
            with open(output_yaml_path, "wt", encoding="utf-8") as yaml_file:
                yaml.dump(
                    mutpred_data, yaml_file, default_flow_style=False, sort_keys=False
                )
            if self.__catalog.features and not mut_features.empty:
                mut_features.to_csv(features_path, sep="\t", index=False)

            if self.__catalog.mechanisms and self.get_positions():
                mut_mechanisms.to_csv(mechanisms_path, sep="\t", index=False)

    def sort_by_p_value(self, data_dict):
        """
        Sorts the dictionary entries based on the P-value from lowest to highest.

        Args:
            data_dict (dict): Dictionary containing mutation effect data.

        Returns:
            list: A list of tuples sorted by P-value, where each tuple consists of
                (feature_name, feature_data).
        """
        # Convert dictionary to a list of tuples with (feature_name, feature_data)
        items = list(data_dict.items())

        # Sort by the P-value in ascending order
        sorted_items = sorted(
            items, key=lambda x: x[1]["Posterior Probability"], reverse=True
        )

        return sorted_items

    def string_format_mechanism(self, cur_mechanism):
        """
        Formats molecular mechanism data into a readable string representation.

        Args:
            cur_mechanism (pd.DataFrame): DataFrame containing mechanism-related data with
                                        columns for property, probability, p-value,
                                        affected position, and mutation type.

        Returns:
            str: A formatted string summarizing the mechanisms with relevant probabilities
                and affected positions.
        """
        cur_mechanism = cur_mechanism[
            (cur_mechanism["P-value"] < 1)
            & (cur_mechanism["Posterior Probability"] >= 0.01)
            & (cur_mechanism["Property"] != "Motifs")
        ]
        cur_mechanism = cur_mechanism.sort_values("P-value")

        formatted_mechanisms = []

        def mech_format(row):
            mech = row["Property"]
            post_pr = row["Posterior Probability"]
            pval = row["P-value"]
            eff_pos = row["Effected Position"]
            mut_type = row["Type"]

            if mech in self.Region and mech in self.Altered:
                mech_formatted = f"Altered {mech} (Pr = {round(post_pr, 2)}| P = {round(pval, 2)})"  # pylint: disable=C0301
                return mech_formatted
            if mech in self.Altered:
                mech_formatted = f"Altered {mech} at {eff_pos} (Pr = {round(post_pr, 2)}| P = {round(pval, 2)})"  # pylint: disable=C0301
                return mech_formatted
            if mech in self.Region:
                mech_formatted = f"{mut_type} of {mech} (Pr = {round(post_pr, 2)} | P = {round(pval, 2)})"  # pylint: disable=C0301
                return mech_formatted

            mech_formatted = f"{mut_type} of {mech} at {eff_pos} (Pr = {round(post_pr, 2)} | P = {round(pval, 2)})"  # pylint: disable=C0301
            return mech_formatted

        formatted_mechanisms = [mech_format(row) for i, row in cur_mechanism.iterrows()]

        formatted_mechanisms = "; ".join(formatted_mechanisms)

        return formatted_mechanisms

    def standard_mutpred2_output_format(
        self,
        prot_id,
        mutation,
        score,
        mechanisms,
        motif,
        remarks,
    ):
        """
        Formats extracted MutPred2 mutation data into a standardized output format.

        Args:
            mutation (str): Mutation identifier.
            score (float): MutPred2 score associated with the mutation.
            mechanisms (dict): Dictionary containing molecular mechanisms and their probabilities.
            motif (str): Motif information related to the mutation.
            remarks (str): Additional remarks or notes on the mutation.

        Returns:
            dict: A dictionary containing the formatted MutPred2 output data.
        """

        formatted_mechanisms = []
        sorted_mechanisms = self.sort_by_p_value(mechanisms.sort_values(""))

        for _, info in sorted_mechanisms:
            formatted_mechanisms = [self.string_format_mechanism(info)]

        data = {
            "ID": prot_id,
            "Substitution": mutation,
            "MutPred2 score": score,
            "Molecular mechanisms with Pr >= 0.01 and P < 1.00": formatted_mechanisms,
            "Motif information": motif,
            "Remarks": remarks,
        }

        return data

    def already_processed(self, sequence_hashes, substitutions):
        """
        Checks if the given mutations have already been processed and stored in the catalog.

        Args:
            sequence_hashes (list): List of sequence hash values.
            substitutions (list): List of mutation substitutions.

        Returns:
            bool: True if all mutations have been processed and stored, otherwise False.
        """

        for i, mutation in enumerate(substitutions):

            pos = re.search(r"\d+", mutation).group()
            _, alt = mutation.split(pos)

            try:
                mutation_path = os.path.join(
                    self.__catalog.catalog_location,
                    "scores",
                    sequence_hashes[i],
                    pos,
                    alt,
                )
            except IndexError:
                print(self.get_job_path())
                print(substitutions)
                print(sequence_hashes)
                exit()
            scores_path = os.path.join(mutation_path, "output.yaml")
            mechs_path = os.path.join(mutation_path, "mechanisms.csv")
            feats_path = os.path.join(mutation_path, "features.csv")

            scores_status = os.path.exists(scores_path)

            if self.__catalog.mechanisms:
                mech_status = os.path.exists(mechs_path)
            else:
                mech_status = True

            if self.__catalog.features:
                feats_status = os.path.exists(feats_path)
            else:
                feats_status = True

            if not (scores_status and mech_status and feats_status):
                return False

        return True

    def process_job(self):
        """
        Processes a MutPred2 job by extracting and cataloging mutation-related data.

        Args:
            catalog (Catalog): An instance of the Catalog class to store job data.

        Returns:
            pd.DataFrame: A DataFrame containing extracted mutation information,
                        including substitutions, mechanisms, motifs, and remarks.
        """

        # Retrieve MutPred2 scores for the identified mutations
        mutpred2_scores = self.get_mutpred2_scores()

        # If no scores found, assume this is an incomplete job and end function
        if not mutpred2_scores:
            return

        # Retrieve unique sequence hashes for the input sequences
        sequence_hashes = self.get_sequence_hash()

        # Extract mutations (substitutions) identified in the job
        substitutions = self.get_mutations()

        if self.already_processed(sequence_hashes, substitutions):
            return

        # Collect mechanisms if cataloging mechanisms is enabled
        if self.__catalog.mechanisms and self.get_positions():
            mechanisms = Mechanisms.collect_mechanisms(self.__catalog, self)
            # string_formated_mechanisms = [
            #    self.string_format_mechanism(mech) for mech in mechanisms
            # ]
        else:
            # If mechanisms are not collected, fill with None values
            mechanisms = [None for i in substitutions]
            # string_formated_mechanisms = ["-" for i in substitutions]

        # Extract motif information related to the mutations
        motifs = self.get_motifs()

        # Retrieve additional remarks or notes related to the mutations
        remarks = self.get_notes()

        # Collect feature data if cataloging features is enabled
        if self.__catalog.features:
            features = self.get_features()
        else:
            # If features are not collected, fill with None values
            features = [pd.DataFrame() for sub in substitutions]

        # Iterate over each mutation and store relevant data in the catalog
        for i, sub in enumerate(substitutions):
            self.write_mutation_to_catalog(
                sub,  # Substitution mutation
                mutpred2_scores[i],  # MutPred2 score for the mutation
                sequence_hashes[i],  # Hash of the sequence where mutation occurs
                mechanisms[i],  # Mechanistic insights related to the mutation
                motifs[i],  # Motif-related information
                remarks[i],  # Remarks or notes related to the mutation
                features[i],  # Additional extracted features (if available)
            )

        # print(
        #     pd.DataFrame(
        #         {
        #             "Substitution": substitutions,
        #             "MutPred2 Scores": mutpred2_scores,
        #             "Mechanisms": string_formated_mechanisms,
        #             "Motifs": motifs,
        #             "Remarks": remarks,
        #         }
        #     )
        # )
        # exit()
