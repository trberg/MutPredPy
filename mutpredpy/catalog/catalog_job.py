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

        self.__num_files = len(
            [
                position_file
                for position_file in os.listdir(job_path)
                if re.search(r"positions_pu_\d+.mat", position_file)
            ]
        )
        self.__file_sizes = [
            sio.loadmat(position_file).get("positions_pu").shape[0]
            for position_file in [
                f"{job_path}/output.txt.positions_pu_{i+1}.mat"
                for i in range(self.__num_files)
            ]
        ]
        self._positions = self.get_positions()

    def get_job_path(self):
        """
        Retrieves the job directory path.

        Returns:
            str: The absolute path of the job directory.
        """
        return self.__job_path

    def get_num_files(self):
        """
        Retrieves the number of files for each output file.

        Returns:
            int: The number of files per output type.
        """
        return self.__num_files

    def get_file_sizes(self):
        """
        Retrieves the size of each file output type.

        Returns:
            list: The size of each file output type.
        """
        return self.__file_sizes

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
        mutations = os.path.join(self.get_job_path(), "output.txt.substitutions.mat")
        muts_df = sio.loadmat(f"{mutations}").get("substitutions")
        mutations = [
            str(muts[0])
            for mutation_list in muts_df.flatten()
            for muts in mutation_list.flatten()
        ]
        return mutations

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
            add_feature_weights(feature_array, feature_column_info)
            for feature_array in all_features
        ]

        return output_features

    def get_sequences(self):
        """
        Extracts protein sequences from the MutPred2 output file.

        Returns:
            list: A list of sequence strings extracted from the job's output.
        """
        sequences = sio.loadmat(
            os.path.join(self.get_job_path(), "output.txt.sequences.mat")
        ).get("sequences")
        sequence_strings = [
            str(item)
            for sublist in sequences
            for arr in sublist
            for item in np.atleast_1d(arr)
        ]

        return sequence_strings

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

        sequences = self.get_sequences()

        position_files = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.positions_pu_{i+1}.mat")
            )
            for i in range(self.get_num_files())
        ]

        positions = np.array(
            [
                np.array([f"{sequences[i][pos - 1]}{pos}" for pos in dim] + ["-"])
                for i, pos_files in enumerate(position_files)
                for dim in pos_files["positions_pu"]
            ]
        )
        return positions

    def get_motifs(self):
        """
        Extracts motif-related information from MutPred2 output files.

        Returns:
            list: A list of motif strings extracted from the job's output.
        """
        motifs = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.motif_info_{i+1}.mat")
            ).get("motif_info")
            for i in range(self.get_num_files())
        ]

        motif_strings = [
            f"Altered Motifs ({str(motif_string[0])})"
            for motif_file in motifs
            for motif_arry in motif_file
            for motif_string in motif_arry
        ]

        return motif_strings

    def get_notes(self):
        """
        Extracts and maps notes related to mutations from MutPred2 output files.

        Returns:
            list: A list of formatted remark strings corresponding to mutation-related notes.
        """
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
            if self.__catalog.features:
                mut_features.to_csv(features_path, sep="\t", index=False)

            if self.__catalog.mechanisms:
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
        items = [(feature, details) for feature, details in data_dict.items()]

        # Sort by the P-value in ascending order
        sorted_items = sorted(
            items, key=lambda x: x[1]["Posterior Probability"], reverse=True
        )

        return sorted_items

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
        sorted_mechanisms = self.sort_by_p_value(mechanisms)

        for mech, info in sorted_mechanisms:
            if info["P-value"] < 1 and info["Posterior Probability"] >= 0.01:
                mech = self.updated_property_list.get(mech, mech)

                if mech == "Motifs":
                    pass
                    # print(f"Altered {mech}")
                elif mech in self.Region and mech in self.Altered:
                    mech_formatted = f"Altered {mech} (Pr = {round(info['Posterior Probability'], 2)}| P = {round(info['P-value'], 2)})"  # pylint: disable=C0301
                    formatted_mechanisms.append(mech_formatted)
                elif mech in self.Altered:
                    mech_formatted = f"Altered {mech} at {info['Effected Position']} (Pr = {round(info['Posterior Probability'], 2)}| P = {round(info['P-value'], 2)})"  # pylint: disable=C0301
                    formatted_mechanisms.append(mech_formatted)
                elif mech in self.Region:
                    mech_formatted = f"{info['Type']} of {mech} (Pr = {round(info['Posterior Probability'], 2)} | {round(info['P-value'], 2)})"  # pylint: disable=C0301
                    formatted_mechanisms.append(mech_formatted)
                else:
                    mech_formatted = f"{info['Type']} of {mech} at {info['Effected Position']} (Pr = {round(info['Posterior Probability'], 2)} | {round(info['P-value'], 2)})"  # pylint: disable=C0301
                    formatted_mechanisms.append(mech_formatted)

        data = {
            "ID": prot_id,
            "Substitution": mutation,
            "MutPred2 score": score,
            "Molecular mechanisms with Pr >= 0.01 and P < 1.00": "; ".join(
                formatted_mechanisms
            ),
            "Motif information": motif,
            "Remarks": remarks,
        }

        return data

    def process_job(self, catalog):
        """
        Processes a MutPred2 job by extracting and cataloging mutation-related data.

        Args:
            catalog (Catalog): An instance of the Catalog class to store job data.

        Returns:
            pd.DataFrame: A DataFrame containing extracted mutation information,
                        including substitutions, mechanisms, motifs, and remarks.
        """
        # Retrieve protein IDs associated with the job
        protein_ids = self.get_ids()

        # Retrieve unique sequence hashes for the input sequences
        sequence_hashes = self.get_sequence_hash()

        # Extract mutations (substitutions) identified in the job
        substitutions = self.get_mutations()

        # Retrieve MutPred2 scores for the identified mutations
        mutpred2_scores = self.get_mutpred2_scores()

        # Collect mechanisms if cataloging mechanisms is enabled
        if self.__catalog.mechanisms:
            mechanisms = Mechanisms.collect_mechanisms(catalog, self)
        else:
            # If mechanisms are not collected, fill with None values
            mechanisms = [None for i in range(substitutions)]

        # Extract motif information related to the mutations
        motifs = self.get_motifs()

        # Retrieve additional remarks or notes related to the mutations
        remarks = self.get_notes()

        # Collect feature data if cataloging features is enabled
        if self.__catalog.features:
            features = self.get_features()
        else:
            # If features are not collected, fill with None values
            features = [None for i in range(substitutions)]

        # Initialize a list to store MutPred2 output data
        mutpred2_output = []

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
