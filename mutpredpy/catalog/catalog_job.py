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
                np.array([f"{sequences[i][pos - 1]}{pos}" for pos in dim] + ["A0"])
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
            str(motif_string[0])
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
        }
        if self.__catalog.mechanisms:
            mutpred_data["Mechanisms"] = mut_mechanisms

        mutpred_data["Motifs"] = mut_motifs
        mutpred_data["Notes"] = mut_remarks

        mutation_path = os.path.join(
            self.__catalog.catalog_location, "scores", sequence_hash, pos, alt
        )

        if self.__catalog.dry_run:
            logger.dry_run("Would've created catalog directory %s", mutation_path)

        else:
            os.makedirs(mutation_path, exist_ok=True)

        output_yaml_path = os.path.join(mutation_path, "output.yaml")
        mechanisms_yaml_path = os.path.join(mutation_path, "mechanisms.yaml")
        features_output = os.path.join(mutation_path, "features.csv")

        if self.__catalog.dry_run:
            logger.dry_run("Would've created the file %s", output_yaml_path)
            logger.dry_run("Would've created the file %s", mechanisms_yaml_path)
            logger.dry_run("Would've created the file %s", features_output)

        else:
            with open(output_yaml_path, "wt", encoding="utf-8") as yaml_file:
                yaml.dump(
                    mutpred_data, yaml_file, default_flow_style=False, sort_keys=False
                )
            if self.__catalog.features:
                mut_features.to_csv(features_output, sep="\t", index=False)

    # def standard_mutpred2_output_format(mutation, score, ):
    ##ID,Substitution,MutPred2 score,Molecular mechanisms with Pr >= 0.01 and P < 1.00,Motif information,Remarks

    def process_job(self, catalog):
        """
        Processes a MutPred2 job by extracting and cataloging mutation-related data.

        Args:
            catalog (Catalog): An instance of the Catalog class to store job data.

        Returns:
            pd.DataFrame: A DataFrame containing extracted mutation information,
                        including substitutions, mechanisms, motifs, and remarks.
        """

        sequence_hashes = self.get_sequence_hash()
        substitutions = self.get_mutations()

        mutpred2_scores = self.get_mutpred2_scores()

        if self.__catalog.mechanisms:
            mechanisms = Mechanisms.collect_mechanisms(catalog, self)
        else:
            mechanisms = [None for i in range(substitutions)]

        motifs = self.get_motifs()
        remarks = self.get_notes()

        if self.__catalog.features:
            features = self.get_features()
        else:
            features = [None for i in range(substitutions)]

        for i, sub in enumerate(substitutions):
            self.write_mutation_to_catalog(
                sub,
                mutpred2_scores[i],
                sequence_hashes[i],
                mechanisms[i],
                motifs[i],
                remarks[i],
                features[i],
            )

        return None
