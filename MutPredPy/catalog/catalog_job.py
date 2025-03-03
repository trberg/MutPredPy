import os
import re
import scipy.io as sio
import numpy as np
import pandas as pd

from .mechanisms import Mechanisms


class Catalog_Job:

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

        super().__init__(catalog.job_dir, catalog.mechanisms, catalog.dry_run)

        self.__job_path = job_path
        self.__num_files = len(
            [
                position_file
                for position_file in os.listdir(job_path)
                if re.search(r"positions_pu_\d+.mat", position_file)
            ]
        )

    def get_job_path(self):
        return self.__job_path

    def get_mutations(self):

        mutations = os.path.join(self.get_job_path(), "output.txt.substitutions.mat")
        muts_df = sio.loadmat(f"{mutations}").get("substitutions")
        mutations = np.array(
            [
                muts[0]
                for mutation_list in muts_df.flatten()
                for muts in mutation_list.flatten()
            ]
        )
        return mutations

    def get_features(self):
        pass

    def get_mechanism_types(self):
        pass

    def get_sequences(self):
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

    def get_positions(self):

        if hasattr(self, "_positions"):
            return self._positions

        sequences = self.get_sequences()

        position_files = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.positions_pu_{i+1}.mat")
            )
            for i in range(self.__num_files)
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

        motifs = [
            sio.loadmat(
                os.path.join(self.get_job_path(), f"output.txt.motif_info_{i+1}.mat")
            ).get("motif_info")
            for i in range(self.__num_files)
        ]

        motif_strings = [
            motif_string[0]
            for motif_file in motifs
            for motif_arry in motif_file
            for motif_string in motif_arry
        ]

        return motif_strings

    def get_notes(self):

        note_indices = np.vstack(
            [
                sio.loadmat(
                    os.path.join(self.get_job_path(), f"output.txt.notes_{i+1}.mat")
                ).get("notes")
                for i in range(self.__num_files)
            ]
        )

        def map_indices_to_strings(index_lists):
            remarks = [self.remarks_list[i] for i in index_lists if i == 1]
            if len(remarks) == 0:
                return "-"
            return ". ".join(remarks)

        remarks = [map_indices_to_strings(note) for note in note_indices]

        return remarks

    def process_job(self, catalog):

        #
        # exit()

        substitutions = self.get_mutations()
        print(substitutions.shape)
        mechanisms = Mechanisms.collect_mechanisms(catalog, self)
        print(len(mechanisms))
        motifs = self.get_motifs()
        print(len(motifs))
        remarks = self.get_notes()
        print(len(remarks))
        exit()
        job_catalog_info = pd.DataFrame(
            {
                "Substitution": substitutions,
                "Mechanisms": mechanisms,
                "Motifs": motifs,
                "Remarks": remarks,
            }
        )
        print(job_catalog_info)
        exit()
        return job_catalog_info
