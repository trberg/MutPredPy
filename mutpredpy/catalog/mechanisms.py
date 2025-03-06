"""
Mechanisms Module for MutPredPy

This module handles the extraction, processing, and analysis of molecular mechanisms 
from MutPred2 output files, including computing p-values and integrating mechanism data.
"""

import os
import re
from operator import itemgetter
import scipy.io as sio
import numpy as np
import pandas as pd
from numba import njit, prange


# JIT-compiled p-value computation (Numba for ultra-fast execution)
@njit(parallel=True)
def compute_p_values(
    loss_props,
    gain_props,
    null_gain_distributions,
    null_loss_distributions,
    n: int,
):
    """
    Compute p-values using parallelized Numba function.

    Args:
        loss_props (np.ndarray): Array of loss property values.
        gain_props (np.ndarray): Array of gain property values.
        null_gain_distributions (np.ndarray): Null distribution for gain properties.
        null_loss_distributions (np.ndarray): Null distribution for loss properties.
        n (int): Number of samples in null distribution.

    Returns:
        np.ndarray: Computed p-values.
    """

    rows = int(loss_props.shape[0])
    cols = int(loss_props.shape[1])

    # Ensure rows and cols are positive integers
    assert rows > 0 and cols > 0, "Shape of loss_props must be valid integers"

    p_values = np.empty((rows, cols), dtype=np.float64)
    mech_gains = np.empty((rows, cols), dtype=np.bool_)

    for i in prange(rows):  # pylint: disable=not-an-iterable
        for j in prange(cols):  # pylint: disable=not-an-iterable
            null_gain_col = null_gain_distributions[:, j]  # Get column as array
            null_loss_col = null_loss_distributions[:, j]  # Get column as array
            if gain_props[i, j] > loss_props[i, j]:
                p_values[i, j] = np.sum(null_gain_col >= gain_props[i, j]) / n
                mech_gains[i, j] = 1
            else:
                p_values[i, j] = np.sum(null_loss_col >= loss_props[i, j]) / n
                mech_gains[i, j] = 0
    return p_values, mech_gains


class Mechanisms:
    """
    Handles the extraction and processing of molecular mechanisms from MutPred2 output files.

    This class provides methods to collect, format, and analyze mechanisms associated with
    mutations, including the computation of p-values and integration of probability scores.

    Attributes:
        None
    """

    @staticmethod
    def collect_mechanisms(catalog, catalog_job):
        """
        Collects molecular mechanisms from MutPred2 output files.

        Args:
            job (str): Path to the job directory.

        Returns:
            pd.DataFrame: A DataFrame containing mechanisms and their associated probabilities.
        """

        job = catalog_job.get_job_path()

        propx_pu = os.path.join(job, "output.txt.propX_pu_1.mat")
        property_scores = os.path.join(job, "output.txt.prop_scores_pu_1.mat")
        property_pvalues = os.path.join(job, "output.txt.prop_pvals_pu_1.mat")

        if os.path.exists(property_scores) and os.path.exists(property_pvalues):

            properties = Mechanisms.get_property_scores_and_pvalues(
                catalog, catalog_job, job
            )

            ## Test that properties are being calculated correctly
            test = False
            if test:
                propx = Mechanisms.vectorized_get_properties_from_propx(
                    catalog, catalog_job
                )
                for i, mechs in enumerate(properties):
                    for mech in mechs:
                        if mechs[mech] != propx[i][mech]:
                            print(mechs[mech])
                            print(propx[i][mech])

            return properties

        if os.path.exists(propx_pu):

            properties = Mechanisms.vectorized_get_properties_from_propx(
                catalog, catalog_job
            )

            return properties

        return None

    @staticmethod
    def get_property_scores_and_pvalues(catalog, catalog_job, job):
        """
        Retrieves property scores and p-values for mutations from MutPred2 output files.

        Args:
            job (str): Path to the job directory.

        Returns:
            pd.DataFrame: Dataframe containing mutations, scores, and associated mechanisms.

        1. Collect position from the output.txt.position_pu.mat files. This will be included
        on a per mechanism basis as "Impacted Amino Acid?"
        2. Create a prop_types for gain or loss.
        3. Track notes
        4. Track motifs
        """
        # catalog.property_names
        ## Collect predicted property scores
        property_scores = np.array(
            [
                score_array
                for p in range(catalog_job.get_num_files())
                for score_array in sio.loadmat(
                    os.path.join(job, f"output.txt.prop_scores_pu_{p+1}.mat")
                ).get("prop_scores_pu")
            ]
        )

        property_pvalues = np.array(
            [
                score_array
                for p in range(catalog_job.get_num_files())
                for score_array in sio.loadmat(
                    os.path.join(job, f"output.txt.prop_pvals_pu_{p+1}.mat")
                ).get("prop_pvals_pu")
            ]
        )

        property_types = np.array(
            [
                score_array - 1
                for p in range(catalog_job.get_num_files())
                for score_array in sio.loadmat(
                    os.path.join(job, f"output.txt.prop_types_pu_{p+1}.mat")
                ).get("prop_types_pu")
            ]
        )

        positions = catalog_job.get_positions()

        ## Format combination as a combined dictionary
        mechanisms = [
            pd.DataFrame(
                {
                    "Property": catalog.property_names,
                    "Posterior Probability": property_scores[i],
                    "P-value": property_pvalues[i],
                    "Effected Position": positions[i],
                    "Type": [
                        "Gain" if mech_type else "Loss"
                        for mech_type in property_types[i]
                    ],
                }
            )
            for i in range(property_scores.shape[0])
        ]

        return mechanisms

    @staticmethod
    def vectorized_get_properties_from_propx(catalog, catalog_job):
        """
        Extracts and processes property scores and p-values for given mutations.

        Args:
            job (str): Path to job directory.
            mutations (list): List of mutation identifiers.

        Returns:
            pd.DataFrame: Dataframe containing processed mutation data.
        """

        # Load all propX files
        job = catalog_job.get_job_path()

        propx_pu = np.array(
            [
                propx_array
                for p in range(catalog_job.get_num_files())
                for propx_array in sio.loadmat(
                    os.path.join(job, f"output.txt.propX_pu_{p+1}.mat")
                ).get("propX_pu")
            ]
        )

        # Identify loss/gain column indices
        loss_mask_idx = np.array(
            [
                i
                for i, col in enumerate(catalog.neutral_property_cols)
                if re.search(r"_loss$|Stability", col)
            ]
        )
        gain_mask_idx = np.array(
            [
                i
                for i, col in enumerate(catalog.neutral_property_cols)
                if re.search(r"_gain$|Stability", col)
            ]
        )

        # Convert to NumPy for fast operations
        loss_props = propx_pu[:, loss_mask_idx]
        gain_props = propx_pu[:, gain_mask_idx]

        # Compute max score across loss and gain
        max_scores = np.maximum(loss_props, gain_props)

        p_values, mech_gains = compute_p_values(
            loss_props,
            gain_props,
            catalog.null_gain_distributions_array,
            catalog.null_loss_distributions_array,
            catalog.n,
        )

        positions = catalog_job.get_positions()

        mechanisms = [
            pd.DataFrame(
                {
                    "Property": catalog.property_names,
                    "Posterior Probability": max_scores[i],
                    "P-value": p_values[i],
                    "Effected Position": positions[i],
                    "Type": [
                        "Gain" if mech_type else "Loss" for mech_type in mech_gains[i]
                    ],
                }
            )
            for i in range(max_scores.shape[0])
        ]

        return mechanisms

    @staticmethod
    def mechanisms_to_dict(cur_array: np.array, value_name, property_columns):
        """
        Convert a row of data into a nested dictionary structure.

        Args:
            row (pd.Series): A Pandas Series containing data, possibly with a "Substitution" column.
            value_name (str): The key under which values will be stored in the nested dictionary.

        Returns:
            dict: A dictionary where each key is a column name from `row`,
                and its value is another dictionary with `value_name` as the key
                and the corresponding row value as the value.

        Notes:
            - If the "Substitution" column exists, it is dropped before processing.
        """

        properties = {
            property_columns[i]: {value_name: r} for i, r in enumerate(cur_array)
        }

        return properties

    @staticmethod
    def mechanisms_np_to_dict(pr_array, pval_array, mech_gains, positions, cols):
        """
        Convert a NumPy array of property values into a dictionary.

        Args:
            row (np.ndarray): Array of property values.
            value_name (str): Name for the value field in the dictionary.

        Returns:
            dict: Dictionary with property names as keys and value_name as sub-key.
        """

        # Create a dictionary where each property maps to its corresponding value
        properties = {
            col: {
                "Posterior Probability": float(pr_array[i]),
                "P-value": float(pval_array[i]),
                "Effected Position": str(positions[i]),
                "Type": "Gain" if mech_gains[i] else "Loss",
            }
            for i, col in enumerate(cols)
        }

        return properties

    @staticmethod
    def combine_mechs(scores, pvals):
        """
        Combines two dictionaries where one contains mechanism scores and the other contains
        mechanism p-values.

        Args:
            row (pd.Series): A Pandas Series with 'Mechanisms_pr' and 'Mechanisms_p' dictionaries.

        Returns:
            dict: A dictionary combining probabilities and p-values for each mechanism.
        """

        merged = []

        def combine_dicts(score, pval):
            cur_dict = {}
            cur_dict.update(score)
            cur_dict.update(pval)
            return cur_dict

        for i, score_dict in enumerate(scores):
            pvals_dict = pvals[i]

            current_dict = {
                i: combine_dicts(score_dict[i], pvals_dict[i])
                for i in pvals_dict.keys()
            }
            merged.append(current_dict)

        return merged
