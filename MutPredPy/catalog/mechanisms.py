import os
import re
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

    @staticmethod
    def collect_mechanisms(catalog, job):
        """
        Collects molecular mechanisms from MutPred2 output files.

        Args:
            job (str): Path to the job directory.

        Returns:
            pd.DataFrame: A DataFrame containing mechanisms and their associated probabilities.
        """

        positions = catalog.get_positions()

        propx_pu = os.path.join(job, "output.txt.propX_pu_1.mat")
        property_scores = os.path.join(job, "output.txt.prop_scores_pu_1.mat")
        property_pvalues = os.path.join(job, "output.txt.prop_pvals_pu_1.mat")

        if os.path.exists(property_scores) and os.path.exists(property_pvalues):

            properties = Mechanisms.get_property_scores_and_pvalues(catalog, job)
            # propX = self.vectorized_get_properties_from_propX(catalog, job)

            # def evaluate(row):
            #    for k,v in row["Mechanisms_x"].items():
            #        if not row["Mechanisms_x"][k] == row["Mechanisms_y"][k]:
            #            print (k, row["Mechanisms_x"][k], row["Mechanisms_y"][k])
            #            exit()

            # testing = properties.merge(propX, on="Substitution", how="inner")
            # testing.apply(lambda row: evaluate(row), axis=1)
            # print (testing)
            # exit()

            return properties

        elif os.path.exists(propx_pu):

            properties = Mechanisms.vectorized_get_properties_from_propx(
                catalog, positions, job
            )

            return properties

    @staticmethod
    def get_property_scores_and_pvalues(catalog, job):
        """
        Retrieves property scores and p-values for mutations from MutPred2 output files.

        Args:
            job (str): Path to the job directory.

        Returns:
            pd.DataFrame: Dataframe containing mutations, scores, and associated mechanisms.

        1. Collect position from the output.txt.position_pu.mat files. This will be included on a per
        mechanism basis as "Impacted Amino Acid?"
        2. Create a prop_types for gain or loss.
        3. Track notes
        4. Track motifs
        """
        ## Collect predicted property scores
        prop_score_files = [
            prop_file
            for prop_file in os.listdir(job)
            if re.search(r"output.txt.prop_scores_pu_\d+.mat", prop_file)
        ]
        labeled_property_scores = pd.concat(
            [
                pd.DataFrame(
                    data=sio.loadmat(
                        os.path.join(job, f"output.txt.prop_scores_pu_{p+1}.mat")
                    ).get("prop_scores_pu"),
                    columns=catalog.property_names,
                )
                for p in range(len(prop_score_files))
            ]
        ).reset_index(drop=True)

        labeled_property_scores["Mechanisms"] = labeled_property_scores.apply(
            lambda row: Mechanisms.mechanisms_to_dict(row, "Posterior Probability"),
            axis=1,
        )
        labeled_property_scores = labeled_property_scores.filter(["Mechanisms"])

        ## Collect predicted property p-values
        prop_pvalues_files = [
            prop_file
            for prop_file in os.listdir(job)
            if re.search(r"output.txt.prop_pvals_pu_\d+.mat", prop_file)
        ]
        labeled_property_pvalues = pd.concat(
            [
                pd.DataFrame(
                    data=sio.loadmat(
                        os.path.join(job, f"output.txt.prop_pvals_pu_{p+1}.mat")
                    ).get("prop_pvals_pu"),
                    columns=catalog.property_names,
                )
                for p in range(len(prop_pvalues_files))
            ]
        ).reset_index(drop=True)

        labeled_property_pvalues["Mechanisms"] = labeled_property_pvalues.apply(
            lambda row: Mechanisms.mechanisms_to_dict(row, "P-value"), axis=1
        )
        labeled_property_pvalues = labeled_property_pvalues.filter(["Mechanisms"])

        ## Combine scores and p-values
        labeled_properies = labeled_property_scores.merge(
            labeled_property_pvalues,
            how="inner",
            suffixes=["_pr", "_p"],
        )

        ## Format combination as a combined dictionary
        labeled_properies["Mechanisms"] = labeled_properies.apply(combine_mechs, axis=1)
        labeled_properies = labeled_properies.filter(["Mechanisms"])

        return labeled_properies

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
        propx_pu_files = [
            os.path.join(job, f)
            for f in os.listdir(job)
            if re.search(r"output.txt.propX_pu_\d+.mat", f)
        ]

        # Load and concatenate all matrices
        prop_data = [sio.loadmat(f)["propX_pu"] for f in propx_pu_files]
        predicted_properties = pd.DataFrame(
            np.vstack(prop_data), columns=catalog.neutral_property_cols
        )

        # Identify loss/gain columns efficiently
        loss_mask = predicted_properties.columns.str.contains(r"_loss$|Stability")
        loss_columns = predicted_properties.columns[loss_mask]
        gain_columns = loss_columns.str.replace("_loss", "_gain")

        # Convert to NumPy for fast operations
        loss_props = predicted_properties[loss_columns].to_numpy()
        gain_props = predicted_properties[gain_columns].to_numpy()

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
        prop_columns = [col.replace("_loss", "") for col in loss_columns]
        mechanisms = [
            mechanisms_np_to_dict(
                max_scores[i], p_values[i], mech_gains[i], positions[i], prop_columns
            )
            for i in range(max_scores.shape[0])
        ]

        return mechanisms

    @staticmethod
    def mechanisms_to_dict(row, value_name):
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
        try:
            cur_row = row.drop("Substitution")
        except KeyError:
            cur_row = row

        properties = {i: {value_name: r} for i, r in cur_row.items()}

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
                "Effected Position": float(positions[i]),
                "Type": "Gain" if mech_gains[i] else "Loss",
            }
            for i, col in enumerate(cols)
        }

        return properties

    @staticmethod
    def combine_mechs(row):
        """
        Combines two dictionaries where one contains mechanism scores and the other contains
        mechanism p-values.

        Args:
            row (pd.Series): A Pandas Series with 'Mechanisms_pr' and 'Mechanisms_p' dictionaries.

        Returns:
            dict: A dictionary combining probabilities and p-values for each mechanism.
        """
        cur_dict = {}
        for k in row["Mechanisms_pr"].keys():
            cur_dict[k] = {}
            cur_dict[k].update(row["Mechanisms_pr"][k])
            cur_dict[k].update(row["Mechanisms_p"][k])

        return cur_dict
