"""
Input Processing Module for MutPredPy

This module provides functions for parsing, validating, and formatting input data
from various sources, including VEP and dbNSFP, to ensure compatibility with MutPred2.
"""

import re
import pandas as pd

from ..utils import collect_value, AminoAcidMap  # Importing from utils
from ..fasta import Protein


def is_vep_format(df):
    """
    Detects whether a pandas DataFrame follows the original VEP output format.

    Checks:
    - Mandatory VEP columns are present.
    - 'Extra' column exists for additional annotations.

    Parameters:
    df (pd.DataFrame): The DataFrame to check.

    Returns:
    bool: True if the DataFrame is in VEP format, False otherwise.
    """
    required_columns = {
        "#Uploaded_variation",
        "Location",
        "Allele",
        "Gene",
        "Feature",
        "Feature_type",
        "Consequence",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Extra",
    }
    return required_columns.issubset(set(df.columns))


def is_dbnsfp_format(df):
    """
    Detects whether a pandas DataFrame follows the dbNSFP output format.

    Checks:
    - Mandatory dbNSFP columns are present.

    Parameters:
    df (pd.DataFrame): The DataFrame to check.

    Returns:
    bool: True if the DataFrame is in dbNSFP format, False otherwise.
    """
    required_columns = {
        "#chr",
        "pos(1-based)",
        "ref",
        "alt",
        "aaref",
        "aaalt",
        "aapos",
        "rs_dbSNP",
        "genename",
        "Ensembl_geneid",
        "Ensembl_transcriptid",
        "Ensembl_proteinid",
        "HGVSp_ANNOVAR",
        "HGVSp_VEP",
        "VEP_canonical",
    }
    return required_columns.issubset(set(df.columns))


def parse_vep_output(df):
    """
    Parses a VEP output file, extracts the required fields,
    and returns a cleaned DataFrame with only the specified columns.

    Extracts:
    - #Uploaded_variation, Location, Allele, Gene, Consequence
    - Protein_position, Amino_acids, HGVSp, HGVSc, ENST, ENSP

    Parameters:
    vep_file (str): Path to the VEP output file (TSV format)

    Returns:
    pd.DataFrame: Parsed VEP data with only the specified columns.
    """

    # Check if DataFrame follows VEP format
    if not is_vep_format(df):
        raise ValueError("The input file does not match the expected VEP format.")

    # Add extracted values to DataFrame
    cols = ["ENST", "ENSP", "HGVSp", "HGVSc"]
    df[cols] = df["Extra"].apply(
        lambda x: pd.Series(
            index=["ENST", "ENSP", "HGVSp", "HGVSc"],
            data=collect_value(x, ["ENST", "ENSP", "HGVSp", "HGVSc"]),
        )
    )

    # Select and return only the required columns
    columns_to_keep = [
        "ENSP",
        "ENST",
        "Gene",
        "Consequence",
        "Protein_position",
        "Amino_acids",
        "HGVSp",
        "HGVSc",
    ]

    # Ensure all required columns exist in the DataFrame before selecting them
    df_filtered = df[[col for col in columns_to_keep if col in df.columns]]

    return df_filtered


def detect_proteins_transcripts_and_amino_acids(df):
    """
    Detects:
    1. HGVSp
    2. ENSP
    3. ENST
    4. ENSG
    5. Whether all three amino acid substitution columns (Reference AA, Alternative AA, Position)
        are present.
    6. Whether the detected Ensembl IDs are versioned.
    7. Whether multiple columns matched a given category.

    Parameters:
    df (pd.DataFrame): Input DataFrame.

    Returns:
    dict: A dictionary containing:
        - Booleans for Ensembl ID detections
        - Single column name where each type was found (or None if not found)
        - Booleans indicating whether each detected Ensembl ID is versioned
        - A `multiple` flag indicating if more than one column matched
    """
    patterns = Protein.protein_transcript_regex

    # Version detection regex
    versioned_pattern = re.compile(r"\.\d+")

    # Amino Acid Regex Patterns from AminoAcidMap
    aa_patterns = {
        "single_letter": re.compile(AminoAcidMap.AA_SINGLE_LETTER_REGEX),
        "three_letter": re.compile(AminoAcidMap.AA_THREE_LETTER_REGEX, re.IGNORECASE),
        "alternative": re.compile(AminoAcidMap.AA_ALTERNATIVE_REGEX, re.IGNORECASE),
        "combined": re.compile(AminoAcidMap.AA_COMBINED_REGEX, re.IGNORECASE),
        "full": re.compile(AminoAcidMap.AA_FULL_NOTATION, re.IGNORECASE),
        "numeric": re.compile(r"^(?:[0-3]?\d{1,4}|[1-3]\d{0,3})$"),
    }

    # Initialize results
    results = {
        k: {"found": False, "column": None, "versioned": False, "multiple": False}
        for k in patterns
    }
    results["has_amino_acid_columns"] = {
        "found": False,
        "columns": {"reference": None, "alternative": None, "position": None},
        "multiple": {"reference": False, "alternative": False, "position": False},
    }

    def update_results(category, col, versioned):
        """Helper function to update detection results for Ensembl IDs."""
        if results[category]["column"] is not None:
            results[category]["multiple"] = True
        results[category]["found"] = True
        results[category]["column"] = col
        results[category]["versioned"] = versioned

    def update_amino_acid_results(col):

        if results["has_amino_acid_columns"]["columns"]["reference"] is None:
            results["has_amino_acid_columns"]["columns"]["reference"] = col
        elif results["has_amino_acid_columns"]["columns"]["alternative"] is None:
            results["has_amino_acid_columns"]["columns"]["alternative"] = col
        else:
            results["has_amino_acid_columns"]["multiple"]["reference"] = True
            results["has_amino_acid_columns"]["multiple"]["alternative"] = True

        if (
            results["has_amino_acid_columns"]["columns"]["reference"]
            and results["has_amino_acid_columns"]["columns"]["reference"]
            and results["has_amino_acid_columns"]["columns"]["position"]
        ):
            results["has_amino_acid_columns"]["found"] = True

    def update_protein_position(col):
        if results["has_amino_acid_columns"]["columns"]["position"] is not None:
            results["has_amino_acid_columns"]["multiple"]["position"] = True
        results["has_amino_acid_columns"]["columns"]["position"] = col

        if (
            results["has_amino_acid_columns"]["columns"]["reference"]
            and results["has_amino_acid_columns"]["columns"]["reference"]
            and results["has_amino_acid_columns"]["columns"]["position"]
        ):
            results["has_amino_acid_columns"]["found"] = True

    threshold = len(df) * 0.9
    for col in df.columns:
        col_values = df[col].astype(str).dropna().str.strip()

        # Check for Protein and Transcript IDs
        for category, pattern in patterns.items():
            if col_values.str.match(pattern).sum() > threshold:
                update_results(
                    category, col, col_values.str.contains(versioned_pattern).any()
                )

        # Check for Full Mutation Notation columns
        if col_values.str.fullmatch(aa_patterns["full"]).sum() > threshold:
            update_amino_acid_results(col)
            update_amino_acid_results(col)
            update_protein_position(col)

        # Check for columns where reference and alternative are in the same column
        elif col_values.str.fullmatch(aa_patterns["combined"]).sum() > threshold:
            update_amino_acid_results(col)
            update_amino_acid_results(col)

        # Check for single, triple, and mixed Amino Acid columns
        elif (
            col_values.str.fullmatch(aa_patterns["single_letter"]).sum() > threshold
            or col_values.str.fullmatch(aa_patterns["three_letter"]).sum() > threshold
            or col_values.str.fullmatch(aa_patterns["alternative"]).sum() > threshold
        ):
            update_amino_acid_results(col)

        # Check for protein position columns
        elif col_values.str.fullmatch(aa_patterns["numeric"]).sum() > threshold:
            update_protein_position(col)

    # Check if all amino acid columns are present
    results["has_amino_acid_columns"]["found"] = all(
        results["has_amino_acid_columns"]["columns"].values()
    )

    return results


def collect_aa_substitution(row, col_mapping):
    """
    Constructs an amino acid substitution string from a given row.

    Args:
        row (pd.Series): A row containing reference, position, and alternative amino acids.
        col_mapping (dict): Dictionary mapping column names.

    Returns:
        str: Formatted amino acid substitution string (e.g., "A123T").
    """
    return f"{row[col_mapping['reference']]}\
        {row[col_mapping['position']]}\
            {row[col_mapping['alternative']]}"


def finalize_processed_input(df, validation_results):
    """
    Finalizes the processed input by extracting relevant protein and mutation information.

    Args:
        df (pd.DataFrame): The input DataFrame containing mutation data.
        validation_results (dict): Dictionary containing validation results
                                   for protein and transcript identification.

    Returns:
        pd.DataFrame: Processed DataFrame with extracted and formatted mutation substitutions.
    """
    processed_output = pd.DataFrame()

    for i, results in validation_results.items():
        # print(results)
        # exit()
        ## Identify proteins and transcripts

        if results["found"] and i != "has_amino_acid_columns":
            processed_output[i] = df[results["column"]]

        if results["found"] and (
            (i == "ENSP" and validation_results["HGVSp"]["column"] == results["column"])
            or (
                i == "ENST"
                and validation_results["HGVSc"]["column"] == results["column"]
            )
        ):
            processed_output[i] = (
                processed_output[results["column"]].str.split(":").str[0]
            )

        ## Handle amino acid substitutions
        elif i == "has_amino_acid_columns":

            ## If HGVSp IDs are found, ignore amino acid columns and just use HGVSp annotated mutations
            if validation_results["HGVSp"]["found"]:
                processed_output["Substitution"] = (
                    df[validation_results["HGVSp"]["column"]]
                    .str.split(":")
                    .str[-1]
                    .str.split("p.")
                    .str[-1]
                    .apply(
                        lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet")
                    )
                )

            ## Handle amino acid columns
            elif results["found"]:
                ref_col = results["columns"]["reference"]
                alt_col = results["columns"]["alternative"]
                pos_col = results["columns"]["position"]

                ## If ref, alt, and pos columns are all different, concat them together
                if alt_col not in (ref_col, pos_col):
                    processed_output["Substitution"] = df.apply(
                        lambda row, ref_col=ref_col, pos_col=pos_col, alt_col=alt_col: f"{row[ref_col]}{row[pos_col]}{row[alt_col]}"  # pylint:disable=C0301
                    ).apply(
                        lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet")
                    )

                ## If ref and alt are the same, try and parse them out.
                elif ref_col == alt_col and alt_col != pos_col:
                    processed_output["Substitution"] = df.apply(
                        lambda row, ref_col=ref_col, pos_col=pos_col, alt_col=alt_col: (
                            f'{row[ref_col].split("/")[0]}\
                                {row[pos_col]}\
                                    {row[alt_col].split("/")[0]}'
                            if "/" in row[ref_col]
                            else f"{row[ref_col]}{row[pos_col]}{row[alt_col]}"
                        )
                    ).apply(
                        lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet")
                    )

                ## If all three columns are the same, treat the column as a mutation
                ## column and convert to single letter notation.
                elif ref_col == alt_col and alt_col == pos_col:

                    processed_output["Substitution"] = df[ref_col].apply(
                        lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet")
                    )

                ## If none of the above, create a null mutation column
                else:
                    processed_output["Substitution"] = None
            else:
                processed_output["Substitution"] = None

    ## Deduplicate
    processed_output = processed_output.drop_duplicates()

    return processed_output


def check_validation_results(results, all_possible):
    """
    Validates the detection of essential columns in the input data.

    Args:
        results (dict): Dictionary containing validation results for detected protein,
                        transcript, and amino acid substitution columns.
        all_possible (bool): Flag indicating whether all possible amino acid substitutions
                             should be generated.

    Raises:
        Exception: If no suitable protein or transcript column is found.
        Exception: If no amino acid substitution information is found and `all_possible` is False.
    """
    protein_columns_found = [
        p for p, v in results.items() if v["found"] and p != "has_amino_acid_columns"
    ]
    amino_acid_columns_found = results["has_amino_acid_columns"]["found"]

    if not protein_columns_found:
        # If no suitable columns are found, raise an error
        raise ValueError(
            "No suitable column found with the required protein or transcript information \
            in the input file. Requires ENSP, ENST, RefSeq, or Uniprot"
        )

    if (
        not all_possible
        and not amino_acid_columns_found
        and not results["HGVSp"]["found"]
    ):
        raise ValueError(
            "No amino acid substitution information found in input file. Use the --all-possible \
                flag to generate all possible amino acid substitutions in each protein/transcript."
        )


def process_input(df, canonical, all_possible):
    """
    Processes input mutation data to extract relevant protein and amino acid
    substitution information.

    Args:
        df (pd.DataFrame): The input DataFrame containing mutation data.
        canonical (bool): Flag indicating whether to restrict processing to canonical isoforms.
        all_possible (bool): Flag indicating whether to generate all possible amino
            acid substitutions.

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: Processed mutation data.
            - dict: Validation results for detected identifiers and amino acid substitutions.

    Raises:
        Exception: If no valid protein or transcript information is detected.
        Exception: If no amino acid substitution information is found and `all_possible` is False.
    """
    if is_vep_format(df):
        processed_input_df = parse_vep_output(df)
    else:
        processed_input_df = df.copy()

    ## Collect validation results

    results = detect_proteins_transcripts_and_amino_acids(processed_input_df)

    check_validation_results(results, all_possible=all_possible)

    processed_input_df = finalize_processed_input(processed_input_df, results)

    if all_possible:
        ## Deduplicate all unique proteins/transcripts
        output = processed_input_df.drop_duplicates()

        ## Null all mutations
        output.loc[:, "Substitution"] = None

    else:
        ## Deduplicate
        output = processed_input_df.drop_duplicates()

        ## ===========================================
        ## TODO: Add canonical=True sequence filtering here
        canonical_status = canonical  # pylint: disable=W0612
        ## ===========================================

        ## currently filtering non-missense out here
        output["Mutation Type"] = output["Substitution"].apply(
            AminoAcidMap.mutation_type
        )
        output = output[output["Mutation Type"] == "Substitution"]

    return output, results
