import pandas as pd
import re
from ..utils import collect_value, AminoAcidMap  # Importing from utils

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
        "#Uploaded_variation", "Location", "Allele", "Gene", "Feature",
        "Feature_type", "Consequence", "Protein_position", "Amino_acids", "Codons", "Extra"
    }
    return required_columns.issubset(set(df.columns))


def is_dbNSFP_format(df):
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
        "#chr","pos(1-based)","ref","alt","aaref","aaalt","aapos","rs_dbSNP","genename","Ensembl_geneid",
        "Ensembl_transcriptid","Ensembl_proteinid","HGVSp_ANNOVAR","HGVSp_VEP","VEP_canonical"
    }
    return required_columns.issubset(set(df.columns))


def is_VCF_format(df):
    """
    Detects whether a pandas DataFrame follows the VCF format.

    Checks:
    - Mandatory VCF columns are present.

    Parameters:
    df (pd.DataFrame): The DataFrame to check.

    Returns:
    bool: True if the DataFrame is in VCF format, False otherwise.
    """
    required_columns = {
        "#CHROM","POS","REF","ALT"
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
    cols = ["ENST","ENSP","HGVSp","HGVSc"]
    df[cols] = df["Extra"].apply(lambda x: pd.Series(index=["ENST","ENSP","HGVSp","HGVSc"], data=collect_value(x, ["ENST","ENSP","HGVSp","HGVSc"])))

    # Select and return only the required columns
    columns_to_keep = [
        "ENSP", "ENST", "Gene", "Consequence",
        "Protein_position", "Amino_acids", "HGVSp", "HGVSc"
    ]

    # Ensure all required columns exist in the DataFrame before selecting them
    df_filtered = df[[col for col in columns_to_keep if col in df.columns]]

    return df_filtered


def parse_dbNSFP_output(df):

    # Select and return only the required columns
    columns_to_keep = [
        "Ensembl_geneid","Ensembl_transcriptid","Ensembl_proteinid","aaref","aaalt","aapos","rs_dbSNP","genename",
        "HGVSp_ANNOVAR","HGVSp_VEP","VEP_canonical"
    ]

    # Ensure all required columns exist in the DataFrame before selecting them
    df_filtered = df[[col for col in columns_to_keep if col in df.columns]]

    return df_filtered


def detect_ensembl_and_amino_acids(df):
    """
    Detects:
    1. HGVSp
    2. ENSP
    3. ENST
    4. ENSG
    5. Whether all three amino acid substitution columns (Reference AA, Alternative AA, Position) are present.
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
    # Ensembl ID Regex Patterns
    patterns = {
        "HGVSp": re.compile(r"ENSP\d{9,}(\.\d+)?:p.[A-Za-z]+\d+[A-Za-z]+"),
        "HGVSc": re.compile(r"ENST\d{9,}(\.\d+)?:c.\d+[NGCTAngcta]+>[NGCTAngcta]+"),
        "ENSP": re.compile(r"ENSP\d{9,}(\.\d+)?$"),
        "ENST": re.compile(r"ENST\d{9,}(\.\d+)?$"),
        "Uniprot": re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
        "ENSG": re.compile(r"ENSG\d{9,}(\.\d+)?$")
    }
    
    # Version detection regex
    versioned_pattern = re.compile(r"\.\d+")

    # Amino Acid Regex Patterns from AminoAcidMap
    aa_patterns = {
        "single_letter": re.compile(AminoAcidMap.AA_SINGLE_LETTER_REGEX),
        "three_letter": re.compile(AminoAcidMap.AA_THREE_LETTER_REGEX, re.IGNORECASE),
        "alternative": re.compile(AminoAcidMap.AA_ALTERNATIVE_REGEX, re.IGNORECASE),
        "combined": re.compile(AminoAcidMap.AA_COMBINED_REGEX, re.IGNORECASE),
        "full": re.compile(AminoAcidMap.AA_FULL_NOTATION, re.IGNORECASE),
        "numeric": re.compile(r"^\d+$")
    }

    # Initialize results
    results = {
        "HGVSp": {"found": False, "column": None, "versioned": False, "multiple": False},
        "HGVSc": {"found": False, "column": None, "versioned": False, "multiple": False},
        "ENSP": {"found": False, "column": None, "versioned": False, "multiple": False},
        "ENST": {"found": False, "column": None, "versioned": False, "multiple": False},
        "Uniprot": {"found": False, "column": None, "versioned": False, "multiple": False},
        "ENSG": {"found": False, "column": None, "versioned": False, "multiple": False},
        "has_amino_acid_columns": {
            "found": False,
            "columns": {"reference": None, "alternative": None, "position": None},
            "multiple": {"reference": False, "alternative": False, "position": False}
        }
    }

    def update_results(category, col, versioned):
        """Helper function to update detection results for Ensembl IDs."""
        if results[category]["column"] is not None:
            results[category]["multiple"] = True
        results[category]["found"] = True
        results[category]["column"] = col
        results[category]["versioned"] = versioned

    def update_amino_acid_results(col):
        results["has_amino_acid_columns"]["found"] = True
        if results["has_amino_acid_columns"]["columns"]["reference"] is None:
            results["has_amino_acid_columns"]["columns"]["reference"] = col
        elif results["has_amino_acid_columns"]["columns"]["alternative"] is None:
            results["has_amino_acid_columns"]["columns"]["alternative"] = col
        else:
            results["has_amino_acid_columns"]["multiple"]["reference"] = True
            results["has_amino_acid_columns"]["multiple"]["alternative"] = True
    
    def update_protein_position(col):
        if results["has_amino_acid_columns"]["columns"]["position"] is not None:
            results["has_amino_acid_columns"]["multiple"]["position"] = True
        results["has_amino_acid_columns"]["columns"]["position"] = col

    
    for col in df.columns:
        col_values = df[col].astype(str).dropna().str.strip()
        
        # Check for Ensembl IDs
        for category, pattern in patterns.items():
            if col_values.str.match(pattern).any():
                update_results(category, col, col_values.str.contains(versioned_pattern).any())
        
        # Check for Full Mutation Notation columns
        if col_values.str.fullmatch(aa_patterns["full"]).any():
            update_amino_acid_results(col)
            update_amino_acid_results(col)
            update_protein_position(col)
        
        # Check for columns where reference and alternative are in the same column
        elif col_values.str.fullmatch(aa_patterns["combined"]).any():
            update_amino_acid_results(col)
            update_amino_acid_results(col)

        # Check for single, triple, and mixed Amino Acid columns
        elif col_values.str.fullmatch(aa_patterns["single_letter"]).any() or \
            col_values.str.fullmatch(aa_patterns["three_letter"]).any() or \
            col_values.str.fullmatch(aa_patterns["alternative"]).any():
            update_amino_acid_results(col)

        # Check for protein position columns
        elif col_values.str.fullmatch(aa_patterns["numeric"]).any():
            update_protein_position(col)

    # Check if all amino acid columns are present
    results["has_amino_acid_columns"]["found"] = all(results["has_amino_acid_columns"]["columns"].values())
    
    return results


def valid_protein_information(df):
    """
    Selects columns to keep from the `detect_ensembl_and_amino_acids` output 
    based on the following priority:
    1. Versioned ENSP with `has_amino_acid_columns` found
    2. Versioned ENSP with mutations
    3. Versioned ENST with `has_amino_acid_columns` found
    4. Unversioned ENSP with amino acid substitution information
    5. Unversioned ENST with amino acid substitution information
    6. Throws an error if none exist.

    Parameters:
    detect_results (dict): Output from `detect_ensembl_and_amino_acids`.

    Returns:
    dict: A dictionary with column names for ENSP/ENST and amino acid columns.
          Example: {"id_column": "ENSP_column", "reference": "ref_column", "alternative": "alt_column", "position": "pos_column"}
    """

    ## results of detect_ensembl_and_amino_acids
    validation_results = detect_ensembl_and_amino_acids(df)

    # Define the priority order with tuples (id_type, versioned, must_have_aa_columns)
    priorities = [
        ("HGVSp", True, True),  # Versioned ENSP with `has_amino_acid_columns`
        ("HGVSp", True, False),  # Versioned ENSP without checking for amino acid columns
        ("HGVSc", True, True),  # Versioned ENST with `has_amino_acid_columns`
        ("HGVSc", True, False),  # Versioned ENST without checking for amino acid columns
        ("ENST", True, True),  # Versioned ENST with `has_amino_acid_columns`
        ("ENSP", False, True),  # Unversioned ENSP with amino acid substitution
        ("ENST", False, True),  # Unversioned ENST with amino acid substitution
        ("ENSP", True, False),
        ("ENST", True, False),
        ("ENSP", False, False),
        ("ENST", False, False)
    ]

    # Iterate over each priority and check if it meets the criteria
    for id_type, versioned, must_have_aa_columns in priorities:
        if validation_results[id_type]["found"] and validation_results[id_type]["versioned"] == versioned:
            if must_have_aa_columns and not validation_results["has_amino_acid_columns"]["found"]:
                continue  # Skip this priority if amino acid columns are required but not found
            
            # Return the relevant column names
            return {
                "id_type": id_type,
                "id_column": validation_results[id_type]["column"],
                "reference": validation_results["has_amino_acid_columns"]["columns"]["reference"],
                "alternative": validation_results["has_amino_acid_columns"]["columns"]["alternative"],
                "position": validation_results["has_amino_acid_columns"]["columns"]["position"]
            }

    return False


def finalize_processed_input(df, validation_results):

    processed_output = pd.DataFrame()
    column_mapping = dict()

    if validation_results["id_type"] == "HGVSp":
        processed_output["ENSP"] = df[validation_results["id_column"]].str.split(":").str[0]
        processed_output["Substitution"] = df[validation_results["id_column"]].str.split(":").str[-1].str.split("p.").str[-1].apply(AminoAcidMap.mutation_mapping)
        
        column_mapping["id_column"] = "ENSP"
        column_mapping["mutation_column"] = "Substitution"

    processed_output = processed_output.drop_duplicates()

    return processed_output, column_mapping


def process_input(df, canonical, all_possible):

    if is_vep_format(df):
        processed_input_df = parse_vep_output(df)

    elif is_dbNSFP_format(df):
        processed_input_df = parse_dbNSFP_output(df)

    else:
        processed_input_df = df.copy()

    ## Collect validation results 

    results = valid_protein_information(processed_input_df)

    if is_VCF_format(df) and not valid_protein_information(df):
        # If no suitable columns are found, raise an error
        raise ValueError("No suitable column set found with the required amino acid substitution information.")
    
    ## ===========================================
    ## TODO: Handle the all_possible=True scenario
    ## ===========================================
    
    ## Deduplicate
    processed_input_df = processed_input_df.drop_duplicates()

    output, col_mapping = finalize_processed_input(processed_input_df, results)

    ## ===========================================
    ## TODO: Add canonical=True sequence filtering here
    ## ===========================================

    ## currently filtering non-missense out here
    output["Mutation Type"] = output["Substitution"].apply(AminoAcidMap.mutation_type)
    substitutions = output[output["Mutation Type"]=="Substitution"]
    
    return substitutions, col_mapping

