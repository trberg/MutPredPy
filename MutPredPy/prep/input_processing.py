import pandas as pd
import re
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



def detect_proteins_transcripts_and_amino_acids(df):
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
        "numeric": re.compile(r"^(?:[0-3]?\d{1,4}|[1-3]\d{0,3})$")
    }

    # Initialize results
    results = {k: {"found": False, "column": None, "versioned": False, "multiple": False} for k in patterns.keys()}
    results["has_amino_acid_columns"] = {
        "found": False, 
        "columns": {"reference": None, "alternative": None, "position": None}, 
        "multiple": {"reference": False, "alternative": False, "position": False}
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

        if results["has_amino_acid_columns"]["columns"]["reference"] and \
                results["has_amino_acid_columns"]["columns"]["reference"] and \
                results["has_amino_acid_columns"]["columns"]["position"]:
            results["has_amino_acid_columns"]["found"] = True
    
    def update_protein_position(col):
        if results["has_amino_acid_columns"]["columns"]["position"] is not None:
            results["has_amino_acid_columns"]["multiple"]["position"] = True
        results["has_amino_acid_columns"]["columns"]["position"] = col

        if results["has_amino_acid_columns"]["columns"]["reference"] and  \
                results["has_amino_acid_columns"]["columns"]["reference"] and \
                results["has_amino_acid_columns"]["columns"]["position"]:
            results["has_amino_acid_columns"]["found"] = True


    threshold = len(df)*0.9
    for col in df.columns:
        col_values = df[col].astype(str).dropna().str.strip()
        
        # Check for Protein and Transcript IDs
        for category, pattern in patterns.items():
            if col_values.str.match(pattern).sum() > threshold:
                update_results(category, col, col_values.str.contains(versioned_pattern).any())
        
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
        elif col_values.str.fullmatch(aa_patterns["single_letter"]).sum() > threshold or \
            col_values.str.fullmatch(aa_patterns["three_letter"]).sum() > threshold or \
            col_values.str.fullmatch(aa_patterns["alternative"]).sum() > threshold:
            update_amino_acid_results(col)

        # Check for protein position columns
        elif col_values.str.fullmatch(aa_patterns["numeric"]).sum() > threshold:
            update_protein_position(col)

    # Check if all amino acid columns are present
    results["has_amino_acid_columns"]["found"] = all(results["has_amino_acid_columns"]["columns"].values())
    
    return results


def collect_AA_substitution(row, col_mapping):

    return f"{row[col_mapping['reference']]}{row[col_mapping['position']]}{row[col_mapping['alternative']]}"


def finalize_processed_input(df, validation_results):

    processed_output = pd.DataFrame()
    
    for i,results in validation_results.items():

        ## Identify proteins and transcripts
        if results["found"] and i != "has_amino_acid_columns":
            processed_output[i] = df[results['column']]

        ## Handle amino acid substitutions
        elif i == "has_amino_acid_columns":

            ## If HGVSp IDs are found, ignore amino acid columns and just HGVSp annotated mutations
            if validation_results["HGVSp"]["found"]:
                processed_output["Substitution"] = df[validation_results["id_column"]].str.split(":").str[-1].str.split("p.").str[-1].apply(lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet"))

            ## Handle amino acid columns
            elif results["found"]:
                ref_col = results["columns"]["reference"]
                alt_col = results["columns"]["alternative"]
                pos_col = results["columns"]["position"]

                ## If ref, alt, and pos columns are all different, concat them together
                if ref_col != alt_col and alt_col != pos_col:
                    processed_output["Substitution"] = df.apply(lambda row: f'{row[ref_col]}{row[pos_col]}{row[alt_col]}').apply(lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet"))
                
                ## If ref and alt are the same, try and parse them out.
                elif ref_col == alt_col and alt_col != pos_col:
                    processed_output["Substitution"] = df.apply(
                        lambda row: f'{row[ref_col].split("/")[0]}{row[pos_col]}{row[alt_col].split("/")[0]}' if "/" in row[ref_col] \
                            else f'{row[ref_col]}{row[pos_col]}{row[alt_col]}'
                        ).apply(lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet"))
                    
                ## If all three columns are the same, treat the column as a mutation column and convert to single letter notation.
                elif ref_col == alt_col and alt_col == pos_col:
                    
                    processed_output["Substitution"] = df[ref_col].apply(lambda x: AminoAcidMap.mutation_mapping(x, direction="singlet"))
                
                ## If none of the above, create a null mutation column
                else:
                    processed_output["Substitution"] = None
            else:
                processed_output["Substitution"] = None
    
    ## Deduplicate
    processed_output = processed_output.drop_duplicates()

    return processed_output


def check_validation_results(results, all_possible):

    protein_columns_found = [p for p,v in results.items() if v["found"] and p != "has_amino_acid_columns"]
    amino_acid_columns_found = results["has_amino_acid_columns"]["found"]

    if not protein_columns_found:
        # If no suitable columns are found, raise an error
        raise Exception("No suitable column found with the required protein or transcript information in the input file. Requires ENSP, ENST, RefSeq, or Uniprot")
    
    if not all_possible and not amino_acid_columns_found:
        raise Exception("No amino acid substitution information found in input file. Use the --all-possible flag to generate all possible amino acid substitutions in each protein/transcript.")

def process_input(df, canonical, all_possible):

    if is_vep_format(df):
        processed_input_df = parse_vep_output(df)

    else:
        processed_input_df = df.copy()

    ## Collect validation results 

    results = detect_proteins_transcripts_and_amino_acids(processed_input_df)

    check_validation_results(results, all_possible=all_possible)
    
    
    ## ===========================================
    ## TODO: Handle the all_possible=True scenario
    ## ===========================================
    
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
        ## ===========================================

        ## currently filtering non-missense out here
        output["Mutation Type"] = output["Substitution"].apply(AminoAcidMap.mutation_type)
        output = output[output["Mutation Type"]=="Substitution"]
    
    return output, results

