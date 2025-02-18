import pandas as pd
import os
import re
from Bio import SeqIO
import importlib.resources as pkg_resources
import logging
logger = logging.getLogger()

from ..computing import lsf
from ..utils import utils as u

class Protein:

    # Protein ID Regex Patterns
    protein_transcript_regex = {
        "HGVSp":    re.compile(r"\bENSP\d{9,}(\.\d+)?:p.[A-Za-z]+\d+[A-Za-z]+\b"),
        "HGVSc":    re.compile(r"\bENST\d{9,}(\.\d+)?:c.\d+[NGCTAngcta]+>[NGCTAngcta]+\b"),
        "ENSP":     re.compile(r"\bENSP\d{9,}(\.\d+)?\b"),
        "ENST":     re.compile(r"\bENST\d{9,}(\.\d+)?\b"),
        "ENSG":     re.compile(r"\bENSG\d{9,}(\.\d+)?\b"),
        "Uniprot":  re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9])){1,2}(-\d{1,2})?\b"),
        "RefSeq_p": re.compile(r"\b[ANYXW]P\_\d+(\.\d+)?\b"),
        "RefSeq_r": re.compile(r"\b[NX]M\_\d+(\.\d+)?\b"),
        "RefSeq_g": re.compile(r"\b[AN][CGTWZ]\_\d+(\.\d+)?\b")

    }

    @staticmethod
    def identify_ids(info: str):
        """
        Identifies various protein-related IDs from the input string using predefined regex patterns.

        Args:
            info (str): The input string to search for IDs (e.g., protein accessions, gene IDs, etc.).

        Returns:
            dict: A dictionary where keys are the types of IDs (e.g., 'HGVSp', 'ENSP', 'Uniprot') 
                and values are the first matching ID found in the input string or None if not found.
        """

        # Initialize an empty dictionary to store the results
        results = dict()

        # Iterate over each ID type and its corresponding regex pattern from the class attribute
        for id_type,pattern in Protein.protein_transcript_regex.items():
            # Search for the pattern in the provided string
            id_match = pattern.search(info)

            # If a match is found, store the matched ID; otherwise, store None
            results[id_type] = id_match.group(0) if id_match else None

        # Return the dictionary containing the results for all ID types
        return results



def check_sequences(row, col_mapping):
    """Validate if the given sequence matches the expected reference amino acids at mutation positions.

    Args:
        row (pd.Series): Row containing 'mutation' and 'sequence'.

    Returns:
        tuple: (status (bool), sequence errors (list), mutation errors (list))
    """

    mutations = row[col_mapping["mutation_column"]].split(" ")
    sequence = row["sequence"]

    try:
        mutation_positions = [int(mut[1:-1]) for mut in mutations]  # Extract numeric positions
        ref_AA = [mut[0] for mut in mutations]  # Extract reference amino acids
    except (ValueError, IndexError):
        return False, ["Invalid mutation format"], []

    errors = {"Sequence Errors": [], "Mutation Errors": []}
    max_position = max(mutation_positions, default=0)  # Get the highest mutation position

    # Sequence length validation
    if len(sequence) < max_position:
        errors["Sequence Errors"].append(f"Sequence length ({len(sequence)}) is shorter than max mutation position ({max_position}).")

    if len(sequence) < 30:
        errors["Sequence Errors"].append("Sequence length is too short (<30).")

    # Invalid characters in sequence
    if "*" in sequence:
        errors["Sequence Errors"].append("Sequence contains stop codon (*).")

    if "U" in sequence:
        errors["Sequence Errors"].append("Sequence contains selenocysteine (U), which is not handled.")

    # Validate that reference amino acids match the sequence
    for i, pos in enumerate(mutation_positions):
        try:
            if ref_AA[i] != sequence[pos - 1]:
                errors["Mutation Errors"].append(f"Mismatch at {mutations[i]}: expected {ref_AA[i]}, found {sequence[pos - 1]}.")
        except:
            print (i, pos)
            #print (ref_AA)
            print (len(sequence))
            #print (sequence[588])
            print (row)
            exit()

    return (not errors["Mutation Errors"] and not errors["Sequence Errors"]), errors["Sequence Errors"], errors["Mutation Errors"]


def clean_FASTA_sequence(sequence):
    """Cleans FASTA sequence by removing invalid characters.

    Replaces:
    - Stop codon (*) with alanine (A)
    - Selenocysteine (U) with alanine (A)
    - Removes leading 'X' (unknown residue)

    Args:
        sequence (str): Protein sequence.

    Returns:
        str: Cleaned sequence.
    """
    try:
        return sequence.replace("*", "A").replace("U", "A").lstrip("X")
    except AttributeError:
        print (sequence)
        print (type(sequence))
        exit()


def data_quality_check(prepare, data, col_mapping):

    data[["status","Sequence Errors", "Mutation Errors"]] = data.apply(lambda x: pd.Series(check_sequences(x, col_mapping)), axis=1)

    ## Check if any sequences are invalid
    not_passed = data[data["status"]==False]
    if not not_passed.empty:
        
        ## Log errors
        log_error_message = f"{not_passed['Substitution'].apply(lambda x: len(x.split(' '))).sum()} mutations in {len(not_passed)} proteins failed the quality check and have been removed." 
        logger.warning(f"{log_error_message} (additional info in {prepare.get_log_folder()}/errors.log)")
        
        if prepare.dry_run:
            logger.dry_run(f"Would've written logs to {prepare.get_log_folder()}/errors.log)")
        else:
            u.check_pass_error_logging(not_passed, prepare.get_log_folder(), col_mapping, log_error_message)
        
        ## Drop variants that do not pass the quality checks
        data = data[data["status"]]

        if data.empty:
            raise Exception("No sequences or mutations passed data quality check.")

    return data

def sequence_quality_check(data):
    """
    Cleans existing sequences in preparation for MutPred2
    Checks that the quality of the sequences in ready for MutPred2
    """
    data["sequence"] = data["sequence"].apply(lambda x: clean_FASTA_sequence(x))

    return data


def alignment_score(data, col_mapping):
    """Calculate alignment score between reference and actual sequence.

    Args:
        data (pd.Series): Row containing 'mutation' and 'sequence'.

    Returns:
        float: Fraction of correctly aligned residues.
    """
    
    mutations = data[col_mapping["mutation_column"]].split(" ")
    positions = [int(re.search(r'\d+', mut).group()) - 1 for mut in mutations]

    return sum(1 for pos, ref in zip(positions, mutations) if data["sequence"][pos] == ref[0]) / len(positions)


def get_fasta_location(assembly="GRCh38"):

    if assembly.lower() in ["hg19","grch37"]:
        cur_assembly = "GRCh37"
    else:
        cur_assembly = "GRCh38"

    pep_file = pkg_resources.files("mutpredpy").joinpath(f"pkg_resources/Homo_sapiens.{cur_assembly}.combined.pep.all.fa")
    
    if not os.path.exists(pep_file):
        raise FileNotFoundError(f"Reference FASTA file not found: {pep_file}")
    
    return pep_file


def read_mutpred_input_fasta(faa_file):
    """Reads input FASTA file and extracts sequence information.

    Args:
        faa_file (str): Path to the FASTA file.

    Returns:
        pd.DataFrame: Parsed FASTA sequences with metadata.
    """
    if not os.path.exists(faa_file):
        raise FileNotFoundError(f"FASTA file not found: {faa_file}")

    output = []
    temp = {"ID": "", "Substitution": "", "sequence": ""}

    with open(faa_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                if temp["ID"]:
                    output.append(temp.copy())
                header = line[1:].strip().split(" ")
                temp.update({"ID": header[0], "Substitution": ",".join(header[1:]), "sequence": ""})
            else:
                temp["sequence"] += line.strip()
        if temp["ID"]:
            output.append(temp)

    df = pd.DataFrame(output).drop_duplicates()
    df["num_substitution"] = df["Substitution"].apply(lambda x: len(x.split(",")))
    return df


def read_fasta(fasta_path):

    output = []
    records = SeqIO.parse(fasta_path, "fasta")
    for rec in records:

        temp = Protein.identify_ids(rec.description)
        temp["sequence"] = str(rec.seq)
        output.append(temp)
    
    output_df = pd.DataFrame(output)

    memory_function = lsf.memory_estimate_function()
    output_df["Memory Estimate (MB)"] = output_df["sequence"].apply(lambda x: memory_function(len(x)))

    time_function = lsf.time_estimate_function()
    output_df["Time per Mutation (hrs)"] = output_df["sequence"].apply(lambda x: time_function(len(x)))

    threshold = 0.8
    output_df = output_df.drop(columns=[col for col in output_df.columns if output_df[col].isnull().mean() > threshold])

    return output_df



def collect_ensembl_fasta(assembly):
    """Parses the ensembl protein FASTA file and extracts sequence metadata.

    Args:
        primary (str, optional): Primary key for identifying proteins.
        drop_dups (bool, optional): Whether to remove duplicate sequences. Defaults to True.

    Returns:
        pd.DataFrame: DataFrame containing protein ID, gene information, and sequence details.
    """

    # Collect the ensembl fasta file for the designated assembly
    if assembly.lower() in ["hg19","grch37"]:
        cur_assembly = "GRCh37"
    else:
        cur_assembly = "GRCh38"

    pep_file = pkg_resources.files("mutpredpy").joinpath(f"pkg_resources/Homo_sapiens.{cur_assembly}.combined.pep.all.fa")
    
    if not os.path.exists(pep_file):
        raise FileNotFoundError(f"Reference FASTA file not found: {pep_file}")
    
    output_df = read_fasta(pep_file)

    return output_df


def filter_non_missense():
    pass