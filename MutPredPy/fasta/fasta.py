import pandas as pd
import numpy as np
import os
import re
import logging
import importlib.resources as pkg_resources



# Configure logging for debugging and tracking
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


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
        if ref_AA[i] != sequence[pos - 1]:
            errors["Mutation Errors"].append(f"Mismatch at {mutations[i]}: expected {ref_AA[i]}, found {sequence[pos - 1]}.")

    return not errors["Mutation Errors"], errors["Sequence Errors"], errors["Mutation Errors"]


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
    return sequence.replace("*", "A").replace("U", "A").lstrip("X")


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


def memory_estimate_function():
    """Load pre-trained memory usage model.

    Returns:
        np.poly1d: Polynomial function to estimate memory usage.
    """
    path = pkg_resources.files("mutpredpy").joinpath("pkg_resources/memory_usage.npy")
    #path = os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/memory_usage.npy"))
    return np.poly1d(np.load(path))


def time_estimate_function():
    """Load pre-trained time estimate function.

    Returns:
        np.poly1d: Polynomial function to estimate runtime.
    """
    path = pkg_resources.files("mutpredpy").joinpath("pkg_resources/sequence_time.npy")
    #path = os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/sequence_time.npy"))
    return np.poly1d(np.load(path))


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


def collect_ensembl_fasta(primary, assembly, versioned=True, drop_dups=True):
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
    

    output = []
    temp = dict()

    with open(pep_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):  # Start of a new protein entry
                if temp.get("ENSP"):
                    output.append(temp.copy())  # Store previous entry

                header_info = line[1:].strip().split(" ")
                keys = {k.split(":")[0]: k.split(":")[1] for k in header_info if ":" in k}

                # Extract relevant information from the header
                temp["ENSP"] = header_info[0]
                temp["ENSG"] = keys.get("gene", "")
                temp["ENST"] = keys.get("transcript", "")
                temp["sequence"] = ""  # Reset sequence for new entry
            else:
                temp["sequence"] += line.strip()  # Append sequence

        # Add the last parsed sequence
        if temp.get("ENSP"):
            output.append(temp)

    # Convert to DataFrame
    output_df = pd.DataFrame(output)

    # Drop duplicate sequences based on primary key and sequence content
    if drop_dups:
        output_df = output_df.drop_duplicates(subset=[primary, "sequence"], keep="first")

    #if not versioned:
    #    # Extract unversioned base protein ID
    #    output_df[primary] = output_df[primary].str.split(".").str[0]

    # Compute sequence lengths
    output_df["length"] = output_df["sequence"].apply(len)

    # Estimate memory and computation time
    memory_function = memory_estimate_function()
    output_df["Memory Estimate (MB)"] = output_df["sequence"].apply(lambda x: memory_function(len(x)))

    time_function = time_estimate_function()
    output_df["Time per Mutation (hrs)"] = output_df["sequence"].apply(lambda x: time_function(len(x)))

    filtered_output_df = output_df.filter(items=[primary, "sequence", "Memory Estimate (MB)", "Time per Mutation (hrs)"])

    return filtered_output_df


def filter_non_missense():
    pass