import hashlib
import yaml
import re
import os
import logging
logger = logging.getLogger()

import importlib.resources as pkg_resources


# Class for amino acid mappings (3-letter to 1-letter and vice versa)
class AminoAcidMap:
    """A class to store mappings between three-letter and one-letter amino acid codes."""
    THREE_TO_ONE = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Sec": "U", "Thr": "T", "Trp": "W", "Tyr": "Y",
        "Val": "V", "Ter": "X", "%3D": "="  # %3D represents "=" in URL encoding
    }
    ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}  # Reverse mapping (1-letter to 3-letter)

    # Regex Patterns
    SINGLE_AA = r"ACDEFGHIKLMNPQRSTVWY" #UX=
    TRIPLETS_AA = r"Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Sec|Thr|Trp|Tyr|Val"

    AA_SINGLE_LETTER_REGEX = f"^[{SINGLE_AA}]$"
    AA_THREE_LETTER_REGEX = f"^(?i)({TRIPLETS_AA})$"
    AA_COMBINED_REGEX = f"^[{SINGLE_AA}]/[{SINGLE_AA}]$"
    AA_ALTERNATIVE_REGEX = (
        f"^([{SINGLE_AA}]|{TRIPLETS_AA})"
        f"(\s*,\s*[{SINGLE_AA}]|{TRIPLETS_AA})*$"
    )
    AA_FULL_NOTATION = f"^(?:[{SINGLE_AA}]|(?:{TRIPLETS_AA}))\d+(?:[{SINGLE_AA}]|(?:{TRIPLETS_AA}))$"

    SINGLE_LETTER_AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    
    @staticmethod
    def map_amino_acid(code):
        """Map a single amino acid code (1-letter to 3-letter or vice versa)."""
        if len(code) == 1:  # Single to triple
            return AminoAcidMap.ONE_TO_THREE.get(code, code)
        elif len(code) == 3:  # Triple to single
            return AminoAcidMap.THREE_TO_ONE.get(code.capitalize(), code)
        return code
    

    @staticmethod
    def mutation_mapping(x, direction="auto"):
        """
        If direction=auto, Automatically detect and convert a mutation string between single-letter and three-letter amino acid notation.
        if direction=singlet, converts a mutation string to single-letter amino acid notation
        if direction=triplet, converts a mutation string to three-letter amino acid notation

        Args:
            x (str): Mutation string (e.g., "A123T" or "Ala123Thr").

        Returns:
            str: Mutation in the opposite notation format.
        """
        try:
            mut = x
            loc = re.search(r'\d+', mut).group()  # Extract numeric position
            ref, alt = mut.split(loc)  # Split reference and alternate residues

            # Detect if the input is single-letter or three-letter notation
            if len(ref) == 1 and len(alt) == 1:  # Single-letter input
                if direction == "singlet":
                    ref_mapped = ref
                    alt_mapped = alt
                else:
                    ref_mapped = AminoAcidMap.map_amino_acid(ref)
                    alt_mapped = AminoAcidMap.map_amino_acid(alt)
            elif len(ref) == 3 and len(alt) == 3:
                if direction == "triplet":
                    ref_mapped = ref
                    alt_mapped = alt
                else:
                    ref_mapped = AminoAcidMap.map_amino_acid(ref)
                    alt_mapped = AminoAcidMap.map_amino_acid(alt)
            else:
                raise ValueError("Inconsistent mutation format detected.")

            return f"{ref_mapped}{loc}{alt_mapped}"
        
        except (AttributeError, ValueError) as e:
            return "."
        

    @staticmethod
    def mutation_type(mutation):
        """
        Annotates a mutation with its type according to HGVS nomenclature.

        Args:
            mutation (str): Mutation in standard HGVS format, e.g., "A25T", "Ala25Thr", "c.123A>G", "NP_0123456.1:p.Arg97fs", etc.

        Returns:
            str: Mutation type annotation (e.g., "Substitution", "Deletion", "Insertion", "Delins", "Duplication", "Frameshift", "Repeated Sequence").
        """
        
        """
        TODO: THIS WILL NEED TO BE ADJUSTED FOR FUTURE MUTATION TYPES
        """
        
        # Define patterns for different mutation types
        patterns = {
            "Substitution": re.compile(f"^[{AminoAcidMap.SINGLE_AA}]\d+[{AminoAcidMap.SINGLE_AA}]$|^(?:{AminoAcidMap.TRIPLETS_AA})\d+(?:{AminoAcidMap.TRIPLETS_AA})$"),
            #"Deletion": re.compile(r"^c\.\d+_?\d*del[A-Z]*$"),
            #"Insertion": re.compile(r"^c\.\d+_?\d*ins[A-Z]+$"),
            #"Delins": re.compile(r"^c\.\d+_?\d*delins[A-Z]+$"),
            #"Duplication": re.compile(r"^c\.\d+_?\d*dup[A-Z]+$"),
            #"Frameshift": re.compile(r"^(c\.\d+[A-Z]>[A-Z]fs$|NP_\d+\.\d+:p\.[A-Za-z]{3}\d+fs(?:Ter\d+)?$)"),  # Combined cDNA and protein frameshift
            #"Repeated Sequence": re.compile(r"^c\.\d+_?\d*\([A-Z]+\)\d+$")  # Repeated sequence like c.123_124(A)4
        }
        
        # Check the mutation against each pattern
        for mutation_type, pattern in patterns.items():
            if pattern.match(mutation):
                return mutation_type
        
        # Return "Unknown" if no pattern matches
        return "Unknown"
    

    @staticmethod
    def all_possible_mutations(sequence):
        
        mutations = []
        for i in range(len(sequence)):
            pos = i + 1
            ref = sequence[i]
            for alt in AminoAcidMap.SINGLE_AA:
                if ref != alt:
                    mutation = f"{ref}{pos}{alt}"
                    mutations.append(mutation)

        return " ".join(mutations)
    

    @staticmethod
    def get_all_possible_mutations(df, col_mapping):

        df[col_mapping["mutation_column"]] = df["sequence"].apply(AminoAcidMap.all_possible_mutations)

        return df


def get_seq_hash(sequence):
    return hashlib.md5(sequence.encode(encoding='utf-8')).hexdigest()


def collect_value(x, keys):
    """
    Extracts a value corresponding to a given key from a semicolon-separated string of key-value pairs.
    Supports both single string input and a list of such strings.

    Parameters:
    x (str): A string containing key-value pairs in the format "key=value" separated by semicolons.
    val (str or list): The key or list of keys whose corresponding value should be returned.

    Returns:
    str or list of str: 
        - If input is a string, returns the value associated with the key or "." if not found.
        - If input is a list, returns a list of values corresponding to each string entry in `x`.

    Raises:
    ValueError: If the input is neither a string nor a list of strings, or if any key-value pair is malformed.
    """

    key_values = dict(item.split("=") for item in x.split(";") if "=" in item)

    keys_of_interest = {k: key_values.get(k, None) for k in keys}
    
    return keys_of_interest



def write_log_errors(errors, log_dir):

    output_log = f"{log_dir}/errors.log"

    if os.path.exists(output_log):
        with open(output_log, "a") as log:
            log.write("\n".join(errors))
    else:
        with open(output_log, "w") as log:
            log.write("\n".join(errors))


def check_pass_error_logging(error_df, log_dir, col_mapping, error_message):

    errors = []
    errors.append(error_message)
    for _,row in error_df.iterrows():
        errors.append(f"{col_mapping['id_column']}: {row[col_mapping['id_column']]}")
        errors.append(f"\tSubstitutions: {row[col_mapping['mutation_column']]}")
        
        if len(row["Sequence Errors"]) > 0:
            errors.append(f"\tSequence Errors: {', '.join(row['Sequence Errors'])}")

        errors.append("\tMutation Errors:")
        for error in row["Mutation Errors"]:
            errors.append(f"    - {error}")
        errors.append("------------------------------------")

    write_log_errors(errors, log_dir)


def missing_sequence_error_logging(error_df, log_dir, col_mapping, error_message):
    errors = []
    errors.append(error_message)
    for _,row in error_df.iterrows():
        errors.append(f"{col_mapping['id_column']}: {row[col_mapping['id_column']]}")
        errors.append(f"\tSubstitutions: {row[col_mapping['mutation_column']]}")
        errors.append(f"\tSequence Error: No matching sequence found in FASTA file.")
        errors.append("------------------------------------")

    write_log_errors(errors, log_dir)

    
DRY_RUN_LEVEL = 25
def dry_run(self, message, *args, **kwargs):
    """Log a message at the DRY RUN level."""
    if self.isEnabledFor(DRY_RUN_LEVEL):
        self._log(DRY_RUN_LEVEL, message, args, **kwargs)


def set_logger_level(log_level):

    # Add Dry Run Logging level
    logging.addLevelName(DRY_RUN_LEVEL, "DRY RUN")
    logging.Logger.dry_run = dry_run

    # Get the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Check if handlers exist (to avoid adding multiple handlers)
    if not root_logger.hasHandlers():
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        root_logger.addHandler(handler)
    else:
        # Update the format of existing handlers
        for handler in root_logger.handlers:
            handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))


def create_directory(dir_path, dry_run, logged_status):
    """
    Create a directory if it does not already exist.

    Parameters:
    - dir_path (str): Path of the directory to be created.
    - dry_run (bool): If True, only logs the action without actually creating the directory.
    - logged_status (bool): Indicates if the action has already been logged.

    Returns:
    - tuple: Directory path and updated logged_status.
    """
    if not os.path.exists(f"{dir_path}"):
        if dry_run:
            # Log the action if it hasn't been logged yet
            if not logged_status:
                logger.dry_run(f"Would've created {dir_path}.")
                logged_status = True
        else:
            # Create the directory if dry_run is False
            os.makedirs(f"{dir_path}")
    
    return dir_path, logged_status


def configurations():
    """
    Load configuration settings from the config.yaml file.

    Returns:
    - dict: Configuration settings for the active location.
    """
    # Load the config.yaml file
    config_path = pkg_resources.files("mutpredpy").joinpath("config.yaml")
    with open(config_path, 'r') as config:
        configs = yaml.safe_load(config)

    # Get the current active location setting
    cur_location = configs.get('active_location', 'local')

    # Retrieve configuration for the active location
    location_config = configs['locations'].get(cur_location, {})

    return location_config


def mutpred2_path():
    """
    Retrieve the path for MutPred2 from the configuration file.

    Returns:
    - str: Path to the MutPred2 executable.

    Raises:
    - Exception: If the path is not defined in the configuration.
    """
    # Get configuration settings
    location_configs = configurations()

    # Retrieve the MutPred2 path
    mutpred_path = location_configs.get('mutpred2_path', "")

    # Raise an error if the path is not defined
    if not mutpred_path:
        raise Exception("No MutPred2 path designated in config.yaml.")
    elif not os.path.exists(mutpred_path):
        raise Exception(f"MutPred2 path {mutpred_path} not found. Change in config.yaml.")
    else:
        return mutpred_path


def catalog_directory():
    """
    Retrieve the path for the catalog directory from the configuration file.

    Returns:
    - str: Path to the catalog directory.

    Raises:
    - Exception: If the catalog path is not defined in the configuration.
    """
    # Get configuration settings
    location_configs = configurations()

    # Retrieve the catalog path
    catalog_path = location_configs.get("catalog", "")

    # Raise an error if the path is not defined
    if not catalog_path:
        raise Exception("No catalog path designated in config.yaml")
    elif not os.path.exists(catalog_path):
        raise Exception(f"Catalog path {catalog_path} not found. Change in config.yaml.")
    else:
        return catalog_path