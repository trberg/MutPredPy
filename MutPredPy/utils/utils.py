import hashlib
import yaml
import re
import logging


# Configure logging for debugging and tracking
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


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
    SINGLE_AA = "ACDEFGHIKLMNPQRSTVWY" #UX=
    TRIPLETS_AA = "Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Sec|Thr|Trp|Tyr|Val"

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
    def mutation_mapping(x):
        """Automatically detect and convert a mutation string between single-letter and three-letter amino acid notation.

        Args:
            x (str): Mutation string (e.g., "A123T" or "Ala123Thr").

        Returns:
            str: Mutation in the opposite notation format.
        """
        try:
            mut = x #.replace("%3D","=")
            loc = re.search(r'\d+', mut).group()  # Extract numeric position
            ref, alt = mut.split(loc)  # Split reference and alternate residues

            # Detect if the input is single-letter or three-letter notation
            if len(ref) == 1:  # Single-letter input
                ref_mapped = AminoAcidMap.map_amino_acid(ref)
            elif len(ref) == 3:
                ref_mapped = AminoAcidMap.map_amino_acid(ref.capitalize())
            elif len(ref) == 0:
                ref_mapped = ""
            else:
                raise ValueError("Inconsistent mutation format detected.")

            if len(alt) == 1:
                alt_mapped = AminoAcidMap.map_amino_acid(alt)
            elif len(alt) == 3:  # Three-letter input
                alt_mapped = AminoAcidMap.map_amino_acid(alt.capitalize())
            elif len(alt) == 0:
                alt_mapped = ""
            else:
                raise ValueError("Inconsistent mutation format detected.")

            return f"{ref_mapped}{loc}{alt_mapped}"
        
        except (AttributeError, ValueError) as e:
            #logging.error(f"Invalid mutation format: '{mut}' - {e}")
            pass
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
