"""
Initialization module for the fasta package.

This module imports essential functions and classes for handling and processing 
FASTA sequences in MutPredPy and defines the public API of the package.
"""

from .fasta import (
    check_sequences,
    clean_fasta_sequence,
    alignment_score,
    read_mutpred_input_fasta,
    Protein,
    read_fasta,
)

__all__ = [
    "check_sequences",
    "clean_fasta_sequence",
    "alignment_score",
    "read_mutpred_input_fasta",
    "Protein",
    "read_fasta",
]
