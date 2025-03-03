import pytest # type: ignore
import pandas as pd
import os
import tempfile
from unittest.mock import patch, MagicMock
import importlib.resources as pkg_resources

from MutPredPy.prep import Prepare

# Load test datasets from test_data directory
TEST_DATA_DIR = "tests/test_data"

def load_test_data(filename):
    return pkg_resources.files("mutpredpy").joinpath(f"tests/test_data/{filename}")


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield tmpdirname


# Test Prepare Class
'''
@patch("prep.u.create_directory")
def test_prepare_class(mock_create_dir, temp_dir):
    mock_create_dir.return_value = (temp_dir, True)
    prepare = Prepare(input="test_input.txt", working_dir=temp_dir, time=24, dry_run=False, 
                      canonical=False, all_possible=False, verbose=True, users=1, fasta="")
    assert prepare.get_working_dir() == os.path.abspath(temp_dir)

# Test Prepare Methods Using Test Datasets
def test_prepare_with_datasets(temp_dir):
    test_files = [
        "test.ENSP_ENST_ENSG.txt", "test.ENSP_unversioned.txt", "test.ENSP.txt",
        "test.refseq_unversioned.txt", "test.refseq.txt", "test.vep_110.txt"
    ]
    
    for test_file in test_files:
        prepare = Prepare(input=test_file, working_dir=temp_dir, time=24, dry_run=False, 
                          canonical=False, all_possible=False, verbose=True, users=1, fasta="")
        assert not prepare.get_input().empty
'''
def test_versioned_ensids(temp_dir):
    test_file = "test.ENSP_ENST_ENSG.txt"
    
    prepare = Prepare(input=load_test_data(test_file), working_dir=temp_dir, time=24, dry_run=True, 
                          canonical=False, all_possible=False, verbose=True, users=1, fasta="")
    assert not prepare.get_input().empty
    assert prepare.get_base() == "test.ENSP_ENST_ENSG"
'''
@patch("prep.Prepare.add_sequences")
def test_prepare_add_sequences(mock_add_sequences, temp_dir):
    test_files = [
        "test.ENSP_ENST_ENSG.txt", "test.ENSP_unversioned.txt", "test.ENSP.txt",
        "test.refseq_unversioned.txt", "test.refseq.txt", "test.vep_110.txt"
    ]
    
    for test_file in test_files:
        df = load_test_data(test_file)
        mock_add_sequences.return_value = (df, {"mutation_column": "Substitution"})
        prepare = Prepare(input=test_file, working_dir=temp_dir, time=24, dry_run=False, 
                          canonical=False, all_possible=False, verbose=True, users=1, fasta="")
        data, col_mapping = prepare.add_sequences(df, {})
        assert "Substitution" in data.columns
'''

# Run the tests if executed directly
if __name__ == "__main__":
    pytest.main()
