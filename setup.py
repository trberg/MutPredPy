"""Setup"""
import os
from setuptools import setup


# figure out the version
about = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "MutPredPy", "__version__.py")) as f:
    exec(f.read(), about)
# HACK: this must be kept because __init__ imports the discussion
# modules which import requests which has to be installed first.
setup(version=about["__version__"], package_data={'MutPredPy':['resources/Homo_sapiens.GRCh38.combined.pep.all.fa','resources/Homo_sapiens.GRCh37.combined.pep.all.fa','resources/sequence_time.npy','resources/memory_usage.npy']})