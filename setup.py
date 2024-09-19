"""Setup"""
import os
from setuptools import setup

# attr: MutPredPy.__version__

# figure out the version
#about = {}
#here = os.path.abspath(os.path.dirname(__file__))
#with open(os.path.join(here, "MutPredPy", "__version__.py")) as f:
#    exec(f.read(), about)
# HACK: this must be kept because __init__ imports the discussion
# modules which import requests which has to be installed first.

setup(name="MutPredPy", version="1.1.0", package_data={"MutPredPy":['MutPredPy/resources/Homo_sapiens.GRCh38.combined.pep.all.fa','MutPredPy/resources/Homo_sapiens.GRCh37.combined.pep.all.fa','MutPredPy/resources/sequence_time.npy','MutPredPy/resources/memory_usage.npy']})
