"""
DeepSea eDNA Processing Pipeline

A comprehensive pipeline for processing environmental DNA sequences from deep-sea samples.
Includes quality control, ASV generation, and advanced feature extraction.
"""

__version__ = "0.1.0"
__author__ = "Deep Divers Team"
__email__ = "team@deepdivers.org"

from . import utils
from . import module1_qc_asv
from . import module2_features
from . import cli

__all__ = [
    "utils",
    "module1_qc_asv", 
    "module2_features",
    "cli"
]