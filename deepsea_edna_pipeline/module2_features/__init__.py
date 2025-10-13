"""
Module 2: Feature Extraction

This module handles advanced feature extraction from ASV sequences including
k-mer analysis, MinHash LSH embedding, and Frequency Chaos Game Representation (FCGR).
"""

from .pipeline_features import FeaturePipeline
from .kmer.kmer_extractor import KmerExtractor
from .kmer.minhash_lsh import MinHashLSH
from .fcgr.fcgr_generator import FCGRGenerator

__all__ = [
    'FeaturePipeline',
    'KmerExtractor',
    'MinHashLSH', 
    'FCGRGenerator'
]