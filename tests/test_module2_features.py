"""
Test suite for Module 2: Feature extraction.
"""

import unittest
import numpy as np
from pathlib import Path
import tempfile

from deepsea_edna.module2_features.kmer.kmer_extractor import KmerExtractor
from deepsea_edna.module2_features.fcgr.fcgr_generator import FCGRGenerator


class TestModule2Features(unittest.TestCase):
    """Test cases for Module 2 components."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_sequence = "ATCGATCGATCGATCGATCGATCG"
        self.kmer_extractor = KmerExtractor(k=4)
        self.fcgr_generator = FCGRGenerator(k=4, size=64)
    
    def test_kmer_extraction(self):
        """Test k-mer extraction from sequence."""
        kmer_vector = self.kmer_extractor.extract_kmers_from_sequence(
            self.test_sequence, normalize=True
        )
        
        self.assertIsInstance(kmer_vector, np.ndarray)
        self.assertEqual(len(kmer_vector), 4**4)  # 4^k possible k-mers
        self.assertAlmostEqual(np.sum(kmer_vector), 1.0, places=6)  # Normalized
    
    def test_fcgr_generation(self):
        """Test FCGR generation from sequence."""
        fcgr = self.fcgr_generator.sequence_to_fcgr(self.test_sequence)
        
        self.assertIsInstance(fcgr, np.ndarray)
        self.assertEqual(fcgr.shape, (64, 64))
        self.assertGreaterEqual(np.sum(fcgr), 0)  # Non-negative values


class TestUtils(unittest.TestCase):
    """Test utility functions."""
    
    def test_kmer_diversity_calculation(self):
        """Test k-mer diversity metrics."""
        from deepsea_edna.module2_features.kmer.kmer_extractor import KmerExtractor
        
        extractor = KmerExtractor(k=3)
        kmer_counts = {'ATG': 10, 'GCT': 5, 'TAA': 3}
        
        diversity = extractor.calculate_kmer_diversity(kmer_counts)
        
        self.assertIn('richness', diversity)
        self.assertIn('shannon', diversity)
        self.assertEqual(diversity['richness'], 3)


if __name__ == '__main__':
    unittest.main()