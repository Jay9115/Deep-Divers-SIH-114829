"""
Test suite for Module 1: QC and ASV generation.
"""

import unittest
from pathlib import Path
import tempfile
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from deepsea_edna.module1_qc_asv.scripts.marker_detect import MarkerDetector
from deepsea_edna.module1_qc_asv.scripts.cutadapt_wrapper import CutadaptWrapper
from deepsea_edna.utils.qc_metrics import calculate_basic_stats


class TestModule1QCAsv(unittest.TestCase):
    """Test cases for Module 1 components."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.test_sequences = [
            SeqRecord(Seq("ATCGATCGATCG"), id="seq1", description="test sequence 1"),
            SeqRecord(Seq("GCTAGCTAGCTA"), id="seq2", description="test sequence 2")
        ]
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_marker_detector(self):
        """Test marker detection functionality."""
        detector = MarkerDetector()
        
        # Test with mock sequences
        detected_marker, stats = detector.detect_marker(self.test_sequences)
        
        self.assertIsInstance(detected_marker, str)
        self.assertIsInstance(stats, dict)
    
    def test_basic_stats_calculation(self):
        """Test basic sequence statistics calculation."""
        stats = calculate_basic_stats(self.test_sequences)
        
        self.assertEqual(stats['total_sequences'], 2)
        self.assertEqual(stats['mean_length'], 12.0)
        self.assertEqual(stats['min_length'], 12)
        self.assertEqual(stats['max_length'], 12)


if __name__ == '__main__':
    unittest.main()