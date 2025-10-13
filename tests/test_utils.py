"""
Test suite for utility functions.
"""

import unittest
from pathlib import Path
import tempfile
import logging

from deepsea_edna.utils.logging_utils import setup_logger
from deepsea_edna.utils.system_check import check_command_availability
from deepsea_edna.utils.parallel import ParallelProcessor


class TestUtils(unittest.TestCase):
    """Test utility functions."""
    
    def test_logger_setup(self):
        """Test logger initialization."""
        logger = setup_logger('test_logger', console_output=True)
        
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, 'test_logger')
    
    def test_command_availability(self):
        """Test command availability check."""
        # Test with a command that should exist
        python_available = check_command_availability('python')
        self.assertIsInstance(python_available, bool)
        
        # Test with a command that shouldn't exist
        fake_available = check_command_availability('nonexistent_command_12345')
        self.assertFalse(fake_available)
    
    def test_parallel_processor(self):
        """Test parallel processing functionality."""
        processor = ParallelProcessor(n_processes=2)
        
        # Test simple function mapping
        def square(x):
            return x ** 2
        
        results = processor.map_function(square, [1, 2, 3, 4])
        expected = [1, 4, 9, 16]
        
        self.assertEqual(results, expected)


if __name__ == '__main__':
    unittest.main()