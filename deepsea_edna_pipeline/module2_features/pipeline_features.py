"""
Feature extraction pipeline orchestrator for Module 2.
"""

from pathlib import Path
from typing import Dict, List, Optional
import logging
import yaml
import time

from ..utils.logging_utils import create_run_logger, log_step_start, log_step_complete
from ..utils.io_utils import create_output_structure
from .kmer.kmer_extractor import KmerExtractor
from .kmer.minhash_lsh import MinHashLSHAnalyzer
from .fcgr.fcgr_generator import FCGRGenerator


class FeaturePipeline:
    """
    Complete feature extraction pipeline for ASV sequences.
    """
    
    def __init__(self, config_file: Optional[Path] = None):
        if config_file is None:
            config_file = Path(__file__).parent / "config.yaml"
        self.config = self._load_config(config_file)
        
    def _load_config(self, config_file: Path) -> Dict:
        """Load configuration from YAML file."""
        try:
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            return self._get_default_config()
    
    def _get_default_config(self) -> Dict:
        """Get default configuration."""
        return {
            'kmer': {
                'k': 6,
                'normalize': True,
                'n_jobs': 4
            },
            'fcgr': {
                'k': 8,
                'size': 256,
                'normalize': 'frequency'
            },
            'minhash': {
                'num_perm': 128,
                'threshold': 0.7
            }
        }
    
    def run_feature_extraction(self, input_fasta: Path, output_dir: Path,
                             run_id: Optional[str] = None) -> Dict:
        """
        Run complete feature extraction pipeline.
        
        Args:
            input_fasta: Input FASTA file with ASV sequences
            output_dir: Output directory
            run_id: Optional run identifier
            
        Returns:
            Dictionary with extraction results
        """
        if run_id is None:
            run_id = f"features_{int(time.time())}"
        
        # Set up output structure and logging
        output_structure = create_output_structure(output_dir, run_id)
        logger = create_run_logger(run_id, output_structure['logs'])
        
        logger.info(f"Starting feature extraction: {run_id}")
        
        results = {}
        
        # Extract k-mer features
        logger.info("=== K-mer Feature Extraction ===")
        kmer_config = self.config['kmer']
        kmer_extractor = KmerExtractor(k=kmer_config['k'], logger=logger)
        
        features, seq_ids = kmer_extractor.extract_features_from_fasta(
            input_fasta, 
            normalize=kmer_config['normalize'],
            n_jobs=kmer_config['n_jobs']
        )
        
        # Save k-mer features
        kmer_output = output_structure['features'] / "kmer_features.npz"
        kmer_extractor.save_features(features, seq_ids, kmer_output)
        results['kmer_features'] = str(kmer_output)
        
        # Generate FCGR representations
        logger.info("=== FCGR Generation ===")
        fcgr_config = self.config['fcgr']
        fcgr_generator = FCGRGenerator(
            k=fcgr_config['k'], 
            size=fcgr_config['size'], 
            logger=logger
        )
        
        fcgr_generator.batch_generate_fcgr(
            input_fasta,
            output_structure['features'] / "fcgr",
            normalize=fcgr_config['normalize']
        )
        results['fcgr_features'] = str(output_structure['features'] / "fcgr" / "fcgr_features.npz")
        
        # Build MinHash LSH index
        logger.info("=== MinHash LSH Index ===")
        minhash_config = self.config['minhash']
        lsh_analyzer = MinHashLSHAnalyzer(
            num_perm=minhash_config['num_perm'],
            threshold=minhash_config['threshold'],
            logger=logger
        )
        
        # Add sequences to LSH index
        from Bio import SeqIO
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_str = str(record.seq).upper()
            kmers = set(seq_str[i:i+kmer_config['k']] 
                       for i in range(len(seq_str) - kmer_config['k'] + 1))
            lsh_analyzer.add_sequence(record.id, kmers)
        
        # Save LSH index
        lsh_output = output_structure['features'] / "minhash_lsh.pkl"
        lsh_analyzer.save_index(lsh_output)
        results['lsh_index'] = str(lsh_output)
        
        logger.info("Feature extraction completed successfully")
        
        return results