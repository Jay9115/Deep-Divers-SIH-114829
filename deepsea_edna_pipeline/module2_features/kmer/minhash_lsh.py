"""
MinHash LSH implementation for scalable sequence similarity search.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import logging
from datasketch import MinHashLSH, MinHash
import mmh3
import pickle


class MinHashLSHAnalyzer:
    """
    MinHash LSH for efficient similarity search in large sequence datasets.
    """
    
    def __init__(self, num_perm: int = 128, threshold: float = 0.7,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize MinHash LSH analyzer.
        
        Args:
            num_perm: Number of permutations for MinHash
            threshold: Jaccard similarity threshold for LSH
            logger: Logger instance
        """
        self.num_perm = num_perm
        self.threshold = threshold
        self.logger = logger or logging.getLogger(__name__)
        
        self.lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
        self.minhashes = {}
        
    def create_minhash(self, kmers: Set[str]) -> MinHash:
        """Create MinHash signature from k-mer set."""
        minhash = MinHash(num_perm=self.num_perm)
        for kmer in kmers:
            minhash.update(kmer.encode('utf-8'))
        return minhash
    
    def add_sequence(self, seq_id: str, kmers: Set[str]) -> None:
        """Add sequence to LSH index."""
        minhash = self.create_minhash(kmers)
        self.lsh.insert(seq_id, minhash)
        self.minhashes[seq_id] = minhash
        
    def query_similar(self, kmers: Set[str]) -> List[str]:
        """Find similar sequences using LSH."""
        query_minhash = self.create_minhash(kmers)
        return list(self.lsh.query(query_minhash))
    
    def estimate_jaccard(self, seq_id1: str, seq_id2: str) -> float:
        """Estimate Jaccard similarity between two sequences."""
        if seq_id1 in self.minhashes and seq_id2 in self.minhashes:
            return self.minhashes[seq_id1].jaccard(self.minhashes[seq_id2])
        return 0.0
    
    def save_index(self, output_path: Path) -> None:
        """Save LSH index to file."""
        with open(output_path, 'wb') as f:
            pickle.dump({
                'lsh': self.lsh,
                'minhashes': self.minhashes,
                'num_perm': self.num_perm,
                'threshold': self.threshold
            }, f)
    
    def load_index(self, input_path: Path) -> None:
        """Load LSH index from file."""
        with open(input_path, 'rb') as f:
            data = pickle.load(f)
            self.lsh = data['lsh']
            self.minhashes = data['minhashes']
            self.num_perm = data['num_perm'] 
            self.threshold = data['threshold']