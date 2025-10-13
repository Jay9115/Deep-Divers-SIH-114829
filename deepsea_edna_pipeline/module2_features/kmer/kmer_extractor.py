"""
K-mer feature extraction for sequence analysis.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from collections import defaultdict, Counter
import logging
from Bio import SeqIO
import joblib
from itertools import product


class KmerExtractor:
    """
    Extract k-mer features from DNA sequences for machine learning applications.
    """
    
    def __init__(self, k: int = 6, logger: Optional[logging.Logger] = None):
        """
        Initialize k-mer extractor.
        
        Args:
            k: K-mer size (default: 6)
            logger: Logger instance
        """
        self.k = k
        self.logger = logger or logging.getLogger(__name__)
        self.alphabet = ['A', 'C', 'G', 'T']
        
        # Generate all possible k-mers
        self.all_kmers = [''.join(kmer) for kmer in product(self.alphabet, repeat=k)]
        self.kmer_to_index = {kmer: i for i, kmer in enumerate(self.all_kmers)}
        
        self.logger.info(f"Initialized k-mer extractor with k={k} ({len(self.all_kmers)} possible k-mers)")
    
    def extract_kmers_from_sequence(self, sequence: str, normalize: bool = True) -> np.ndarray:
        """
        Extract k-mer frequencies from a single sequence.
        
        Args:
            sequence: DNA sequence string
            normalize: Whether to normalize frequencies
            
        Returns:
            K-mer frequency vector
        """
        # Clean sequence - keep only valid nucleotides
        clean_seq = ''.join(c.upper() for c in sequence if c.upper() in self.alphabet)
        
        if len(clean_seq) < self.k:
            # Return zero vector if sequence too short
            return np.zeros(len(self.all_kmers))
        
        # Count k-mers
        kmer_counts = np.zeros(len(self.all_kmers))
        
        for i in range(len(clean_seq) - self.k + 1):
            kmer = clean_seq[i:i + self.k]
            if kmer in self.kmer_to_index:
                kmer_counts[self.kmer_to_index[kmer]] += 1
        
        if normalize and np.sum(kmer_counts) > 0:
            kmer_counts = kmer_counts / np.sum(kmer_counts)
        
        return kmer_counts
    
    def extract_features_from_fasta(self, fasta_path: Path, 
                                  normalize: bool = True,
                                  n_jobs: int = 1) -> Tuple[np.ndarray, List[str]]:
        """
        Extract k-mer features from all sequences in a FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            normalize: Whether to normalize frequencies  
            n_jobs: Number of parallel jobs (-1 for all cores)
            
        Returns:
            Tuple of (feature_matrix, sequence_ids)
        """
        self.logger.info(f"Extracting k-mer features from {fasta_path}")
        
        # Read sequences
        sequences = []
        seq_ids = []
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences.append(str(record.seq))
            seq_ids.append(record.id)
        
        if not sequences:
            raise ValueError(f"No sequences found in {fasta_path}")
        
        self.logger.info(f"Processing {len(sequences)} sequences")
        
        # Extract features in parallel
        if n_jobs == 1:
            # Sequential processing
            features = [self.extract_kmers_from_sequence(seq, normalize) 
                       for seq in sequences]
        else:
            # Parallel processing
            features = joblib.Parallel(n_jobs=n_jobs)(
                joblib.delayed(self.extract_kmers_from_sequence)(seq, normalize)
                for seq in sequences
            )
        
        feature_matrix = np.vstack(features)
        
        self.logger.info(f"Extracted feature matrix: {feature_matrix.shape}")
        
        return feature_matrix, seq_ids
    
    def extract_kmer_spectrum(self, sequences: List[str]) -> Dict[str, int]:
        """
        Extract k-mer spectrum (all k-mers and their counts) from sequences.
        
        Args:
            sequences: List of DNA sequences
            
        Returns:
            Dictionary mapping k-mers to their total counts
        """
        kmer_counts = Counter()
        
        for sequence in sequences:
            clean_seq = ''.join(c.upper() for c in sequence if c.upper() in self.alphabet)
            
            for i in range(len(clean_seq) - self.k + 1):
                kmer = clean_seq[i:i + self.k]
                if len(kmer) == self.k and all(c in self.alphabet for c in kmer):
                    kmer_counts[kmer] += 1
        
        return dict(kmer_counts)
    
    def calculate_kmer_diversity(self, kmer_counts: Dict[str, int]) -> Dict[str, float]:
        """
        Calculate k-mer diversity metrics.
        
        Args:
            kmer_counts: Dictionary of k-mer counts
            
        Returns:
            Dictionary with diversity metrics
        """
        if not kmer_counts:
            return {'richness': 0, 'shannon': 0, 'simpson': 0, 'evenness': 0}
        
        counts = np.array(list(kmer_counts.values()))
        total = np.sum(counts)
        
        if total == 0:
            return {'richness': 0, 'shannon': 0, 'simpson': 0, 'evenness': 0}
        
        # Richness (number of unique k-mers)
        richness = len(kmer_counts)
        
        # Shannon diversity
        proportions = counts / total
        shannon = -np.sum(proportions * np.log(proportions + 1e-10))
        
        # Simpson diversity
        simpson = 1 - np.sum(proportions ** 2)
        
        # Evenness
        evenness = shannon / np.log(richness) if richness > 1 else 0
        
        return {
            'richness': richness,
            'shannon': shannon,
            'simpson': simpson,
            'evenness': evenness
        }
    
    def save_features(self, feature_matrix: np.ndarray, seq_ids: List[str],
                     output_path: Path, format: str = 'numpy') -> None:
        """
        Save extracted features to file.
        
        Args:
            feature_matrix: Feature matrix (samples x features)
            seq_ids: Sequence identifiers
            output_path: Output file path
            format: Output format ('numpy', 'csv', 'hdf5')
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if format == 'numpy':
            # Save as compressed numpy arrays
            np.savez_compressed(output_path, 
                              features=feature_matrix, 
                              seq_ids=seq_ids,
                              kmer_names=self.all_kmers)
            
        elif format == 'csv':
            # Save as CSV file
            import pandas as pd
            df = pd.DataFrame(feature_matrix, 
                            index=seq_ids, 
                            columns=self.all_kmers)
            df.to_csv(output_path)
            
        elif format == 'hdf5':
            # Save as HDF5 file
            import h5py
            with h5py.File(output_path, 'w') as f:
                f.create_dataset('features', data=feature_matrix)
                f.create_dataset('seq_ids', data=[s.encode() for s in seq_ids])
                f.create_dataset('kmer_names', data=[k.encode() for k in self.all_kmers])
                f.attrs['k'] = self.k
                f.attrs['num_sequences'] = len(seq_ids)
                f.attrs['num_features'] = len(self.all_kmers)
        
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        self.logger.info(f"Saved k-mer features to {output_path} (format: {format})")
    
    def load_features(self, input_path: Path) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Load previously saved k-mer features.
        
        Args:
            input_path: Path to saved features file
            
        Returns:
            Tuple of (feature_matrix, seq_ids, kmer_names)
        """
        if input_path.suffix == '.npz':
            # Load numpy format
            data = np.load(input_path)
            return data['features'], data['seq_ids'].tolist(), data['kmer_names'].tolist()
            
        elif input_path.suffix == '.csv':
            # Load CSV format
            import pandas as pd
            df = pd.read_csv(input_path, index_col=0)
            return df.values, df.index.tolist(), df.columns.tolist()
            
        elif input_path.suffix in ['.h5', '.hdf5']:
            # Load HDF5 format
            import h5py
            with h5py.File(input_path, 'r') as f:
                features = f['features'][:]
                seq_ids = [s.decode() for s in f['seq_ids'][:]]
                kmer_names = [k.decode() for k in f['kmer_names'][:]]
                return features, seq_ids, kmer_names
        
        else:
            raise ValueError(f"Unsupported file format: {input_path.suffix}")
    
    def compare_kmer_profiles(self, profile1: np.ndarray, profile2: np.ndarray,
                            metric: str = 'cosine') -> float:
        """
        Compare two k-mer profiles using specified distance metric.
        
        Args:
            profile1: First k-mer profile
            profile2: Second k-mer profile  
            metric: Distance metric ('cosine', 'euclidean', 'manhattan', 'js_divergence')
            
        Returns:
            Distance/similarity value
        """
        from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances, manhattan_distances
        
        if metric == 'cosine':
            return cosine_similarity([profile1], [profile2])[0, 0]
        elif metric == 'euclidean':
            return euclidean_distances([profile1], [profile2])[0, 0]
        elif metric == 'manhattan':
            return manhattan_distances([profile1], [profile2])[0, 0]
        elif metric == 'js_divergence':
            # Jensen-Shannon divergence
            def js_divergence(p, q):
                p = p + 1e-10  # Add small constant to avoid log(0)
                q = q + 1e-10
                m = 0.5 * (p + q)
                return 0.5 * np.sum(p * np.log(p / m)) + 0.5 * np.sum(q * np.log(q / m))
            
            # Normalize profiles
            p1 = profile1 / np.sum(profile1) if np.sum(profile1) > 0 else profile1
            p2 = profile2 / np.sum(profile2) if np.sum(profile2) > 0 else profile2
            
            return js_divergence(p1, p2)
        else:
            raise ValueError(f"Unknown metric: {metric}")