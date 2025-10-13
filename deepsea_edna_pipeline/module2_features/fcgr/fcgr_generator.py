"""
Frequency Chaos Game Representation (FCGR) generator for DNA sequences.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
from Bio import SeqIO
import matplotlib.pyplot as plt


class FCGRGenerator:
    """
    Generate Frequency Chaos Game Representation (FCGR) from DNA sequences.
    FCGR provides a 2D image representation of sequence composition patterns.
    """
    
    def __init__(self, k: int = 8, size: int = 256, 
                 logger: Optional[logging.Logger] = None):
        """
        Initialize FCGR generator.
        
        Args:
            k: K-mer size for FCGR (determines resolution)
            size: Output image size (size x size)
            logger: Logger instance
        """
        self.k = k
        self.size = size
        self.logger = logger or logging.getLogger(__name__)
        
        # Define nucleotide coordinates in unit square
        self.coordinates = {
            'A': (0.0, 0.0),
            'C': (1.0, 0.0), 
            'G': (1.0, 1.0),
            'T': (0.0, 1.0)
        }
        
        # Calculate scale factor
        self.scale = size / (2 ** k)
        
    def sequence_to_fcgr(self, sequence: str) -> np.ndarray:
        """
        Convert DNA sequence to FCGR representation.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            2D numpy array representing FCGR
        """
        # Clean sequence
        clean_seq = ''.join(c.upper() for c in sequence if c.upper() in 'ACGT')
        
        if len(clean_seq) < self.k:
            return np.zeros((self.size, self.size))
        
        # Initialize FCGR matrix
        fcgr = np.zeros((self.size, self.size))
        
        # Process k-mers
        for i in range(len(clean_seq) - self.k + 1):
            kmer = clean_seq[i:i + self.k]
            
            # Calculate position using chaos game rules
            x, y = 0.5, 0.5  # Start at center
            
            for nucleotide in kmer:
                if nucleotide in self.coordinates:
                    nx, ny = self.coordinates[nucleotide]
                    # Move halfway to nucleotide corner
                    x = (x + nx) / 2
                    y = (y + ny) / 2
            
            # Convert to matrix indices
            row = int(y * (self.size - 1))
            col = int(x * (self.size - 1))
            
            # Increment frequency
            fcgr[row, col] += 1
        
        return fcgr
    
    def generate_fcgr_from_fasta(self, fasta_path: Path) -> Dict[str, np.ndarray]:
        """
        Generate FCGR representations for all sequences in FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to FCGR arrays
        """
        self.logger.info(f"Generating FCGR representations from {fasta_path}")
        
        fcgr_dict = {}
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            fcgr = self.sequence_to_fcgr(str(record.seq))
            fcgr_dict[record.id] = fcgr
        
        self.logger.info(f"Generated FCGR for {len(fcgr_dict)} sequences")
        
        return fcgr_dict
    
    def normalize_fcgr(self, fcgr: np.ndarray, method: str = 'frequency') -> np.ndarray:
        """
        Normalize FCGR representation.
        
        Args:
            fcgr: FCGR array
            method: Normalization method ('frequency', 'log', 'sqrt')
            
        Returns:
            Normalized FCGR array
        """
        if method == 'frequency':
            # Normalize to frequencies
            total = np.sum(fcgr)
            return fcgr / total if total > 0 else fcgr
            
        elif method == 'log':
            # Log transformation
            return np.log1p(fcgr)
            
        elif method == 'sqrt':
            # Square root transformation
            return np.sqrt(fcgr)
            
        else:
            raise ValueError(f"Unknown normalization method: {method}")
    
    def save_fcgr_image(self, fcgr: np.ndarray, output_path: Path,
                       cmap: str = 'viridis') -> None:
        """
        Save FCGR as image file.
        
        Args:
            fcgr: FCGR array
            output_path: Output image path
            cmap: Matplotlib colormap name
        """
        plt.figure(figsize=(8, 8))
        plt.imshow(fcgr, cmap=cmap, origin='lower')
        plt.title('FCGR Representation')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    
    def batch_generate_fcgr(self, input_fasta: Path, output_dir: Path,
                           save_images: bool = False, 
                           normalize: str = 'frequency') -> None:
        """
        Generate FCGR for all sequences and save results.
        
        Args:
            input_fasta: Input FASTA file
            output_dir: Output directory
            save_images: Whether to save individual FCGR images
            normalize: Normalization method
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate FCGRs
        fcgr_dict = self.generate_fcgr_from_fasta(input_fasta)
        
        # Normalize if requested
        if normalize:
            fcgr_dict = {seq_id: self.normalize_fcgr(fcgr, normalize)
                        for seq_id, fcgr in fcgr_dict.items()}
        
        # Save as numpy array collection
        fcgr_array = np.stack(list(fcgr_dict.values()))
        seq_ids = list(fcgr_dict.keys())
        
        np.savez_compressed(output_dir / "fcgr_features.npz",
                           fcgr=fcgr_array,
                           seq_ids=seq_ids,
                           k=self.k,
                           size=self.size)
        
        # Save individual images if requested
        if save_images:
            image_dir = output_dir / "fcgr_images"
            image_dir.mkdir(exist_ok=True)
            
            for seq_id, fcgr in fcgr_dict.items():
                self.save_fcgr_image(fcgr, image_dir / f"{seq_id}_fcgr.png")
        
        self.logger.info(f"Saved FCGR features to {output_dir}")
    
    def compare_fcgr(self, fcgr1: np.ndarray, fcgr2: np.ndarray,
                    metric: str = 'euclidean') -> float:
        """
        Compare two FCGR representations.
        
        Args:
            fcgr1: First FCGR array
            fcgr2: Second FCGR array
            metric: Distance metric ('euclidean', 'cosine', 'correlation')
            
        Returns:
            Distance/similarity value
        """
        # Flatten arrays
        f1 = fcgr1.flatten()
        f2 = fcgr2.flatten()
        
        if metric == 'euclidean':
            return np.linalg.norm(f1 - f2)
        elif metric == 'cosine':
            dot_product = np.dot(f1, f2)
            norms = np.linalg.norm(f1) * np.linalg.norm(f2)
            return 1 - (dot_product / norms) if norms > 0 else 1
        elif metric == 'correlation':
            return np.corrcoef(f1, f2)[0, 1] if len(f1) > 1 else 0
        else:
            raise ValueError(f"Unknown metric: {metric}")