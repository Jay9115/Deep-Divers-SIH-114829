"""
Chimera filtering using combined reference-based and de novo methods.
"""

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
from Bio import SeqIO


class ChimeraFilter:
    """
    Chimera detection and removal using VSEARCH.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize chimera filter.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
    def remove_chimeras(self,
                       input_fasta: Path,
                       output_dir: Path,
                       reference_db: Optional[Path] = None,
                       abundance_skew: float = 2.0,
                       threads: int = 1) -> Tuple[Path, Dict]:
        """
        Remove chimeric sequences using reference-based and de novo methods.
        
        Args:
            input_fasta: Input FASTA file with sequences
            output_dir: Output directory
            reference_db: Reference database for chimera detection (optional)
            abundance_skew: Abundance skew for de novo detection
            threads: Number of threads to use
            
        Returns:
            Tuple of (non_chimeric_output_path, chimera_stats)
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output files
        nonchimeric_output = output_dir / f"nonchimeric_{input_fasta.name}"
        chimeric_output = output_dir / f"chimeric_{input_fasta.name}"
        
        stats = {}
        
        # Step 1: Reference-based chimera detection (if reference provided)
        if reference_db and reference_db.exists():
            ref_nonchimeric = output_dir / f"ref_nonchimeric_{input_fasta.name}"
            ref_chimeric = output_dir / f"ref_chimeric_{input_fasta.name}"
            
            ref_stats = self._reference_chimera_detection(
                input_fasta, ref_nonchimeric, ref_chimeric, reference_db, threads
            )
            stats['reference_based'] = ref_stats
            
            # Use reference-filtered sequences for de novo detection
            denovo_input = ref_nonchimeric
        else:
            denovo_input = input_fasta
            stats['reference_based'] = {'skipped': True}
        
        # Step 2: De novo chimera detection
        denovo_stats = self._denovo_chimera_detection(
            denovo_input, nonchimeric_output, chimeric_output, 
            abundance_skew, threads
        )
        stats['de_novo'] = denovo_stats
        
        # Calculate overall statistics
        input_count = self._count_sequences(input_fasta)
        output_count = self._count_sequences(nonchimeric_output)
        
        stats['overall'] = {
            'input_sequences': input_count,
            'output_sequences': output_count,
            'chimeric_sequences': input_count - output_count,
            'chimera_rate': (input_count - output_count) / input_count if input_count > 0 else 0
        }
        
        self.logger.info(f"Chimera filtering completed. Retained {output_count}/{input_count} sequences")
        
        return nonchimeric_output, stats
    
    def _reference_chimera_detection(self,
                                   input_fasta: Path,
                                   nonchimeric_output: Path,
                                   chimeric_output: Path,
                                   reference_db: Path,
                                   threads: int) -> Dict:
        """
        Perform reference-based chimera detection.
        """
        self.logger.info("Running reference-based chimera detection")
        
        cmd = [
            'vsearch',
            '--uchime_ref', str(input_fasta),
            '--db', str(reference_db),
            '--nonchimeras', str(nonchimeric_output),
            '--chimeras', str(chimeric_output),
            '--threads', str(threads)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return self._parse_chimera_stats(result.stderr)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Reference-based chimera detection failed: {e.stderr}")
            raise
    
    def _denovo_chimera_detection(self,
                                input_fasta: Path,
                                nonchimeric_output: Path,
                                chimeric_output: Path,
                                abundance_skew: float,
                                threads: int) -> Dict:
        """
        Perform de novo chimera detection.
        """
        self.logger.info("Running de novo chimera detection")
        
        cmd = [
            'vsearch',
            '--uchime_denovo', str(input_fasta),
            '--nonchimeras', str(nonchimeric_output),
            '--chimeras', str(chimeric_output),
            '--abskew', str(abundance_skew),
            '--threads', str(threads)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return self._parse_chimera_stats(result.stderr)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"De novo chimera detection failed: {e.stderr}")
            raise
    
    def _parse_chimera_stats(self, stderr_output: str) -> Dict:
        """
        Parse chimera detection statistics from VSEARCH stderr.
        """
        stats = {}
        
        lines = stderr_output.split('\\n')
        
        for line in lines:
            line = line.strip()
            
            if 'sequences' in line and 'chimeric' in line:
                parts = line.split()
                stats['chimeric_sequences'] = int(parts[0])
            elif 'sequences' in line and 'non-chimeric' in line:
                parts = line.split()
                stats['nonchimeric_sequences'] = int(parts[0])
        
        return stats
    
    def _count_sequences(self, fasta_file: Path) -> int:
        """
        Count sequences in FASTA file.
        """
        try:
            return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        except:
            return 0