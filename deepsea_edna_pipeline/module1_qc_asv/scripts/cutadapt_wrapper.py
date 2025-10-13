"""
Cutadapt wrapper for primer trimming and adapter removal.
"""

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging
import tempfile
from ..utils.logging_utils import log_step_start, log_step_complete


class CutadaptWrapper:
    """
    Wrapper for cutadapt tool to handle primer trimming and adapter removal.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize cutadapt wrapper.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
    def trim_primers(self, 
                    forward_reads: Path,
                    reverse_reads: Optional[Path],
                    forward_primer: str,
                    reverse_primer: str,
                    output_dir: Path,
                    min_length: int = 50,
                    max_error_rate: float = 0.1,
                    min_overlap: int = 3,
                    cores: int = 1) -> Tuple[Path, Optional[Path], Dict]:
        """
        Trim primers from paired-end reads using cutadapt.
        
        Args:
            forward_reads: Path to forward reads FASTQ file
            reverse_reads: Path to reverse reads FASTQ file (None for single-end)
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence  
            output_dir: Output directory
            min_length: Minimum read length after trimming
            max_error_rate: Maximum error rate for primer matching
            min_overlap: Minimum overlap between primer and read
            cores: Number of CPU cores to use
            
        Returns:
            Tuple of (forward_output_path, reverse_output_path, stats_dict)
        """
        log_step_start(self.logger, "primer_trimming", 
                      forward_reads=str(forward_reads),
                      reverse_reads=str(reverse_reads) if reverse_reads else None,
                      forward_primer=forward_primer,
                      reverse_primer=reverse_primer)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output file paths
        output_f = output_dir / f"trimmed_{forward_reads.stem}.fastq"
        output_r = None
        
        if reverse_reads:
            output_r = output_dir / f"trimmed_{reverse_reads.stem}.fastq"
        
        # Build cutadapt command
        cmd = [
            'cutadapt',
            '-g', forward_primer,  # 5' adapter (forward primer)
            '-a', self._reverse_complement(reverse_primer),  # 3' adapter (reverse of reverse primer)
            '--error-rate', str(max_error_rate),
            '--minimum-length', str(min_length), 
            '--overlap', str(min_overlap),
            '--cores', str(cores),
            '--output', str(output_f)
        ]
        
        if reverse_reads:
            # Paired-end mode
            cmd.extend([
                '-G', reverse_primer,  # 5' adapter for reverse read
                '-A', self._reverse_complement(forward_primer),  # 3' adapter for reverse read
                '--paired-output', str(output_r),
                str(forward_reads),
                str(reverse_reads)
            ])
        else:
            # Single-end mode
            cmd.append(str(forward_reads))
        
        # Run cutadapt
        self.logger.info(f"Running cutadapt command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse statistics from cutadapt output
            stats = self._parse_cutadapt_stats(result.stderr)
            
            log_step_complete(self.logger, "primer_trimming", **stats)
            
            return output_f, output_r, stats
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Cutadapt failed: {e.stderr}")
            raise
    
    def _reverse_complement(self, sequence: str) -> str:
        """
        Generate reverse complement of DNA sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            Reverse complement sequence
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                     'N': 'N', 'Y': 'R', 'R': 'Y', 'W': 'W',
                     'S': 'S', 'K': 'M', 'M': 'K', 'B': 'V',
                     'V': 'B', 'D': 'H', 'H': 'D'}
        
        return ''.join(complement.get(base.upper(), base) 
                      for base in reversed(sequence))
    
    def _parse_cutadapt_stats(self, stderr_output: str) -> Dict:
        """
        Parse statistics from cutadapt stderr output.
        
        Args:
            stderr_output: Cutadapt stderr text
            
        Returns:
            Dictionary with parsed statistics
        """
        stats = {}
        
        lines = stderr_output.split('\\n')
        
        for line in lines:
            if 'Total reads processed:' in line:
                stats['total_reads'] = int(line.split(':')[1].strip().replace(',', ''))
            elif 'Reads with adapters:' in line:
                stats['reads_with_adapters'] = int(line.split(':')[1].strip().split()[0].replace(',', ''))
            elif 'Reads that were too short:' in line:
                stats['reads_too_short'] = int(line.split(':')[1].strip().replace(',', ''))
            elif 'Reads written (passing filters):' in line:
                stats['reads_written'] = int(line.split(':')[1].strip().split()[0].replace(',', ''))
        
        return stats
    
    def remove_adapters(self,
                       input_files: Union[Path, List[Path]],
                       adapters: List[str], 
                       output_dir: Path,
                       **kwargs) -> Tuple[Union[Path, List[Path]], Dict]:
        """
        Remove sequencing adapters from reads.
        
        Args:
            input_files: Input FASTQ file(s) 
            adapters: List of adapter sequences to remove
            output_dir: Output directory
            **kwargs: Additional cutadapt parameters
            
        Returns:
            Tuple of (output_file_paths, statistics)
        """
        log_step_start(self.logger, "adapter_removal", 
                      input_files=str(input_files),
                      adapters=adapters)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if isinstance(input_files, Path):
            input_files = [input_files]
        
        output_files = []
        all_stats = {}
        
        for input_file in input_files:
            output_file = output_dir / f"cleaned_{input_file.name}"
            
            cmd = ['cutadapt']
            
            # Add adapters
            for adapter in adapters:
                cmd.extend(['-a', adapter])
            
            # Add standard parameters
            cmd.extend([
                '--minimum-length', str(kwargs.get('min_length', 30)),
                '--quality-cutoff', str(kwargs.get('quality_cutoff', 20)),
                '--cores', str(kwargs.get('cores', 1)),
                '--output', str(output_file),
                str(input_file)
            ])
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                stats = self._parse_cutadapt_stats(result.stderr)
                all_stats[input_file.name] = stats
                output_files.append(output_file)
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Adapter removal failed for {input_file}: {e.stderr}")
                raise
        
        log_step_complete(self.logger, "adapter_removal", **all_stats)
        
        if len(output_files) == 1:
            return output_files[0], all_stats
        else:
            return output_files, all_stats