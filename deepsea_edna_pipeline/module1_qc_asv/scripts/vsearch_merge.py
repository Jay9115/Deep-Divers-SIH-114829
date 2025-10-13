"""
VSEARCH wrapper for merging paired-end reads.
"""

import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple
import logging


class VSearchMerger:
    """
    Wrapper for VSEARCH tool to merge paired-end reads.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize VSEARCH merger.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
    def merge_pairs(self,
                   forward_reads: Path,
                   reverse_reads: Path, 
                   output_dir: Path,
                   min_overlap: int = 10,
                   max_diffs: int = 5,
                   max_length: int = 600,
                   threads: int = 1,
                   **kwargs) -> Tuple[Path, Dict]:
        """
        Merge paired-end reads using VSEARCH.
        
        Args:
            forward_reads: Forward reads FASTQ file
            reverse_reads: Reverse reads FASTQ file
            output_dir: Output directory
            min_overlap: Minimum overlap between reads
            max_diffs: Maximum number of differences in overlap region
            max_length: Maximum length of merged reads
            threads: Number of threads to use
            **kwargs: Additional VSEARCH parameters
            
        Returns:
            Tuple of (merged_output_path, merge_stats)
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output files
        merged_output = output_dir / f"merged_{forward_reads.stem}.fastq"
        unmerged_f = output_dir / f"unmerged_f_{forward_reads.stem}.fastq"
        unmerged_r = output_dir / f"unmerged_r_{reverse_reads.stem}.fastq"
        
        # Build VSEARCH command
        cmd = [
            'vsearch',
            '--fastq_mergepairs', str(forward_reads),
            '--reverse', str(reverse_reads),
            '--fastqout', str(merged_output),
            '--fastqout_notmerged_fwd', str(unmerged_f),
            '--fastqout_notmerged_rev', str(unmerged_r),
            '--fastq_minovlen', str(min_overlap),
            '--fastq_maxdiffs', str(max_diffs),
            '--fastq_maxlen', str(max_length),
            '--threads', str(threads)
        ]
        
        # Add optional parameters
        if kwargs.get('fastq_minlen'):
            cmd.extend(['--fastq_minlen', str(kwargs['fastq_minlen'])])
        
        if kwargs.get('fastq_maxee'):
            cmd.extend(['--fastq_maxee', str(kwargs['fastq_maxee'])])
        
        if kwargs.get('fastq_truncqual'):
            cmd.extend(['--fastq_truncqual', str(kwargs['fastq_truncqual'])])
        
        # Run VSEARCH
        self.logger.info(f"Running VSEARCH merge command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse merge statistics from stderr
            merge_stats = self._parse_merge_stats(result.stderr)
            
            self.logger.info("VSEARCH merging completed successfully")
            self.logger.info(f"Merge statistics: {merge_stats}")
            
            return merged_output, merge_stats
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"VSEARCH merge failed: {e.stderr}")
            raise
    
    def _parse_merge_stats(self, stderr_output: str) -> Dict:
        """
        Parse merge statistics from VSEARCH stderr output.
        
        Args:
            stderr_output: VSEARCH stderr text
            
        Returns:
            Dictionary with merge statistics
        """
        stats = {}
        
        lines = stderr_output.split('\\n')
        
        for line in lines:
            line = line.strip()
            
            if 'pairs read' in line:
                stats['pairs_read'] = int(line.split()[0])
            elif 'pairs merged' in line and '%' in line:
                parts = line.split()
                stats['pairs_merged'] = int(parts[0])
                # Extract percentage
                pct_str = next(part for part in parts if '%' in part)
                stats['merge_rate'] = float(pct_str.rstrip('%)'))
            elif 'pairs not merged' in line:
                stats['pairs_not_merged'] = int(line.split()[0])
        
        return stats
    
    def quality_filter(self,
                      merged_reads: Path,
                      output_dir: Path,
                      max_ee: float = 1.0,
                      min_length: int = 100,
                      max_length: int = 500,
                      **kwargs) -> Tuple[Path, Dict]:
        """
        Apply quality filtering to merged reads.
        
        Args:
            merged_reads: Merged reads FASTQ file
            output_dir: Output directory
            max_ee: Maximum expected errors
            min_length: Minimum read length
            max_length: Maximum read length
            **kwargs: Additional filtering parameters
            
        Returns:
            Tuple of (filtered_output_path, filter_stats)
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filtered_output = output_dir / f"filtered_{merged_reads.name}"
        
        cmd = [
            'vsearch',
            '--fastq_filter', str(merged_reads),
            '--fastqout', str(filtered_output),
            '--fastq_maxee', str(max_ee),
            '--fastq_minlen', str(min_length),
            '--fastq_maxlen', str(max_length)
        ]
        
        # Add optional parameters
        if kwargs.get('fastq_maxns'):
            cmd.extend(['--fastq_maxns', str(kwargs['fastq_maxns'])])
        
        if kwargs.get('fastq_truncqual'):
            cmd.extend(['--fastq_truncqual', str(kwargs['fastq_truncqual'])])
        
        self.logger.info(f"Running VSEARCH filter command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse filter statistics
            filter_stats = self._parse_filter_stats(result.stderr)
            
            self.logger.info("VSEARCH filtering completed successfully")
            self.logger.info(f"Filter statistics: {filter_stats}")
            
            return filtered_output, filter_stats
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"VSEARCH filter failed: {e.stderr}")
            raise
    
    def _parse_filter_stats(self, stderr_output: str) -> Dict:
        """
        Parse filter statistics from VSEARCH stderr output.
        
        Args:
            stderr_output: VSEARCH stderr text
            
        Returns:
            Dictionary with filter statistics
        """
        stats = {}
        
        lines = stderr_output.split('\\n')
        
        for line in lines:
            line = line.strip()
            
            if 'sequences kept' in line and '%' in line:
                parts = line.split()
                stats['sequences_kept'] = int(parts[0])
                # Extract percentage
                pct_str = next(part for part in parts if '%' in part)
                stats['keep_rate'] = float(pct_str.rstrip('%)'))
            elif 'sequences discarded' in line:
                stats['sequences_discarded'] = int(line.split()[0])
        
        return stats