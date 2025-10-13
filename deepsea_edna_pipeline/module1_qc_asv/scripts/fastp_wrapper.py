"""
FastP wrapper for quality control and trimming.
"""

import subprocess
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging


class FastPWrapper:
    """
    Wrapper for fastp tool to perform quality control and trimming.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize fastp wrapper.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
    def run_qc(self,
               input_r1: Path,
               input_r2: Optional[Path] = None,
               output_dir: Path = None,
               quality_threshold: int = 20,
               min_length: int = 50,
               threads: int = 1,
               **kwargs) -> Tuple[Path, Optional[Path], Dict]:
        """
        Run quality control and trimming with fastp.
        
        Args:
            input_r1: Forward reads FASTQ file
            input_r2: Reverse reads FASTQ file (optional)
            output_dir: Output directory
            quality_threshold: Quality score threshold
            min_length: Minimum read length after trimming
            threads: Number of threads to use
            **kwargs: Additional fastp parameters
            
        Returns:
            Tuple of (output_r1_path, output_r2_path, qc_stats)
        """
        if output_dir is None:
            output_dir = input_r1.parent / 'fastp_output'
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output file paths
        output_r1 = output_dir / f"qc_{input_r1.name}"
        output_r2 = None
        if input_r2:
            output_r2 = output_dir / f"qc_{input_r2.name}"
        
        # Report files
        html_report = output_dir / f"fastp_report_{input_r1.stem}.html"
        json_report = output_dir / f"fastp_report_{input_r1.stem}.json"
        
        # Build fastp command
        cmd = [
            'fastp',
            '--in1', str(input_r1),
            '--out1', str(output_r1),
            '--html', str(html_report),
            '--json', str(json_report),
            '--qualified_quality_phred', str(quality_threshold),
            '--length_required', str(min_length),
            '--thread', str(threads)
        ]
        
        # Add paired-end parameters if R2 provided
        if input_r2 and output_r2:
            cmd.extend(['--in2', str(input_r2), '--out2', str(output_r2)])
        
        # Add additional parameters
        if kwargs.get('cut_front', False):
            cmd.extend(['--cut_front', '--cut_front_window_size', str(kwargs.get('cut_front_window', 4))])
        
        if kwargs.get('cut_tail', False):
            cmd.extend(['--cut_tail', '--cut_tail_window_size', str(kwargs.get('cut_tail_window', 4))])
        
        if kwargs.get('trim_poly_g', False):
            cmd.append('--trim_poly_g')
        
        if kwargs.get('trim_poly_x', False):
            cmd.append('--trim_poly_x')
        
        # Run fastp
        self.logger.info(f"Running fastp command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse JSON report
            with open(json_report, 'r') as f:
                qc_stats = json.load(f)
            
            self.logger.info("FastP quality control completed successfully")
            
            return output_r1, output_r2, qc_stats
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastP failed: {e.stderr}")
            raise
        except json.JSONDecodeError as e:
            self.logger.error(f"Failed to parse FastP JSON report: {e}")
            raise
    
    def extract_qc_metrics(self, qc_stats: Dict) -> Dict:
        """
        Extract key QC metrics from fastp JSON output.
        
        Args:
            qc_stats: FastP JSON statistics
            
        Returns:
            Dictionary with key metrics
        """
        metrics = {}
        
        # Before filtering stats
        before = qc_stats.get('summary', {}).get('before_filtering', {})
        metrics['input_total_reads'] = before.get('total_reads', 0)
        metrics['input_total_bases'] = before.get('total_bases', 0)
        metrics['input_q20_rate'] = before.get('q20_rate', 0)
        metrics['input_q30_rate'] = before.get('q30_rate', 0)
        metrics['input_gc_content'] = before.get('gc_content', 0)
        
        # After filtering stats
        after = qc_stats.get('summary', {}).get('after_filtering', {})
        metrics['output_total_reads'] = after.get('total_reads', 0)
        metrics['output_total_bases'] = after.get('total_bases', 0)
        metrics['output_q20_rate'] = after.get('q20_rate', 0)
        metrics['output_q30_rate'] = after.get('q30_rate', 0)
        metrics['output_gc_content'] = after.get('gc_content', 0)
        
        # Calculate filtering stats
        if metrics['input_total_reads'] > 0:
            metrics['reads_passed_rate'] = metrics['output_total_reads'] / metrics['input_total_reads']
            metrics['bases_passed_rate'] = metrics['output_total_bases'] / metrics['input_total_bases']
        else:
            metrics['reads_passed_rate'] = 0
            metrics['bases_passed_rate'] = 0
        
        # Filtering result
        filtering = qc_stats.get('filtering_result', {})
        metrics['reads_with_low_quality'] = filtering.get('low_quality_reads', 0)
        metrics['reads_with_too_many_N'] = filtering.get('too_many_N_reads', 0)
        metrics['reads_too_short'] = filtering.get('too_short_reads', 0)
        
        return metrics
    
    def batch_process(self,
                     input_files: List[Tuple[Path, Optional[Path]]],
                     output_dir: Path,
                     **kwargs) -> List[Tuple[Path, Optional[Path], Dict]]:
        """
        Process multiple file pairs in batch.
        
        Args:
            input_files: List of (R1, R2) file path tuples
            output_dir: Output directory 
            **kwargs: FastP parameters
            
        Returns:
            List of (output_R1, output_R2, stats) tuples
        """
        results = []
        
        for i, (r1, r2) in enumerate(input_files):
            self.logger.info(f"Processing file pair {i+1}/{len(input_files)}: {r1.name}")
            
            sample_output_dir = output_dir / r1.stem
            output_r1, output_r2, stats = self.run_qc(
                r1, r2, sample_output_dir, **kwargs
            )
            
            results.append((output_r1, output_r2, stats))
        
        return results