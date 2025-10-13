"""
Quality control metrics calculation and reporting utilities.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_basic_stats(sequences: List[SeqRecord]) -> Dict[str, float]:
    """
    Calculate basic sequence statistics.
    
    Args:
        sequences: List of SeqRecord objects
        
    Returns:
        Dictionary with basic statistics
    """
    if not sequences:
        return {}
    
    lengths = [len(seq) for seq in sequences]
    
    stats = {
        'total_sequences': len(sequences),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'std_length': np.std(lengths)
    }
    
    return stats


def calculate_quality_stats(sequences: List[SeqRecord]) -> Dict[str, float]:
    """
    Calculate quality score statistics for sequences with quality information.
    
    Args:
        sequences: List of SeqRecord objects with quality scores
        
    Returns:
        Dictionary with quality statistics
    """
    if not sequences or not hasattr(sequences[0], 'letter_annotations'):
        return {}
    
    all_qualities = []
    for seq in sequences:
        if 'phred_quality' in seq.letter_annotations:
            all_qualities.extend(seq.letter_annotations['phred_quality'])
    
    if not all_qualities:
        return {}
    
    stats = {
        'mean_quality': np.mean(all_qualities),
        'median_quality': np.median(all_qualities),
        'min_quality': min(all_qualities),
        'max_quality': max(all_qualities),
        'q25_quality': np.percentile(all_qualities, 25),
        'q75_quality': np.percentile(all_qualities, 75)
    }
    
    return stats


def calculate_gc_content(sequences: List[SeqRecord]) -> Dict[str, float]:
    """
    Calculate GC content statistics.
    
    Args:
        sequences: List of SeqRecord objects
        
    Returns:
        Dictionary with GC content statistics
    """
    if not sequences:
        return {}
    
    gc_contents = []
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        total_bases = len(seq_str)
        if total_bases > 0:
            gc_content = (gc_count / total_bases) * 100
            gc_contents.append(gc_content)
    
    if not gc_contents:
        return {}
    
    stats = {
        'mean_gc_content': np.mean(gc_contents),
        'median_gc_content': np.median(gc_contents),
        'std_gc_content': np.std(gc_contents)
    }
    
    return stats


def generate_length_distribution_plot(sequences: List[SeqRecord], 
                                    output_path: Optional[Path] = None) -> plt.Figure:
    """
    Generate sequence length distribution plot.
    
    Args:
        sequences: List of SeqRecord objects
        output_path: Optional path to save plot
        
    Returns:
        Matplotlib figure object
    """
    lengths = [len(seq) for seq in sequences]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Frequency')
    ax.set_title('Sequence Length Distribution')
    ax.grid(True, alpha=0.3)
    
    # Add summary statistics as text
    stats_text = f'Mean: {np.mean(lengths):.1f} bp\\n'
    stats_text += f'Median: {np.median(lengths):.1f} bp\\n'
    stats_text += f'Total sequences: {len(lengths):,}'
    ax.text(0.75, 0.75, stats_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig


def generate_qc_report(qc_results: Dict, output_path: Path) -> None:
    """
    Generate comprehensive QC report.
    
    Args:
        qc_results: Dictionary containing QC metrics
        output_path: Path to save the report
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write("# Quality Control Report\\n\\n")
        f.write(f"Generated on: {pd.Timestamp.now()}\\n\\n")
        
        for stage, metrics in qc_results.items():
            f.write(f"## {stage.replace('_', ' ').title()}\\n\\n")
            
            if isinstance(metrics, dict):
                for metric, value in metrics.items():
                    if isinstance(value, float):
                        f.write(f"- **{metric.replace('_', ' ').title()}**: {value:.2f}\\n")
                    else:
                        f.write(f"- **{metric.replace('_', ' ').title()}**: {value}\\n")
            f.write("\\n")


def calculate_qc_metrics(sequences: List[SeqRecord], 
                        stage_name: str = "raw") -> Dict[str, Dict]:
    """
    Calculate comprehensive QC metrics for a set of sequences.
    
    Args:
        sequences: List of SeqRecord objects
        stage_name: Name of the processing stage
        
    Returns:
        Dictionary with comprehensive QC metrics
    """
    metrics = {}
    
    metrics['basic_stats'] = calculate_basic_stats(sequences)
    metrics['quality_stats'] = calculate_quality_stats(sequences)
    metrics['gc_content'] = calculate_gc_content(sequences)
    
    return {stage_name: metrics}