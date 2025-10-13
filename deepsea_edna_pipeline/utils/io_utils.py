"""
I/O utilities for reading and writing FASTQ files, metadata, and artifacts.
"""

import os
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Union, Optional
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json
import pickle


def read_fastq(file_path: Union[str, Path], format: str = "fastq") -> List[SeqRecord]:
    """
    Read FASTQ file (supports gzipped files).
    
    Args:
        file_path: Path to FASTQ file
        format: File format (default: "fastq")
        
    Returns:
        List of SeqRecord objects
    """
    file_path = Path(file_path)
    
    if file_path.suffix == '.gz':
        with gzip.open(file_path, "rt") as handle:
            return list(SeqIO.parse(handle, format))
    else:
        return list(SeqIO.parse(file_path, format))


def write_fastq(sequences: List[SeqRecord], output_path: Union[str, Path], 
                compress: bool = False) -> None:
    """
    Write sequences to FASTQ file.
    
    Args:
        sequences: List of SeqRecord objects
        output_path: Output file path
        compress: Whether to gzip the output
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if compress:
        with gzip.open(f"{output_path}.gz", "wt") as handle:
            SeqIO.write(sequences, handle, "fastq")
    else:
        SeqIO.write(sequences, output_path, "fastq")


def read_metadata(file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Read metadata file (TSV format).
    
    Args:
        file_path: Path to metadata file
        
    Returns:
        DataFrame with metadata
    """
    return pd.read_csv(file_path, sep='\t')


def write_metadata(metadata: pd.DataFrame, output_path: Union[str, Path]) -> None:
    """
    Write metadata to TSV file.
    
    Args:
        metadata: DataFrame with metadata
        output_path: Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output_path, sep='\t', index=False)


def save_artifacts(artifacts: Dict, output_dir: Union[str, Path], 
                   run_id: str) -> None:
    """
    Save pipeline artifacts (JSON and pickle formats).
    
    Args:
        artifacts: Dictionary of artifacts to save
        output_dir: Output directory
        run_id: Unique run identifier
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save JSON-serializable artifacts
    json_artifacts = {k: v for k, v in artifacts.items() 
                     if isinstance(v, (str, int, float, list, dict, bool))}
    
    with open(output_dir / f"artifacts_{run_id}.json", 'w') as f:
        json.dump(json_artifacts, f, indent=2)
    
    # Save all artifacts as pickle
    with open(output_dir / f"artifacts_{run_id}.pkl", 'wb') as f:
        pickle.dump(artifacts, f)


def load_artifacts(file_path: Union[str, Path]) -> Dict:
    """
    Load pipeline artifacts from pickle file.
    
    Args:
        file_path: Path to artifacts file
        
    Returns:
        Dictionary of loaded artifacts
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)


def create_output_structure(base_dir: Union[str, Path], run_id: str) -> Dict[str, Path]:
    """
    Create standardized output directory structure.
    
    Args:
        base_dir: Base output directory
        run_id: Unique run identifier
        
    Returns:
        Dictionary mapping structure names to paths
    """
    base_dir = Path(base_dir)
    run_dir = base_dir / run_id
    
    structure = {
        'run_dir': run_dir,
        'qc': run_dir / 'qc',
        'asv': run_dir / 'asv', 
        'features': run_dir / 'features',
        'logs': run_dir / 'logs',
        'reports': run_dir / 'reports'
    }
    
    for path in structure.values():
        path.mkdir(parents=True, exist_ok=True)
    
    return structure