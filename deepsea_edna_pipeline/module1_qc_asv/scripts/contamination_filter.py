"""
Contamination filtering using negative controls.
"""

from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


class ContaminationFilter:
    """
    Filter contaminating sequences identified through negative controls.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize contamination filter.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
    def identify_contaminants(self,
                            asv_table: pd.DataFrame,
                            negative_controls: List[str],
                            prevalence_threshold: float = 0.1,
                            abundance_ratio: float = 0.5) -> Set[str]:
        """
        Identify contaminant ASVs based on negative control analysis.
        
        Args:
            asv_table: ASV abundance table (samples x ASVs)
            negative_controls: List of negative control sample names
            prevalence_threshold: Maximum prevalence in negative controls
            abundance_ratio: Maximum ratio of abundance in negatives vs samples
            
        Returns:
            Set of contaminant ASV identifiers
        """
        self.logger.info("Identifying contaminant sequences from negative controls")
        
        contaminants = set()
        
        # Check if negative controls exist in the data
        available_negatives = [nc for nc in negative_controls if nc in asv_table.index]
        
        if not available_negatives:
            self.logger.warning("No negative controls found in ASV table")
            return contaminants
        
        self.logger.info(f"Found {len(available_negatives)} negative control samples")
        
        # Get sample columns (non-negative controls)
        sample_columns = [col for col in asv_table.columns 
                         if col not in negative_controls]
        
        for asv in asv_table.columns:
            # Calculate prevalence in negative controls
            neg_presence = (asv_table.loc[available_negatives, asv] > 0).sum()
            neg_prevalence = neg_presence / len(available_negatives)
            
            # Calculate abundance ratio
            neg_abundance = asv_table.loc[available_negatives, asv].mean()
            sample_abundance = asv_table.loc[sample_columns, asv].mean() if sample_columns else 0
            
            abundance_ratio_value = (neg_abundance / sample_abundance 
                                   if sample_abundance > 0 else float('inf'))
            
            # Mark as contaminant if exceeds thresholds
            if (neg_prevalence > prevalence_threshold or 
                abundance_ratio_value > abundance_ratio):
                contaminants.add(asv)
                
                self.logger.debug(f"Contaminant identified: {asv} "
                                f"(prevalence: {neg_prevalence:.3f}, "
                                f"abundance_ratio: {abundance_ratio_value:.3f})")
        
        self.logger.info(f"Identified {len(contaminants)} contaminant ASVs")
        
        return contaminants
    
    def filter_sequences(self,
                        input_fasta: Path,
                        contaminant_ids: Set[str],
                        output_dir: Path) -> Tuple[Path, Dict]:
        """
        Remove contaminant sequences from FASTA file.
        
        Args:
            input_fasta: Input FASTA file
            contaminant_ids: Set of contaminant sequence identifiers
            output_dir: Output directory
            
        Returns:
            Tuple of (filtered_fasta_path, filter_stats)
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filtered_output = output_dir / f"decontaminated_{input_fasta.name}"
        contaminant_output = output_dir / f"contaminants_{input_fasta.name}"
        
        clean_sequences = []
        contaminant_sequences = []
        
        # Process sequences
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in contaminant_ids:
                contaminant_sequences.append(record)
            else:
                clean_sequences.append(record)
        
        # Write filtered sequences
        SeqIO.write(clean_sequences, filtered_output, "fasta")
        SeqIO.write(contaminant_sequences, contaminant_output, "fasta")
        
        # Generate statistics
        stats = {
            'input_sequences': len(clean_sequences) + len(contaminant_sequences),
            'clean_sequences': len(clean_sequences),
            'contaminant_sequences': len(contaminant_sequences),
            'contamination_rate': (len(contaminant_sequences) / 
                                 (len(clean_sequences) + len(contaminant_sequences))
                                 if (len(clean_sequences) + len(contaminant_sequences)) > 0 else 0)
        }
        
        self.logger.info(f"Contamination filtering completed. "
                        f"Removed {len(contaminant_sequences)} contaminant sequences")
        
        return filtered_output, stats
    
    def filter_asv_table(self,
                        asv_table: pd.DataFrame,
                        contaminant_ids: Set[str],
                        output_path: Path) -> Tuple[pd.DataFrame, Dict]:
        """
        Remove contaminant ASVs from abundance table.
        
        Args:
            asv_table: ASV abundance table
            contaminant_ids: Set of contaminant ASV identifiers
            output_path: Output file path for filtered table
            
        Returns:
            Tuple of (filtered_table, filter_stats)
        """
        # Filter out contaminant columns
        clean_columns = [col for col in asv_table.columns 
                        if col not in contaminant_ids]
        
        filtered_table = asv_table[clean_columns].copy()
        
        # Save filtered table
        output_path.parent.mkdir(parents=True, exist_ok=True)
        filtered_table.to_csv(output_path, sep='\\t')
        
        # Generate statistics
        stats = {
            'input_asvs': len(asv_table.columns),
            'clean_asvs': len(filtered_table.columns),
            'removed_asvs': len(contaminant_ids),
            'input_reads': asv_table.sum().sum(),
            'clean_reads': filtered_table.sum().sum()
        }
        
        stats['reads_removed_rate'] = ((stats['input_reads'] - stats['clean_reads']) / 
                                     stats['input_reads'] if stats['input_reads'] > 0 else 0)
        
        return filtered_table, stats
    
    def run_complete_decontamination(self,
                                   asv_table_path: Path,
                                   asv_fasta_path: Path,
                                   negative_controls: List[str],
                                   output_dir: Path,
                                   **kwargs) -> Dict:
        """
        Run complete decontamination pipeline.
        
        Args:
            asv_table_path: Path to ASV abundance table
            asv_fasta_path: Path to ASV FASTA file
            negative_controls: List of negative control sample names
            output_dir: Output directory
            **kwargs: Additional filtering parameters
            
        Returns:
            Dictionary with decontamination results
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load ASV table
        asv_table = pd.read_csv(asv_table_path, sep='\\t', index_col=0)
        
        # Identify contaminants
        contaminants = self.identify_contaminants(
            asv_table, negative_controls,
            kwargs.get('prevalence_threshold', 0.1),
            kwargs.get('abundance_ratio', 0.5)
        )
        
        # Filter FASTA file
        clean_fasta, fasta_stats = self.filter_sequences(
            asv_fasta_path, contaminants, output_dir
        )
        
        # Filter ASV table
        clean_table, table_stats = self.filter_asv_table(
            asv_table, contaminants, output_dir / "decontaminated_asv_table.tsv"
        )
        
        # Combined results
        results = {
            'contaminant_asvs': list(contaminants),
            'clean_fasta_path': clean_fasta,
            'clean_table_path': output_dir / "decontaminated_asv_table.tsv",
            'fasta_stats': fasta_stats,
            'table_stats': table_stats
        }
        
        return results