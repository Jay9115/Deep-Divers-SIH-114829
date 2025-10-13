"""
Automatic marker detection for 18S and COI genes.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
from Bio import SeqIO
from Bio.Seq import Seq
import re


class MarkerDetector:
    """
    Automatically detect 18S or COI markers in sequence data.
    """
    
    def __init__(self, primers_dir: Optional[Path] = None, 
                 logger: Optional[logging.Logger] = None):
        """
        Initialize marker detector.
        
        Args:
            primers_dir: Directory containing primer JSON files
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        if primers_dir is None:
            primers_dir = Path(__file__).parent.parent / "primers"
        
        self.primers_dir = primers_dir
        self.primer_sets = self._load_primer_sets()
        
    def _load_primer_sets(self) -> Dict:
        """Load primer sets from JSON files."""
        primer_sets = {}
        
        primer_files = {
            '18S': self.primers_dir / "18s_primers.json",
            'COI': self.primers_dir / "coi_primers.json"
        }
        
        for marker, file_path in primer_files.items():
            try:
                if file_path.exists():
                    with open(file_path, 'r') as f:
                        primer_sets[marker] = json.load(f)
                    self.logger.info(f"Loaded {marker} primer set from {file_path}")
                else:
                    self.logger.warning(f"Primer file not found: {file_path}")
            except Exception as e:
                self.logger.error(f"Error loading primer set {marker}: {e}")
        
        return primer_sets
    
    def detect_marker(self, 
                     sequences: List[SeqIO.SeqRecord],
                     min_matches: int = 10,
                     match_threshold: float = 0.8) -> Tuple[str, Dict]:
        """
        Detect the most likely marker gene in sequence data.
        
        Args:
            sequences: List of sequence records to analyze
            min_matches: Minimum number of sequences that must match
            match_threshold: Minimum fraction of sequences that must match
            
        Returns:
            Tuple of (detected_marker, match_statistics)
        """
        self.logger.info(f"Detecting marker gene in {len(sequences)} sequences")
        
        marker_scores = {}
        
        for marker, primer_data in self.primer_sets.items():
            score = self._score_marker_match(sequences, primer_data, marker)
            marker_scores[marker] = score
            
            self.logger.info(f"{marker} match score: {score['match_rate']:.3f} "
                           f"({score['matches']}/{len(sequences)} sequences)")
        
        # Determine best marker
        best_marker = None
        best_score = 0
        
        for marker, score in marker_scores.items():
            if (score['matches'] >= min_matches and 
                score['match_rate'] >= match_threshold and
                score['match_rate'] > best_score):
                best_marker = marker
                best_score = score['match_rate']
        
        if best_marker:
            self.logger.info(f"Detected marker: {best_marker} "
                           f"(confidence: {best_score:.3f})")
        else:
            self.logger.warning("Could not confidently detect marker gene")
            best_marker = "unknown"
        
        return best_marker, marker_scores
    
    def _score_marker_match(self, sequences: List[SeqIO.SeqRecord], 
                          primer_data: Dict, marker: str) -> Dict:
        """
        Score how well sequences match a specific marker's primers.
        """
        if f"{marker.lower()}_primers" not in primer_data:
            return {'matches': 0, 'match_rate': 0.0, 'primer_matches': {}}
        
        primers_info = primer_data[f"{marker.lower()}_primers"]
        primers = primers_info.get('primers', [])
        
        total_matches = 0
        primer_matches = {}
        
        for primer in primers:
            primer_seq = primer['sequence']
            primer_name = primer['name']
            
            # Convert IUPAC codes to regex
            regex_pattern = self._iupac_to_regex(primer_seq)
            
            matches = 0
            for seq_record in sequences:
                seq_str = str(seq_record.seq).upper()
                
                # Search for primer in sequence (both strands)
                if (re.search(regex_pattern, seq_str) or 
                    re.search(regex_pattern, str(Seq(seq_str).reverse_complement()))):
                    matches += 1
            
            primer_matches[primer_name] = {
                'matches': matches,
                'match_rate': matches / len(sequences) if sequences else 0
            }
            
            # Use best primer match for overall score
            total_matches = max(total_matches, matches)
        
        return {
            'matches': total_matches,
            'match_rate': total_matches / len(sequences) if sequences else 0,
            'primer_matches': primer_matches
        }
    
    def _iupac_to_regex(self, sequence: str) -> str:
        """
        Convert IUPAC nucleotide codes to regex pattern.
        """
        iupac_codes = {
            'R': '[AG]',
            'Y': '[CT]', 
            'S': '[GC]',
            'W': '[AT]',
            'K': '[GT]',
            'M': '[AC]',
            'B': '[CGT]',
            'D': '[AGT]',
            'H': '[ACT]',
            'V': '[ACG]',
            'N': '[ACGT]'
        }
        
        regex_pattern = sequence.upper()
        for iupac, regex in iupac_codes.items():
            regex_pattern = regex_pattern.replace(iupac, regex)
        
        return regex_pattern
    
    def suggest_primers(self, detected_marker: str) -> Dict:
        """
        Suggest appropriate primers for detected marker.
        
        Args:
            detected_marker: Detected marker gene
            
        Returns:
            Dictionary with suggested primer information
        """
        if detected_marker not in self.primer_sets:
            return {}
        
        marker_data = self.primer_sets[detected_marker]
        marker_key = f"{detected_marker.lower()}_primers"
        
        if marker_key not in marker_data:
            return {}
        
        primers_info = marker_data[marker_key]
        
        suggestions = {
            'marker': detected_marker,
            'description': primers_info.get('description', ''),
            'recommended_pairs': primers_info.get('primer_pairs', []),
            'processing_parameters': marker_data.get('processing_parameters', {})
        }
        
        return suggestions
    
    def analyze_file(self, input_file: Path, 
                    sample_size: int = 1000) -> Tuple[str, Dict]:
        """
        Analyze a FASTA/FASTQ file to detect marker gene.
        
        Args:
            input_file: Path to sequence file
            sample_size: Number of sequences to sample for analysis
            
        Returns:
            Tuple of (detected_marker, analysis_results)
        """
        self.logger.info(f"Analyzing file: {input_file}")
        
        # Determine file format
        if input_file.suffix.lower() in ['.fq', '.fastq']:
            file_format = 'fastq'
        else:
            file_format = 'fasta'
        
        # Sample sequences for analysis
        sequences = []
        with open(input_file, 'r') as handle:
            for i, record in enumerate(SeqIO.parse(handle, file_format)):
                sequences.append(record)
                if len(sequences) >= sample_size:
                    break
        
        if not sequences:
            raise ValueError(f"No sequences found in {input_file}")
        
        self.logger.info(f"Sampled {len(sequences)} sequences for analysis")
        
        # Detect marker
        detected_marker, match_stats = self.detect_marker(sequences)
        
        # Add file-specific information
        analysis_results = {
            'input_file': str(input_file),
            'sequences_analyzed': len(sequences),
            'file_format': file_format,
            'detected_marker': detected_marker,
            'match_statistics': match_stats
        }
        
        if detected_marker != "unknown":
            analysis_results['suggested_primers'] = self.suggest_primers(detected_marker)
        
        return detected_marker, analysis_results