"""
Main ASV pipeline orchestrator for Module 1.
"""

import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging
import yaml
import subprocess
import tempfile

from ..utils.logging_utils import create_run_logger, log_step_start, log_step_complete
from ..utils.io_utils import create_output_structure, save_artifacts
from ..utils.system_check import check_dependencies
from .scripts.marker_detect import MarkerDetector
from .scripts.cutadapt_wrapper import CutadaptWrapper
from .scripts.fastp_wrapper import FastPWrapper
from .scripts.vsearch_merge import VSearchMerger
from .scripts.chimera_filter import ChimeraFilter
from .scripts.contamination_filter import ContaminationFilter


class ASVPipeline:
    """
    Complete ASV generation pipeline from raw reads to clean ASVs.
    """
    
    def __init__(self, config_file: Optional[Path] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize ASV pipeline.
        
        Args:
            config_file: Path to configuration YAML file
            logger: Logger instance
        """
        # Load configuration
        if config_file is None:
            config_file = Path(__file__).parent / "config.yaml"
        
        self.config = self._load_config(config_file)
        
        # Set up logging
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize processing modules
        self.marker_detector = MarkerDetector(logger=self.logger)
        self.cutadapt_wrapper = CutadaptWrapper(logger=self.logger)
        self.fastp_wrapper = FastPWrapper(logger=self.logger)
        self.vsearch_merger = VSearchMerger(logger=self.logger)
        self.chimera_filter = ChimeraFilter(logger=self.logger)
        self.contamination_filter = ContaminationFilter(logger=self.logger)
        
    def _load_config(self, config_file: Path) -> Dict:
        """Load pipeline configuration from YAML file."""
        try:
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            # Return default configuration if file not found
            return self._get_default_config()
    
    def _get_default_config(self) -> Dict:
        """Get default pipeline configuration."""
        return {
            'quality_control': {
                'min_quality': 20,
                'min_length': 50,
                'max_expected_errors': 1.0
            },
            'primer_trimming': {
                'max_error_rate': 0.1,
                'min_overlap': 3
            },
            'read_merging': {
                'min_overlap': 10,
                'max_diffs': 5,
                'max_length': 600
            },
            'chimera_removal': {
                'abundance_skew': 2.0
            },
            'contamination_filter': {
                'prevalence_threshold': 0.1,
                'abundance_ratio': 0.5
            },
            'dada2': {
                'pool': False,
                'error_estimation_samples': 1000000
            }
        }
    
    def run_complete_pipeline(self,
                            input_dir: Path,
                            output_dir: Path,
                            metadata_file: Optional[Path] = None,
                            reference_db: Optional[Path] = None,
                            negative_controls: Optional[List[str]] = None,
                            run_id: Optional[str] = None) -> Dict:
        """
        Run the complete ASV pipeline.
        
        Args:
            input_dir: Directory containing input FASTQ files
            output_dir: Output directory for results
            metadata_file: Optional metadata file
            reference_db: Optional reference database for chimera detection
            negative_controls: List of negative control sample names
            run_id: Optional run identifier
            
        Returns:
            Dictionary with pipeline results and statistics
        """
        # Generate run ID if not provided
        if run_id is None:
            run_id = f"asv_pipeline_{int(time.time())}"
        
        # Set up logging and output structure
        output_structure = create_output_structure(output_dir, run_id)
        run_logger = create_run_logger(run_id, output_structure['logs'])
        
        run_logger.info(f"Starting ASV pipeline run: {run_id}")
        
        # Check system dependencies
        deps_ok, dep_status = check_dependencies(run_logger)
        if not deps_ok:
            raise RuntimeError("Missing critical dependencies. Check logs for details.")
        
        pipeline_start_time = time.time()
        
        try:
            # Step 1: Auto-detect marker gene
            run_logger.info("=== Step 1: Marker Detection ===")
            detected_marker = self._run_marker_detection(
                input_dir, output_structure['qc'], run_logger
            )
            
            # Step 2: Quality control and trimming
            run_logger.info("=== Step 2: Quality Control ===")
            qc_files = self._run_quality_control(
                input_dir, output_structure['qc'], run_logger
            )
            
            # Step 3: Primer trimming (if marker detected)
            run_logger.info("=== Step 3: Primer Trimming ===")
            trimmed_files = self._run_primer_trimming(
                qc_files, detected_marker, output_structure['qc'], run_logger
            )
            
            # Step 4: Read merging (for paired-end data)
            run_logger.info("=== Step 4: Read Merging ===")
            merged_files = self._run_read_merging(
                trimmed_files, output_structure['qc'], run_logger
            )
            
            # Step 5: DADA2 denoising
            run_logger.info("=== Step 5: DADA2 Denoising ===")
            asv_results = self._run_dada2_denoising(
                merged_files, output_structure['asv'], run_logger
            )
            
            # Step 6: Chimera removal
            run_logger.info("=== Step 6: Chimera Removal ===")
            nonchimeric_asvs = self._run_chimera_removal(
                asv_results, reference_db, output_structure['asv'], run_logger
            )
            
            # Step 7: Contamination filtering
            run_logger.info("=== Step 7: Contamination Filtering ===")
            final_results = self._run_contamination_filtering(
                nonchimeric_asvs, negative_controls, output_structure['asv'], run_logger
            )
            
            # Calculate total execution time
            total_time = time.time() - pipeline_start_time
            
            # Compile final results
            pipeline_results = {
                'run_id': run_id,
                'execution_time': total_time,
                'detected_marker': detected_marker,
                'output_structure': {k: str(v) for k, v in output_structure.items()},
                'final_asv_fasta': str(final_results['clean_fasta_path']),
                'final_asv_table': str(final_results['clean_table_path']),
                'pipeline_stats': final_results
            }
            
            # Save pipeline results
            save_artifacts(pipeline_results, output_structure['run_dir'], run_id)
            
            run_logger.info(f"Pipeline completed successfully in {total_time:.2f} seconds")
            
            return pipeline_results
            
        except Exception as e:
            run_logger.error(f"Pipeline failed: {str(e)}")
            raise
    
    def _run_marker_detection(self, input_dir: Path, output_dir: Path, 
                            logger: logging.Logger) -> str:
        """Run marker gene detection step."""
        log_step_start(logger, "marker_detection")
        
        # Find first FASTQ file for analysis
        fastq_files = list(input_dir.glob("*.fastq")) + list(input_dir.glob("*.fq"))
        
        if not fastq_files:
            logger.warning("No FASTQ files found for marker detection")
            return "unknown"
        
        # Analyze first file
        detected_marker, analysis = self.marker_detector.analyze_file(fastq_files[0])
        
        # Save analysis results
        with open(output_dir / "marker_detection_results.yaml", 'w') as f:
            yaml.dump(analysis, f)
        
        log_step_complete(logger, "marker_detection", detected_marker=detected_marker)
        
        return detected_marker
    
    def _run_quality_control(self, input_dir: Path, output_dir: Path,
                           logger: logging.Logger) -> List[Tuple[Path, Optional[Path]]]:
        """Run quality control step."""
        log_step_start(logger, "quality_control")
        
        # Find input file pairs
        file_pairs = self._find_file_pairs(input_dir)
        
        qc_config = self.config['quality_control']
        
        # Run FastP on all file pairs
        qc_results = self.fastp_wrapper.batch_process(
            file_pairs,
            output_dir / "fastp_results",
            quality_threshold=qc_config['min_quality'],
            min_length=qc_config['min_length']
        )
        
        # Extract file paths from results
        qc_files = [(result[0], result[1]) for result in qc_results]
        
        log_step_complete(logger, "quality_control", 
                         processed_samples=len(qc_files))
        
        return qc_files
    
    def _run_primer_trimming(self, qc_files: List[Tuple[Path, Optional[Path]]],
                           detected_marker: str, output_dir: Path,
                           logger: logging.Logger) -> List[Tuple[Path, Optional[Path]]]:
        """Run primer trimming step."""
        log_step_start(logger, "primer_trimming")
        
        if detected_marker == "unknown":
            logger.warning("Skipping primer trimming - marker not detected")
            return qc_files
        
        # Get primer suggestions for detected marker
        primer_suggestions = self.marker_detector.suggest_primers(detected_marker)
        
        if not primer_suggestions or not primer_suggestions.get('recommended_pairs'):
            logger.warning("No primer information available - skipping trimming")
            return qc_files
        
        # Use first recommended primer pair
        primer_pair = primer_suggestions['recommended_pairs'][0]
        
        # Get primer sequences (simplified - would need proper primer lookup)
        forward_primer = "CCAGCASCYGCGGTAATTCC"  # Example 18S primer
        reverse_primer = "ACTTTCGTTCTTGATYRA"
        
        trimmed_files = []
        trim_config = self.config['primer_trimming']
        
        for r1_file, r2_file in qc_files:
            trimmed_r1, trimmed_r2, stats = self.cutadapt_wrapper.trim_primers(
                r1_file, r2_file, forward_primer, reverse_primer,
                output_dir / "trimmed",
                max_error_rate=trim_config['max_error_rate'],
                min_overlap=trim_config['min_overlap']
            )
            
            trimmed_files.append((trimmed_r1, trimmed_r2))
        
        log_step_complete(logger, "primer_trimming", 
                         processed_samples=len(trimmed_files))
        
        return trimmed_files
    
    def _run_read_merging(self, trimmed_files: List[Tuple[Path, Optional[Path]]],
                        output_dir: Path, logger: logging.Logger) -> List[Path]:
        """Run read merging step."""
        log_step_start(logger, "read_merging")
        
        merged_files = []
        merge_config = self.config['read_merging']
        
        for r1_file, r2_file in trimmed_files:
            if r2_file is None:
                # Single-end data - no merging needed
                merged_files.append(r1_file)
                continue
            
            merged_output, merge_stats = self.vsearch_merger.merge_pairs(
                r1_file, r2_file,
                output_dir / "merged",
                min_overlap=merge_config['min_overlap'],
                max_diffs=merge_config['max_diffs'],
                max_length=merge_config['max_length']
            )
            
            merged_files.append(merged_output)
        
        log_step_complete(logger, "read_merging", 
                         processed_samples=len(merged_files))
        
        return merged_files
    
    def _run_dada2_denoising(self, merged_files: List[Path], output_dir: Path,
                           logger: logging.Logger) -> Dict:
        """Run DADA2 denoising step."""
        log_step_start(logger, "dada2_denoising")
        
        # Prepare sample names file
        sample_names = [f.stem.replace("merged_", "").replace("qc_", "") 
                       for f in merged_files]
        
        sample_names_file = output_dir / "sample_names.txt"
        with open(sample_names_file, 'w') as f:
            for name in sample_names:
                f.write(f"{name}\\n")
        
        # Create input directory with standardized names
        dada2_input_dir = output_dir / "dada2_input"
        dada2_input_dir.mkdir(exist_ok=True)
        
        for i, merged_file in enumerate(merged_files):
            standardized_name = f"{sample_names[i]}.fastq"
            target_path = dada2_input_dir / standardized_name
            
            # Copy or symlink file
            if not target_path.exists():
                target_path.symlink_to(merged_file.absolute())
        
        # Run DADA2 R script
        dada2_output_dir = output_dir / "dada2_output"
        dada2_script = Path(__file__).parent / "scripts" / "dada2_runner.R"
        
        cmd = [
            'Rscript', str(dada2_script),
            str(dada2_input_dir),
            str(dada2_output_dir),
            str(sample_names_file)
        ]
        
        logger.info(f"Running DADA2 command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, 
                                  check=True, timeout=3600)
            logger.info("DADA2 completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"DADA2 failed: {e.stderr}")
            raise
        
        asv_results = {
            'asv_fasta': dada2_output_dir / "asvs.fasta",
            'asv_table': dada2_output_dir / "asv_table.tsv",
            'processing_summary': dada2_output_dir / "processing_summary.tsv"
        }
        
        log_step_complete(logger, "dada2_denoising")
        
        return asv_results
    
    def _run_chimera_removal(self, asv_results: Dict, reference_db: Optional[Path],
                           output_dir: Path, logger: logging.Logger) -> Dict:
        """Run chimera removal step."""
        log_step_start(logger, "chimera_removal")
        
        chimera_config = self.config['chimera_removal']
        
        nonchimeric_fasta, chimera_stats = self.chimera_filter.remove_chimeras(
            asv_results['asv_fasta'],
            output_dir / "chimera_filtered",
            reference_db=reference_db,
            abundance_skew=chimera_config['abundance_skew']
        )
        
        asv_results['nonchimeric_fasta'] = nonchimeric_fasta
        asv_results['chimera_stats'] = chimera_stats
        
        log_step_complete(logger, "chimera_removal", **chimera_stats['overall'])
        
        return asv_results
    
    def _run_contamination_filtering(self, asv_results: Dict, 
                                   negative_controls: Optional[List[str]],
                                   output_dir: Path, logger: logging.Logger) -> Dict:
        """Run contamination filtering step."""
        log_step_start(logger, "contamination_filtering")
        
        if not negative_controls:
            logger.warning("No negative controls provided - skipping contamination filtering")
            return {
                'clean_fasta_path': asv_results['nonchimeric_fasta'],
                'clean_table_path': asv_results['asv_table']
            }
        
        contam_config = self.config['contamination_filter']
        
        decontam_results = self.contamination_filter.run_complete_decontamination(
            asv_results['asv_table'],
            asv_results['nonchimeric_fasta'],
            negative_controls,
            output_dir / "decontaminated",
            **contam_config
        )
        
        log_step_complete(logger, "contamination_filtering", 
                         **decontam_results['table_stats'])
        
        return decontam_results
    
    def _find_file_pairs(self, input_dir: Path) -> List[Tuple[Path, Optional[Path]]]:
        """Find paired-end file pairs in input directory."""
        # Simple pattern matching for R1/R2 files
        r1_files = list(input_dir.glob("*R1*.fastq")) + list(input_dir.glob("*_1.fastq"))
        
        file_pairs = []
        
        for r1_file in r1_files:
            # Try to find corresponding R2 file
            r2_pattern = r1_file.name.replace("R1", "R2").replace("_1", "_2")
            r2_file = input_dir / r2_pattern
            
            if r2_file.exists():
                file_pairs.append((r1_file, r2_file))
            else:
                # Single-end file
                file_pairs.append((r1_file, None))
        
        # Add any remaining single files
        all_files = list(input_dir.glob("*.fastq")) + list(input_dir.glob("*.fq"))
        processed_files = {pair[0] for pair in file_pairs}
        processed_files.update({pair[1] for pair in file_pairs if pair[1]})
        
        for fastq_file in all_files:
            if fastq_file not in processed_files:
                file_pairs.append((fastq_file, None))
        
        return file_pairs