"""
Logging utilities with provenance tracking for reproducible analyses.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime
import json
import colorlog


def setup_logger(name: str, log_file: Optional[Path] = None, 
                level: int = logging.INFO, 
                console_output: bool = True) -> logging.Logger:
    """
    Set up a logger with file and console handlers.
    
    Args:
        name: Logger name
        log_file: Optional log file path
        level: Logging level
        console_output: Whether to output to console
        
    Returns:
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # File handler
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    # Console handler with colors
    if console_output:
        console_handler = colorlog.StreamHandler()
        console_handler.setLevel(level)
        
        color_formatter = colorlog.ColoredFormatter(
            '%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            log_colors={
                'DEBUG': 'cyan',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'red,bg_white',
            }
        )
        
        console_handler.setFormatter(color_formatter)
        logger.addHandler(console_handler)
    
    return logger


def log_provenance(logger: logging.Logger, stage: str, 
                  parameters: Dict[str, Any], 
                  input_files: Dict[str, str],
                  output_files: Dict[str, str],
                  execution_time: Optional[float] = None) -> None:
    """
    Log provenance information for reproducibility.
    
    Args:
        logger: Logger instance
        stage: Processing stage name
        parameters: Stage parameters
        input_files: Input file paths and their hashes
        output_files: Output file paths and their hashes
        execution_time: Optional execution time in seconds
    """
    provenance = {
        'timestamp': datetime.now().isoformat(),
        'stage': stage,
        'parameters': parameters,
        'input_files': input_files,
        'output_files': output_files,
        'execution_time_seconds': execution_time
    }
    
    logger.info(f"PROVENANCE: {json.dumps(provenance, indent=2)}")


def log_system_info(logger: logging.Logger) -> None:
    """
    Log system information for debugging and reproducibility.
    
    Args:
        logger: Logger instance
    """
    import platform
    import sys
    
    system_info = {
        'python_version': sys.version,
        'platform': platform.platform(),
        'processor': platform.processor(),
        'python_executable': sys.executable
    }
    
    logger.info(f"SYSTEM_INFO: {json.dumps(system_info, indent=2)}")


def log_step_start(logger: logging.Logger, step_name: str, **kwargs) -> None:
    """
    Log the start of a processing step.
    
    Args:
        logger: Logger instance
        step_name: Name of the processing step
        **kwargs: Additional parameters to log
    """
    logger.info(f"Starting step: {step_name}")
    if kwargs:
        logger.info(f"Parameters: {json.dumps(kwargs, indent=2)}")


def log_step_complete(logger: logging.Logger, step_name: str, 
                     execution_time: Optional[float] = None,
                     **results) -> None:
    """
    Log the completion of a processing step.
    
    Args:
        logger: Logger instance
        step_name: Name of the processing step
        execution_time: Execution time in seconds
        **results: Step results to log
    """
    message = f"Completed step: {step_name}"
    if execution_time is not None:
        message += f" (execution time: {execution_time:.2f}s)"
    
    logger.info(message)
    
    if results:
        logger.info(f"Results: {json.dumps(results, indent=2)}")


def create_run_logger(run_id: str, log_dir: Path, 
                     level: int = logging.INFO) -> logging.Logger:
    """
    Create a logger for a specific pipeline run.
    
    Args:
        run_id: Unique run identifier
        log_dir: Directory for log files
        level: Logging level
        
    Returns:
        Configured logger for the run
    """
    log_file = log_dir / f"{run_id}.log"
    logger = setup_logger(f"deepsea_edna.run.{run_id}", log_file, level)
    
    # Log run initialization
    logger.info(f"Initialized logging for run: {run_id}")
    log_system_info(logger)
    
    return logger