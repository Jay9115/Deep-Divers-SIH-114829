"""
Parallel processing utilities for efficient pipeline execution.
"""

import multiprocessing as mp
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, List, Any, Optional, Union
import logging
from functools import partial
import time


class ParallelProcessor:
    """
    Manager for parallel processing operations in the pipeline.
    """
    
    def __init__(self, n_processes: Optional[int] = None, 
                 logger: Optional[logging.Logger] = None):
        """
        Initialize parallel processor.
        
        Args:
            n_processes: Number of processes to use (default: CPU count)
            logger: Logger instance for progress tracking
        """
        self.n_processes = n_processes or mp.cpu_count()
        self.logger = logger or logging.getLogger(__name__)
        
    def map_function(self, func: Callable, items: List[Any], 
                    chunksize: Optional[int] = None) -> List[Any]:
        """
        Apply function to items in parallel using multiprocessing.
        
        Args:
            func: Function to apply to each item
            items: List of items to process
            chunksize: Size of chunks for processing
            
        Returns:
            List of results
        """
        if len(items) == 0:
            return []
        
        if len(items) == 1 or self.n_processes == 1:
            return [func(item) for item in items]
        
        self.logger.info(f"Processing {len(items)} items with {self.n_processes} processes")
        
        start_time = time.time()
        
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            results = list(executor.map(func, items, chunksize=chunksize))
        
        execution_time = time.time() - start_time
        self.logger.info(f"Parallel processing completed in {execution_time:.2f} seconds")
        
        return results
    
    def map_function_with_progress(self, func: Callable, items: List[Any],
                                  progress_callback: Optional[Callable] = None) -> List[Any]:
        """
        Apply function to items in parallel with progress tracking.
        
        Args:
            func: Function to apply to each item
            items: List of items to process
            progress_callback: Optional callback for progress updates
            
        Returns:
            List of results
        """
        if len(items) == 0:
            return []
        
        results = []
        
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            # Submit all tasks
            future_to_item = {executor.submit(func, item): item for item in items}
            
            completed = 0
            for future in future_to_item:
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1
                    
                    if progress_callback:
                        progress_callback(completed, len(items))
                    
                    if completed % max(1, len(items) // 10) == 0:
                        self.logger.info(f"Progress: {completed}/{len(items)} ({100*completed/len(items):.1f}%)")
                        
                except Exception as e:
                    self.logger.error(f"Error processing item: {e}")
                    results.append(None)
        
        return results
    
    def starmap_function(self, func: Callable, arg_tuples: List[tuple]) -> List[Any]:
        """
        Apply function to argument tuples in parallel.
        
        Args:
            func: Function to apply
            arg_tuples: List of argument tuples
            
        Returns:
            List of results
        """
        if len(arg_tuples) == 0:
            return []
        
        if len(arg_tuples) == 1 or self.n_processes == 1:
            return [func(*args) for args in arg_tuples]
        
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            results = list(executor.map(lambda args: func(*args), arg_tuples))
        
        return results


def chunk_list(lst: List[Any], chunk_size: int) -> List[List[Any]]:
    """
    Split a list into chunks of specified size.
    
    Args:
        lst: List to split
        chunk_size: Size of each chunk
        
    Returns:
        List of chunks
    """
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def optimal_chunk_size(total_items: int, n_processes: int) -> int:
    """
    Calculate optimal chunk size for parallel processing.
    
    Args:
        total_items: Total number of items to process
        n_processes: Number of processes
        
    Returns:
        Optimal chunk size
    """
    if total_items <= n_processes:
        return 1
    
    # Aim for 2-4 chunks per process
    chunks_per_process = 3
    chunk_size = max(1, total_items // (n_processes * chunks_per_process))
    
    return chunk_size


def parallel_file_processing(file_paths: List[str], 
                           process_func: Callable[[str], Any],
                           n_processes: Optional[int] = None,
                           logger: Optional[logging.Logger] = None) -> List[Any]:
    """
    Process multiple files in parallel.
    
    Args:
        file_paths: List of file paths to process
        process_func: Function to process each file
        n_processes: Number of processes to use
        logger: Logger instance
        
    Returns:
        List of processing results
    """
    processor = ParallelProcessor(n_processes, logger)
    return processor.map_function(process_func, file_paths)


def parallel_sequence_processing(sequences: List[Any],
                               process_func: Callable[[Any], Any],
                               batch_size: int = 1000,
                               n_processes: Optional[int] = None,
                               logger: Optional[logging.Logger] = None) -> List[Any]:
    """
    Process sequences in parallel batches.
    
    Args:
        sequences: List of sequences to process
        process_func: Function to process each sequence
        batch_size: Number of sequences per batch
        n_processes: Number of processes to use
        logger: Logger instance
        
    Returns:
        Flattened list of processing results
    """
    if not sequences:
        return []
    
    # Split sequences into batches
    batches = chunk_list(sequences, batch_size)
    
    # Create batch processing function
    def process_batch(batch):
        return [process_func(seq) for seq in batch]
    
    processor = ParallelProcessor(n_processes, logger)
    batch_results = processor.map_function(process_batch, batches)
    
    # Flatten results
    results = []
    for batch_result in batch_results:
        results.extend(batch_result)
    
    return results