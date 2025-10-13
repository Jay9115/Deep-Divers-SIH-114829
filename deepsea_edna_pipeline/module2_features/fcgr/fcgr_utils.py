"""
FCGR utilities for normalization and processing.
"""

import numpy as np
from typing import Dict, List, Optional
import logging


class FCGRUtils:
    """Utility functions for FCGR processing."""
    
    @staticmethod
    def apply_gaussian_filter(fcgr: np.ndarray, sigma: float = 1.0) -> np.ndarray:
        """Apply Gaussian smoothing to FCGR."""
        from scipy.ndimage import gaussian_filter
        return gaussian_filter(fcgr, sigma=sigma)
    
    @staticmethod
    def extract_texture_features(fcgr: np.ndarray) -> Dict[str, float]:
        """Extract texture features from FCGR."""
        # Basic texture metrics
        mean_intensity = np.mean(fcgr)
        std_intensity = np.std(fcgr)
        max_intensity = np.max(fcgr)
        
        # Contrast and homogeneity
        contrast = np.var(fcgr)
        
        return {
            'mean_intensity': mean_intensity,
            'std_intensity': std_intensity,
            'max_intensity': max_intensity,
            'contrast': contrast
        }
    
    @staticmethod
    def compress_fcgr(fcgr: np.ndarray, target_size: int) -> np.ndarray:
        """Compress FCGR to smaller size."""
        from skimage.transform import resize
        return resize(fcgr, (target_size, target_size), anti_aliasing=True)