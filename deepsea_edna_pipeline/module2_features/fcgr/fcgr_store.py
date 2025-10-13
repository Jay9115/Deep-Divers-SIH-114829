"""
FCGR storage utilities for saving and loading FCGR data.
"""

import numpy as np
from pathlib import Path
import h5py
import logging


class FCGRStore:
    """Efficient storage system for FCGR data."""
    
    def __init__(self, storage_path: Path, logger: Optional[logging.Logger] = None):
        self.storage_path = storage_path
        self.logger = logger or logging.getLogger(__name__)
    
    def save_fcgr_collection(self, fcgr_dict: Dict[str, np.ndarray], 
                           metadata: Optional[Dict] = None) -> None:
        """Save collection of FCGR arrays to HDF5 file."""
        self.storage_path.parent.mkdir(parents=True, exist_ok=True)
        
        with h5py.File(self.storage_path, 'w') as f:
            # Save FCGRs
            for seq_id, fcgr in fcgr_dict.items():
                f.create_dataset(f"fcgr/{seq_id}", data=fcgr, compression='gzip')
            
            # Save metadata
            if metadata:
                for key, value in metadata.items():
                    f.attrs[key] = value
    
    def load_fcgr_collection(self) -> Dict[str, np.ndarray]:
        """Load FCGR collection from storage."""
        fcgr_dict = {}
        
        with h5py.File(self.storage_path, 'r') as f:
            fcgr_group = f['fcgr']
            for seq_id in fcgr_group.keys():
                fcgr_dict[seq_id] = fcgr_group[seq_id][:]
        
        return fcgr_dict