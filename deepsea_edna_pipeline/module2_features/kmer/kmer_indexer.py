"""
K-mer indexing utilities for fast storage and retrieval.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import logging
import sqlite3
import json


class KmerIndexer:
    """
    Efficient storage and retrieval system for k-mer data.
    """
    
    def __init__(self, db_path: Path, logger: Optional[logging.Logger] = None):
        """
        Initialize k-mer indexer with SQLite backend.
        
        Args:
            db_path: Path to SQLite database
            logger: Logger instance
        """
        self.db_path = db_path
        self.logger = logger or logging.getLogger(__name__)
        
        # Create database and tables
        self._init_database()
    
    def _init_database(self) -> None:
        """Initialize SQLite database schema."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Create tables
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS sequences (
                    id INTEGER PRIMARY KEY,
                    seq_id TEXT UNIQUE NOT NULL,
                    length INTEGER,
                    metadata TEXT
                )
            ''')
            
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS kmers (
                    id INTEGER PRIMARY KEY,
                    kmer TEXT UNIQUE NOT NULL
                )
            ''')
            
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS kmer_counts (
                    seq_id INTEGER,
                    kmer_id INTEGER,
                    count INTEGER,
                    FOREIGN KEY (seq_id) REFERENCES sequences(id),
                    FOREIGN KEY (kmer_id) REFERENCES kmers(id),
                    PRIMARY KEY (seq_id, kmer_id)
                )
            ''')
            
            # Create indices
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_seq_id ON sequences(seq_id)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_kmer ON kmers(kmer)')
            
            conn.commit()
    
    def add_sequence_kmers(self, seq_id: str, kmer_counts: Dict[str, int],
                          metadata: Optional[Dict] = None) -> None:
        """
        Add k-mer counts for a sequence to the database.
        
        Args:
            seq_id: Sequence identifier
            kmer_counts: Dictionary mapping k-mers to counts
            metadata: Optional metadata dictionary
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Insert sequence
            metadata_json = json.dumps(metadata) if metadata else None
            seq_length = sum(kmer_counts.values()) if kmer_counts else 0
            
            cursor.execute('''
                INSERT OR REPLACE INTO sequences (seq_id, length, metadata)
                VALUES (?, ?, ?)
            ''', (seq_id, seq_length, metadata_json))
            
            # Get sequence ID
            cursor.execute('SELECT id FROM sequences WHERE seq_id = ?', (seq_id,))
            seq_db_id = cursor.fetchone()[0]
            
            # Insert k-mers and counts
            for kmer, count in kmer_counts.items():
                # Insert k-mer if not exists
                cursor.execute('INSERT OR IGNORE INTO kmers (kmer) VALUES (?)', (kmer,))
                
                # Get k-mer ID
                cursor.execute('SELECT id FROM kmers WHERE kmer = ?', (kmer,))
                kmer_id = cursor.fetchone()[0]
                
                # Insert count
                cursor.execute('''
                    INSERT OR REPLACE INTO kmer_counts (seq_id, kmer_id, count)
                    VALUES (?, ?, ?)
                ''', (seq_db_id, kmer_id, count))
            
            conn.commit()
    
    def get_sequence_kmers(self, seq_id: str) -> Dict[str, int]:
        """Retrieve k-mer counts for a sequence."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT k.kmer, kc.count
                FROM sequences s
                JOIN kmer_counts kc ON s.id = kc.seq_id
                JOIN kmers k ON kc.kmer_id = k.id
                WHERE s.seq_id = ?
            ''', (seq_id,))
            
            return dict(cursor.fetchall())
    
    def find_sequences_with_kmer(self, kmer: str) -> List[Tuple[str, int]]:
        """Find all sequences containing a specific k-mer."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT s.seq_id, kc.count
                FROM sequences s
                JOIN kmer_counts kc ON s.id = kc.seq_id
                JOIN kmers k ON kc.kmer_id = k.id
                WHERE k.kmer = ?
            ''', (kmer,))
            
            return cursor.fetchall()
    
    def get_database_stats(self) -> Dict:
        """Get statistics about the k-mer database."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Count sequences
            cursor.execute('SELECT COUNT(*) FROM sequences')
            num_sequences = cursor.fetchone()[0]
            
            # Count unique k-mers
            cursor.execute('SELECT COUNT(*) FROM kmers')
            num_kmers = cursor.fetchone()[0]
            
            # Total k-mer instances
            cursor.execute('SELECT SUM(count) FROM kmer_counts')
            total_kmer_counts = cursor.fetchone()[0] or 0
            
            return {
                'num_sequences': num_sequences,
                'num_unique_kmers': num_kmers,
                'total_kmer_counts': total_kmer_counts,
                'avg_kmers_per_sequence': total_kmer_counts / num_sequences if num_sequences > 0 else 0
            }