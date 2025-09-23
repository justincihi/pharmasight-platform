"""
Comprehensive Compound Discovery Logging System for PharmaSight™
Provides regulatory-compliant logging for all discovered compounds and research references.
"""

import sqlite3
import json
import datetime
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
import os

@dataclass
class DiscoveryLog:
    """Compound discovery log data structure."""
    id: int
    compound_name: str
    smiles: str
    category: str
    confidence_score: float
    ip_status: str
    discovery_date: str
    source_research_topic_id: int
    generated_by: str
    molecular_weight: float
    logp: float
    receptor_activity: Dict[str, float]
    binding_affinity: float
    notes: str

@dataclass
class ResearchReference:
    """Research reference (DOI) log data structure."""
    id: int
    doi: str
    title: str
    authors: List[str]
    journal: str
    publication_year: int
    discovery_log_id: int
    retrieval_date: str

class DiscoveryLogger:
    """Manages compound discovery and research reference logging with SQLite."""
    
    def __init__(self, db_path: str = "/home/ubuntu/pharmasight-platform/data/discovery_logs.db"):
        self.db_path = db_path
        self.ensure_data_directory()
        self.init_database()
    
    def ensure_data_directory(self):
        """Ensure the data directory exists."""
        data_dir = os.path.dirname(self.db_path)
        os.makedirs(data_dir, exist_ok=True)
    
    def init_database(self):
        """Initialize the SQLite database with logging tables."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Create discovery logs table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS discovery_logs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    compound_name TEXT NOT NULL,
                    smiles TEXT NOT NULL UNIQUE,
                    category TEXT,
                    confidence_score REAL NOT NULL,
                    ip_status TEXT,
                    discovery_date TEXT NOT NULL,
                    source_research_topic_id INTEGER,
                    generated_by TEXT,
                    molecular_weight REAL,
                    logp REAL,
                    receptor_activity TEXT,  -- JSON object
                    binding_affinity REAL,
                    notes TEXT,
                    FOREIGN KEY (source_research_topic_id) REFERENCES research_topics (id)
                )
            """)
            
            # Create research references table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS research_references (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    doi TEXT NOT NULL UNIQUE,
                    title TEXT,
                    authors TEXT,  -- JSON array
                    journal TEXT,
                    publication_year INTEGER,
                    discovery_log_id INTEGER,
                    retrieval_date TEXT NOT NULL,
                    FOREIGN KEY (discovery_log_id) REFERENCES discovery_logs (id)
                )
            """)
            
            conn.commit()
    
    def log_discovery(self, compound_name: str, smiles: str, category: str, 
                      confidence_score: float, ip_status: str, 
                      source_research_topic_id: int, generated_by: str = "AI Engine",
                      molecular_weight: float = None, logp: float = None,
                      receptor_activity: Dict[str, float] = None, 
                      binding_affinity: float = None, notes: str = None) -> int:
        """Log a new compound discovery."""
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            try:
                cursor.execute("""
                    INSERT INTO discovery_logs 
                    (compound_name, smiles, category, confidence_score, ip_status, 
                     discovery_date, source_research_topic_id, generated_by, 
                     molecular_weight, logp, receptor_activity, binding_affinity, notes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    compound_name, smiles, category, confidence_score, ip_status,
                    timestamp, source_research_topic_id, generated_by,
                    molecular_weight, logp, json.dumps(receptor_activity) if receptor_activity else None,
                    binding_affinity, notes
                ))
                
                log_id = cursor.lastrowid
                conn.commit()
                return log_id
            
            except sqlite3.IntegrityError:
                # SMILES string already exists, update the log
                cursor.execute("SELECT id FROM discovery_logs WHERE smiles = ?", (smiles,))
                log_id = cursor.fetchone()[0]
                
                # For now, we just return the existing ID. In a real system,
                # you might want to update the existing record if new data is better.
                return log_id
    
    def log_research_reference(self, doi: str, discovery_log_id: int, 
                               title: str = None, authors: List[str] = None, 
                               journal: str = None, publication_year: int = None) -> int:
        """Log a research reference (DOI) associated with a discovery."""
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            try:
                cursor.execute("""
                    INSERT INTO research_references 
                    (doi, title, authors, journal, publication_year, discovery_log_id, retrieval_date)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (
                    doi, title, json.dumps(authors) if authors else None, 
                    journal, publication_year, discovery_log_id, timestamp
                ))
                
                ref_id = cursor.lastrowid
                conn.commit()
                return ref_id
            
            except sqlite3.IntegrityError:
                cursor.execute("SELECT id FROM research_references WHERE doi = ?", (doi,))
                return cursor.fetchone()[0]
    
    def get_discoveries(self, confidence_threshold: float = None, ip_status: str = None, 
                          limit: int = 100, offset: int = 0) -> List[DiscoveryLog]:
        """Get discovery logs with optional filtering."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            query = "SELECT * FROM discovery_logs WHERE 1=1"
            params = []
            
            if confidence_threshold is not None:
                query += " AND confidence_score >= ?"
                params.append(confidence_threshold)
            
            if ip_status:
                query += " AND ip_status = ?"
                params.append(ip_status)
            
            query += " ORDER BY discovery_date DESC LIMIT ? OFFSET ?"
            params.extend([limit, offset])
            
            cursor.execute(query, params)
            
            logs = []
            for row in cursor.fetchall():
                log = DiscoveryLog(
                    id=row[0],
                    compound_name=row[1],
                    smiles=row[2],
                    category=row[3],
                    confidence_score=row[4],
                    ip_status=row[5],
                    discovery_date=row[6],
                    source_research_topic_id=row[7],
                    generated_by=row[8],
                    molecular_weight=row[9],
                    logp=row[10],
                    receptor_activity=json.loads(row[11]) if row[11] else {},
                    binding_affinity=row[12],
                    notes=row[13]
                )
                logs.append(log)
            
            return logs
    
    def get_references_for_discovery(self, discovery_log_id: int) -> List[ResearchReference]:
        """Get all research references for a specific discovery."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute("SELECT * FROM research_references WHERE discovery_log_id = ?", (discovery_log_id,))
            
            refs = []
            for row in cursor.fetchall():
                ref = ResearchReference(
                    id=row[0],
                    doi=row[1],
                    title=row[2],
                    authors=json.loads(row[3]) if row[3] else [],
                    journal=row[4],
                    publication_year=row[5],
                    discovery_log_id=row[6],
                    retrieval_date=row[7]
                )
                refs.append(ref)
            
            return refs

# Global instance
discovery_logger = DiscoveryLogger()

def get_discovery_logger():
    """Get the global discovery logger instance."""
    return discovery_logger

# Test function
if __name__ == "__main__":
    logger = DiscoveryLogger()
    
    # Test logging a discovery
    log_id = logger.log_discovery(
        compound_name="TestCompound-001",
        smiles="C1=CC=C(C=C1)C(=O)O",
        category="Test-2024-A1",
        confidence_score=0.95,
        ip_status="Potential Opportunity",
        source_research_topic_id=1,
        molecular_weight=122.12,
        logp=1.87
    )
    
    print(f"Logged discovery with ID: {log_id}")
    
    # Test logging a reference
    ref_id = logger.log_research_reference(
        doi="10.1021/acs.jcim.0c00123",
        discovery_log_id=log_id,
        title="Test Article",
        authors=["Test Author"],
        journal="Test Journal",
        publication_year=2024
    )
    
    print(f"Logged reference with ID: {ref_id}")
    
    # Test getting discoveries
    discoveries = logger.get_discoveries(confidence_threshold=0.9)
    print(f"High-confidence discoveries: {len(discoveries)}")
    
    for discovery in discoveries:
        print(f"- {discovery.compound_name} (SMILES: {discovery.smiles}, Score: {discovery.confidence_score})")
        refs = logger.get_references_for_discovery(discovery.id)
        for ref in refs:
            print(f"  - Ref: {ref.doi}")
