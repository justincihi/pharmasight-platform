"""
Confidence-Based Data Organization and DOI Tracking System for PharmaSight™
Provides regulatory-compliant organization of discoveries by confidence levels and comprehensive DOI tracking.
"""

import sqlite3
import json
import csv
import datetime
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
import os
from discovery_logging_system import get_discovery_logger

@dataclass
class ConfidenceDataset:
    """Confidence-based dataset structure."""
    confidence_threshold: float
    total_compounds: int
    compounds: List[Dict[str, Any]]
    avg_confidence: float
    ip_opportunities: int
    therapeutic_areas: Dict[str, int]
    generated_date: str

@dataclass
class DOIReference:
    """DOI reference with comprehensive metadata."""
    doi: str
    title: str
    authors: List[str]
    journal: str
    publication_year: int
    abstract: str
    keywords: List[str]
    citation_count: int
    impact_factor: float
    research_area: str
    compound_relevance: str
    retrieval_date: str
    confidence_score: float

class ConfidenceDOIManager:
    """Manages confidence-based data organization and DOI tracking."""
    
    def __init__(self, db_path: str = "/home/ubuntu/pharmasight-platform/data/confidence_doi.db"):
        self.db_path = db_path
        self.ensure_data_directory()
        self.init_database()
        self.discovery_logger = get_discovery_logger()
    
    def ensure_data_directory(self):
        """Ensure the data directory exists."""
        data_dir = os.path.dirname(self.db_path)
        os.makedirs(data_dir, exist_ok=True)
    
    def init_database(self):
        """Initialize the SQLite database with confidence and DOI tracking tables."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Create DOI references table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS doi_references (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    doi TEXT NOT NULL UNIQUE,
                    title TEXT,
                    authors TEXT,  -- JSON array
                    journal TEXT,
                    publication_year INTEGER,
                    abstract TEXT,
                    keywords TEXT,  -- JSON array
                    citation_count INTEGER DEFAULT 0,
                    impact_factor REAL DEFAULT 0.0,
                    research_area TEXT,
                    compound_relevance TEXT,
                    retrieval_date TEXT NOT NULL,
                    confidence_score REAL DEFAULT 0.8
                )
            """)
            
            # Create confidence datasets table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS confidence_datasets (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    confidence_threshold REAL NOT NULL,
                    total_compounds INTEGER,
                    avg_confidence REAL,
                    ip_opportunities INTEGER,
                    therapeutic_areas TEXT,  -- JSON object
                    generated_date TEXT NOT NULL,
                    dataset_type TEXT DEFAULT 'standard'  -- 'high_confidence', 'standard', 'all'
                )
            """)
            
            # Create compound-DOI associations table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS compound_doi_associations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    compound_id INTEGER,
                    doi_id INTEGER,
                    relevance_score REAL DEFAULT 0.8,
                    association_type TEXT,  -- 'discovery', 'validation', 'reference'
                    created_date TEXT NOT NULL,
                    FOREIGN KEY (compound_id) REFERENCES discovery_logs (id),
                    FOREIGN KEY (doi_id) REFERENCES doi_references (id)
                )
            """)
            
            conn.commit()
    
    def add_doi_reference(self, doi: str, title: str = None, authors: List[str] = None,
                         journal: str = None, publication_year: int = None,
                         abstract: str = None, keywords: List[str] = None,
                         citation_count: int = 0, impact_factor: float = 0.0,
                         research_area: str = None, compound_relevance: str = None,
                         confidence_score: float = 0.8) -> int:
        """Add a DOI reference to the tracking system."""
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            try:
                cursor.execute("""
                    INSERT INTO doi_references 
                    (doi, title, authors, journal, publication_year, abstract, keywords,
                     citation_count, impact_factor, research_area, compound_relevance,
                     retrieval_date, confidence_score)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    doi, title, json.dumps(authors) if authors else None,
                    journal, publication_year, abstract, json.dumps(keywords) if keywords else None,
                    citation_count, impact_factor, research_area, compound_relevance,
                    timestamp, confidence_score
                ))
                
                doi_id = cursor.lastrowid
                conn.commit()
                return doi_id
            
            except sqlite3.IntegrityError:
                # DOI already exists, return existing ID
                cursor.execute("SELECT id FROM doi_references WHERE doi = ?", (doi,))
                return cursor.fetchone()[0]
    
    def associate_compound_doi(self, compound_id: int, doi_id: int, 
                              relevance_score: float = 0.8, 
                              association_type: str = "reference") -> int:
        """Associate a compound with a DOI reference."""
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute("""
                INSERT INTO compound_doi_associations 
                (compound_id, doi_id, relevance_score, association_type, created_date)
                VALUES (?, ?, ?, ?, ?)
            """, (compound_id, doi_id, relevance_score, association_type, timestamp))
            
            assoc_id = cursor.lastrowid
            conn.commit()
            return assoc_id
    
    def generate_confidence_dataset(self, confidence_threshold: float = 0.9, 
                                   dataset_type: str = "high_confidence") -> ConfidenceDataset:
        """Generate a confidence-based dataset."""
        discoveries = self.discovery_logger.get_discoveries(confidence_threshold=confidence_threshold)
        
        # Convert discoveries to dict format
        compounds = []
        therapeutic_areas = {}
        ip_opportunities = 0
        total_confidence = 0
        
        for discovery in discoveries:
            compound_dict = {
                'id': discovery.id,
                'name': discovery.compound_name,
                'smiles': discovery.smiles,
                'category': discovery.category,
                'confidence_score': discovery.confidence_score,
                'ip_status': discovery.ip_status,
                'discovery_date': discovery.discovery_date,
                'molecular_weight': discovery.molecular_weight,
                'logp': discovery.logp,
                'receptor_activity': discovery.receptor_activity,
                'binding_affinity': discovery.binding_affinity,
                'notes': discovery.notes
            }
            compounds.append(compound_dict)
            
            # Track therapeutic areas (would come from research topic associations)
            # For now, infer from compound category
            if discovery.category:
                area = discovery.category.split('-')[0] if '-' in discovery.category else 'Unknown'
                therapeutic_areas[area] = therapeutic_areas.get(area, 0) + 1
            
            # Count IP opportunities
            if discovery.ip_status and 'opportunity' in discovery.ip_status.lower():
                ip_opportunities += 1
            
            total_confidence += discovery.confidence_score
        
        avg_confidence = total_confidence / len(discoveries) if discoveries else 0
        
        dataset = ConfidenceDataset(
            confidence_threshold=confidence_threshold,
            total_compounds=len(compounds),
            compounds=compounds,
            avg_confidence=avg_confidence,
            ip_opportunities=ip_opportunities,
            therapeutic_areas=therapeutic_areas,
            generated_date=datetime.datetime.now().isoformat()
        )
        
        # Store dataset metadata
        self._store_dataset_metadata(dataset, dataset_type)
        
        return dataset
    
    def _store_dataset_metadata(self, dataset: ConfidenceDataset, dataset_type: str):
        """Store dataset metadata in the database."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute("""
                INSERT INTO confidence_datasets 
                (confidence_threshold, total_compounds, avg_confidence, ip_opportunities,
                 therapeutic_areas, generated_date, dataset_type)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                dataset.confidence_threshold, dataset.total_compounds, dataset.avg_confidence,
                dataset.ip_opportunities, json.dumps(dataset.therapeutic_areas),
                dataset.generated_date, dataset_type
            ))
            
            conn.commit()
    
    def export_confidence_dataset_json(self, confidence_threshold: float = 0.9, 
                                      output_path: str = None) -> str:
        """Export confidence dataset to JSON file."""
        dataset = self.generate_confidence_dataset(confidence_threshold)
        
        if output_path is None:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"/home/ubuntu/pharmasight-platform/data/confidence_dataset_{confidence_threshold}_{timestamp}.json"
        
        with open(output_path, 'w') as f:
            json.dump(asdict(dataset), f, indent=2)
        
        return output_path
    
    def export_confidence_dataset_csv(self, confidence_threshold: float = 0.9, 
                                     output_path: str = None) -> str:
        """Export confidence dataset to CSV file."""
        dataset = self.generate_confidence_dataset(confidence_threshold)
        
        if output_path is None:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"/home/ubuntu/pharmasight-platform/data/confidence_dataset_{confidence_threshold}_{timestamp}.csv"
        
        with open(output_path, 'w', newline='') as f:
            if dataset.compounds:
                fieldnames = dataset.compounds[0].keys()
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(dataset.compounds)
        
        return output_path
    
    def get_doi_references(self, research_area: str = None, 
                          min_confidence: float = None) -> List[DOIReference]:
        """Get DOI references with optional filtering."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            query = "SELECT * FROM doi_references WHERE 1=1"
            params = []
            
            if research_area:
                query += " AND research_area = ?"
                params.append(research_area)
            
            if min_confidence:
                query += " AND confidence_score >= ?"
                params.append(min_confidence)
            
            query += " ORDER BY confidence_score DESC, publication_year DESC"
            
            cursor.execute(query, params)
            
            references = []
            for row in cursor.fetchall():
                ref = DOIReference(
                    doi=row[1],
                    title=row[2],
                    authors=json.loads(row[3]) if row[3] else [],
                    journal=row[4],
                    publication_year=row[5],
                    abstract=row[6],
                    keywords=json.loads(row[7]) if row[7] else [],
                    citation_count=row[8],
                    impact_factor=row[9],
                    research_area=row[10],
                    compound_relevance=row[11],
                    retrieval_date=row[12],
                    confidence_score=row[13]
                )
                references.append(ref)
            
            return references
    
    def generate_regulatory_report(self, confidence_threshold: float = 0.9) -> Dict[str, Any]:
        """Generate a comprehensive regulatory compliance report."""
        high_confidence_dataset = self.generate_confidence_dataset(confidence_threshold, "high_confidence")
        all_dataset = self.generate_confidence_dataset(0.0, "all")
        doi_references = self.get_doi_references()
        
        report = {
            'report_metadata': {
                'generated_date': datetime.datetime.now().isoformat(),
                'generated_by': 'PharmaSight™ AI Engine',
                'confidence_threshold': confidence_threshold,
                'total_references': len(doi_references)
            },
            'high_confidence_discoveries': {
                'threshold': confidence_threshold,
                'total_compounds': high_confidence_dataset.total_compounds,
                'avg_confidence': high_confidence_dataset.avg_confidence,
                'ip_opportunities': high_confidence_dataset.ip_opportunities,
                'therapeutic_areas': high_confidence_dataset.therapeutic_areas,
                'compounds': high_confidence_dataset.compounds
            },
            'all_discoveries': {
                'total_compounds': all_dataset.total_compounds,
                'avg_confidence': all_dataset.avg_confidence,
                'ip_opportunities': all_dataset.ip_opportunities,
                'therapeutic_areas': all_dataset.therapeutic_areas
            },
            'research_references': [
                {
                    'doi': ref.doi,
                    'title': ref.title,
                    'authors': ref.authors,
                    'journal': ref.journal,
                    'publication_year': ref.publication_year,
                    'research_area': ref.research_area,
                    'confidence_score': ref.confidence_score
                }
                for ref in doi_references
            ],
            'compliance_summary': {
                'total_discoveries_logged': all_dataset.total_compounds,
                'high_confidence_discoveries': high_confidence_dataset.total_compounds,
                'total_references_tracked': len(doi_references),
                'ip_opportunities_identified': high_confidence_dataset.ip_opportunities,
                'regulatory_readiness': 'Complete' if high_confidence_dataset.total_compounds > 0 else 'Pending'
            }
        }
        
        return report
    
    def export_regulatory_report(self, confidence_threshold: float = 0.9, 
                                output_path: str = None) -> str:
        """Export regulatory compliance report to JSON file."""
        report = self.generate_regulatory_report(confidence_threshold)
        
        if output_path is None:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"/home/ubuntu/pharmasight-platform/data/regulatory_report_{timestamp}.json"
        
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        return output_path

# Global instance
confidence_doi_manager = ConfidenceDOIManager()

def get_confidence_doi_manager():
    """Get the global confidence DOI manager instance."""
    return confidence_doi_manager

# Test function
if __name__ == "__main__":
    manager = ConfidenceDOIManager()
    
    # Test adding DOI references
    doi_id = manager.add_doi_reference(
        doi="10.1021/acs.jcim.2024.0123",
        title="Novel Psychedelic Compounds for Mental Health Applications",
        authors=["Smith, J.", "Johnson, A.", "Brown, K."],
        journal="Journal of Chemical Information and Modeling",
        publication_year=2024,
        research_area="Psychedelic Research",
        compound_relevance="High",
        confidence_score=0.95
    )
    
    print(f"Added DOI reference with ID: {doi_id}")
    
    # Test generating confidence dataset
    dataset = manager.generate_confidence_dataset(0.8)
    print(f"Generated dataset with {dataset.total_compounds} compounds")
    print(f"Average confidence: {dataset.avg_confidence:.3f}")
    
    # Test generating regulatory report
    report = manager.generate_regulatory_report(0.9)
    print(f"Generated regulatory report with {len(report['research_references'])} references")
    print(f"Regulatory readiness: {report['compliance_summary']['regulatory_readiness']}")
    
    print("Confidence DOI Manager working successfully!")
