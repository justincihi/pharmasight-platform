"""
Administrator Research Topic Management System for PharmaSight™
Allows dynamic creation and management of research topics/goals
"""

import sqlite3
import json
import datetime
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
import os

@dataclass
class ResearchTopic:
    """Research topic data structure."""
    id: int
    title: str
    description: str
    keywords: List[str]
    priority: int  # 1-5 scale
    status: str  # 'active', 'completed', 'paused'
    created_by: str
    created_date: str
    last_modified: str
    target_compounds: List[str]
    confidence_threshold: float
    ip_focus: bool  # Focus on IP opportunities
    therapeutic_areas: List[str]

class ResearchTopicManager:
    """Manages research topics with SQLite database persistence."""
    
    def __init__(self, db_path: str = "/home/ubuntu/pharmasight-platform/data/research_topics.db"):
        self.db_path = db_path
        self.ensure_data_directory()
        self.init_database()
    
    def ensure_data_directory(self):
        """Ensure the data directory exists."""
        data_dir = os.path.dirname(self.db_path)
        os.makedirs(data_dir, exist_ok=True)
    
    def init_database(self):
        """Initialize the SQLite database with research topics table."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Create research topics table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS research_topics (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    title TEXT NOT NULL,
                    description TEXT NOT NULL,
                    keywords TEXT NOT NULL,  -- JSON array
                    priority INTEGER DEFAULT 3,
                    status TEXT DEFAULT 'active',
                    created_by TEXT NOT NULL,
                    created_date TEXT NOT NULL,
                    last_modified TEXT NOT NULL,
                    target_compounds TEXT,  -- JSON array
                    confidence_threshold REAL DEFAULT 0.8,
                    ip_focus BOOLEAN DEFAULT 1,
                    therapeutic_areas TEXT  -- JSON array
                )
            ''')
            
            # Create research topic history table for audit trail
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS research_topic_history (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    topic_id INTEGER,
                    action TEXT NOT NULL,  -- 'created', 'updated', 'deleted'
                    changes TEXT,  -- JSON of changes
                    modified_by TEXT NOT NULL,
                    timestamp TEXT NOT NULL,
                    FOREIGN KEY (topic_id) REFERENCES research_topics (id)
                )
            ''')
            
            conn.commit()
            
            # Initialize with default topics if empty
            cursor.execute('SELECT COUNT(*) FROM research_topics')
            if cursor.fetchone()[0] == 0:
                self._create_default_topics()
    
    def _create_default_topics(self):
        """Create default research topics."""
        default_topics = [
            {
                'title': 'Novel Psychedelic Therapeutics',
                'description': 'Research novel psychedelic compounds for mental health applications with focus on patent-free opportunities',
                'keywords': ['psychedelics', 'mental health', 'serotonin', '5-HT2A', 'depression', 'PTSD'],
                'priority': 5,
                'therapeutic_areas': ['Mental Health', 'Neurology'],
                'target_compounds': ['psilocybin', 'LSD', 'DMT', 'mescaline'],
                'confidence_threshold': 0.85,
                'ip_focus': True
            },
            {
                'title': 'Ketamine Analog Development',
                'description': 'Develop ketamine analogs with improved safety profiles and reduced abuse potential',
                'keywords': ['ketamine', 'NMDA', 'dissociative', 'anesthesia', 'depression'],
                'priority': 4,
                'therapeutic_areas': ['Anesthesia', 'Mental Health'],
                'target_compounds': ['ketamine', 'esketamine', 'arketamine'],
                'confidence_threshold': 0.80,
                'ip_focus': True
            },
            {
                'title': 'MDMA-Inspired Therapeutics',
                'description': 'Research MDMA-like compounds for PTSD therapy with reduced neurotoxicity',
                'keywords': ['MDMA', 'empathogen', 'PTSD', 'serotonin', 'dopamine', 'therapy'],
                'priority': 4,
                'therapeutic_areas': ['Mental Health', 'Trauma Therapy'],
                'target_compounds': ['MDMA', 'MDA', '6-APB'],
                'confidence_threshold': 0.90,
                'ip_focus': True
            },
            {
                'title': 'Novel Antidepressant Mechanisms',
                'description': 'Explore alternative mechanisms for antidepressant action beyond traditional monoamine targets',
                'keywords': ['antidepressant', 'novel mechanisms', 'glutamate', 'GABA', 'neuroplasticity'],
                'priority': 3,
                'therapeutic_areas': ['Mental Health', 'Neurology'],
                'target_compounds': ['sertraline', 'fluoxetine', 'venlafaxine'],
                'confidence_threshold': 0.75,
                'ip_focus': False
            }
        ]
        
        for topic_data in default_topics:
            self.create_topic(
                title=topic_data['title'],
                description=topic_data['description'],
                keywords=topic_data['keywords'],
                priority=topic_data['priority'],
                created_by='system',
                target_compounds=topic_data['target_compounds'],
                confidence_threshold=topic_data['confidence_threshold'],
                ip_focus=topic_data['ip_focus'],
                therapeutic_areas=topic_data['therapeutic_areas']
            )
    
    def create_topic(self, title: str, description: str, keywords: List[str], 
                    priority: int = 3, created_by: str = 'admin',
                    target_compounds: List[str] = None, confidence_threshold: float = 0.8,
                    ip_focus: bool = True, therapeutic_areas: List[str] = None) -> int:
        """Create a new research topic."""
        
        if target_compounds is None:
            target_compounds = []
        if therapeutic_areas is None:
            therapeutic_areas = []
            
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT INTO research_topics 
                (title, description, keywords, priority, created_by, created_date, 
                 last_modified, target_compounds, confidence_threshold, ip_focus, therapeutic_areas)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                title, description, json.dumps(keywords), priority, created_by,
                timestamp, timestamp, json.dumps(target_compounds),
                confidence_threshold, ip_focus, json.dumps(therapeutic_areas)
            ))
            
            topic_id = cursor.lastrowid
            
            # Log the creation
            cursor.execute('''
                INSERT INTO research_topic_history 
                (topic_id, action, changes, modified_by, timestamp)
                VALUES (?, ?, ?, ?, ?)
            ''', (
                topic_id, 'created', json.dumps({
                    'title': title,
                    'description': description,
                    'keywords': keywords,
                    'priority': priority
                }), created_by, timestamp
            ))
            
            conn.commit()
            return topic_id
    
    def get_all_topics(self, status: str = None) -> List[ResearchTopic]:
        """Get all research topics, optionally filtered by status."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            if status:
                cursor.execute('SELECT * FROM research_topics WHERE status = ? ORDER BY priority DESC, created_date DESC', (status,))
            else:
                cursor.execute('SELECT * FROM research_topics ORDER BY priority DESC, created_date DESC')
            
            topics = []
            for row in cursor.fetchall():
                topic = ResearchTopic(
                    id=row[0],
                    title=row[1],
                    description=row[2],
                    keywords=json.loads(row[3]),
                    priority=row[4],
                    status=row[5],
                    created_by=row[6],
                    created_date=row[7],
                    last_modified=row[8],
                    target_compounds=json.loads(row[9]) if row[9] else [],
                    confidence_threshold=row[10],
                    ip_focus=bool(row[11]),
                    therapeutic_areas=json.loads(row[12]) if row[12] else []
                )
                topics.append(topic)
            
            return topics
    
    def get_topic_by_id(self, topic_id: int) -> Optional[ResearchTopic]:
        """Get a specific research topic by ID."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT * FROM research_topics WHERE id = ?', (topic_id,))
            row = cursor.fetchone()
            
            if row:
                return ResearchTopic(
                    id=row[0],
                    title=row[1],
                    description=row[2],
                    keywords=json.loads(row[3]),
                    priority=row[4],
                    status=row[5],
                    created_by=row[6],
                    created_date=row[7],
                    last_modified=row[8],
                    target_compounds=json.loads(row[9]) if row[9] else [],
                    confidence_threshold=row[10],
                    ip_focus=bool(row[11]),
                    therapeutic_areas=json.loads(row[12]) if row[12] else []
                )
            return None
    
    def update_topic(self, topic_id: int, **kwargs) -> bool:
        """Update a research topic."""
        topic = self.get_topic_by_id(topic_id)
        if not topic:
            return False
        
        # Track changes
        changes = {}
        update_fields = []
        update_values = []
        
        for field, value in kwargs.items():
            if hasattr(topic, field):
                old_value = getattr(topic, field)
                if old_value != value:
                    changes[field] = {'old': old_value, 'new': value}
                    
                    if field in ['keywords', 'target_compounds', 'therapeutic_areas']:
                        update_fields.append(f"{field} = ?")
                        update_values.append(json.dumps(value))
                    else:
                        update_fields.append(f"{field} = ?")
                        update_values.append(value)
        
        if not changes:
            return True  # No changes needed
        
        # Add last_modified timestamp
        timestamp = datetime.datetime.now().isoformat()
        update_fields.append("last_modified = ?")
        update_values.append(timestamp)
        update_values.append(topic_id)
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Update the topic
            query = f"UPDATE research_topics SET {', '.join(update_fields)} WHERE id = ?"
            cursor.execute(query, update_values)
            
            # Log the update
            cursor.execute('''
                INSERT INTO research_topic_history 
                (topic_id, action, changes, modified_by, timestamp)
                VALUES (?, ?, ?, ?, ?)
            ''', (
                topic_id, 'updated', json.dumps(changes),
                kwargs.get('modified_by', 'admin'), timestamp
            ))
            
            conn.commit()
            return True
    
    def delete_topic(self, topic_id: int, deleted_by: str = 'admin') -> bool:
        """Delete a research topic (soft delete by changing status)."""
        topic = self.get_topic_by_id(topic_id)
        if not topic:
            return False
        
        timestamp = datetime.datetime.now().isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Soft delete by changing status
            cursor.execute('''
                UPDATE research_topics 
                SET status = 'deleted', last_modified = ?
                WHERE id = ?
            ''', (timestamp, topic_id))
            
            # Log the deletion
            cursor.execute('''
                INSERT INTO research_topic_history 
                (topic_id, action, changes, modified_by, timestamp)
                VALUES (?, ?, ?, ?, ?)
            ''', (
                topic_id, 'deleted', json.dumps({'status': 'deleted'}),
                deleted_by, timestamp
            ))
            
            conn.commit()
            return True
    
    def search_topics(self, query: str) -> List[ResearchTopic]:
        """Search topics by title, description, or keywords."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            
            search_pattern = f"%{query.lower()}%"
            cursor.execute('''
                SELECT * FROM research_topics 
                WHERE status != 'deleted' AND (
                    LOWER(title) LIKE ? OR 
                    LOWER(description) LIKE ? OR 
                    LOWER(keywords) LIKE ?
                )
                ORDER BY priority DESC, created_date DESC
            ''', (search_pattern, search_pattern, search_pattern))
            
            topics = []
            for row in cursor.fetchall():
                topic = ResearchTopic(
                    id=row[0],
                    title=row[1],
                    description=row[2],
                    keywords=json.loads(row[3]),
                    priority=row[4],
                    status=row[5],
                    created_by=row[6],
                    created_date=row[7],
                    last_modified=row[8],
                    target_compounds=json.loads(row[9]) if row[9] else [],
                    confidence_threshold=row[10],
                    ip_focus=bool(row[11]),
                    therapeutic_areas=json.loads(row[12]) if row[12] else []
                )
                topics.append(topic)
            
            return topics
    
    def get_topic_history(self, topic_id: int) -> List[Dict[str, Any]]:
        """Get the history of changes for a topic."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT action, changes, modified_by, timestamp 
                FROM research_topic_history 
                WHERE topic_id = ? 
                ORDER BY timestamp DESC
            ''', (topic_id,))
            
            history = []
            for row in cursor.fetchall():
                history.append({
                    'action': row[0],
                    'changes': json.loads(row[1]) if row[1] else {},
                    'modified_by': row[2],
                    'timestamp': row[3]
                })
            
            return history
    
    def get_active_research_goals(self) -> List[Dict[str, Any]]:
        """Get active research goals formatted for the research engine."""
        active_topics = self.get_all_topics(status='active')
        
        goals = []
        for topic in active_topics:
            goal = {
                'id': topic.id,
                'title': topic.title,
                'description': topic.description,
                'keywords': topic.keywords,
                'priority': topic.priority,
                'target_compounds': topic.target_compounds,
                'confidence_threshold': topic.confidence_threshold,
                'ip_focus': topic.ip_focus,
                'therapeutic_areas': topic.therapeutic_areas,
                'search_criteria': {
                    'min_confidence': topic.confidence_threshold,
                    'focus_areas': topic.therapeutic_areas,
                    'patent_status': 'patent_free' if topic.ip_focus else 'any'
                }
            }
            goals.append(goal)
        
        return goals

# Global instance
research_topic_manager = ResearchTopicManager()

def get_research_topic_manager():
    """Get the global research topic manager instance."""
    return research_topic_manager

# Test function
if __name__ == "__main__":
    manager = ResearchTopicManager()
    
    # Test creating a topic
    topic_id = manager.create_topic(
        title="Test Novel Compounds",
        description="Testing the research topic management system",
        keywords=["test", "novel", "compounds"],
        priority=2,
        created_by="test_user"
    )
    
    print(f"Created topic with ID: {topic_id}")
    
    # Test getting all topics
    topics = manager.get_all_topics()
    print(f"Total topics: {len(topics)}")
    
    for topic in topics:
        print(f"- {topic.title} (Priority: {topic.priority}, Status: {topic.status})")
