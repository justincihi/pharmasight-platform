"""
Research Article Database for PharmaSight™
Simple in-memory database for storing research articles
"""

from typing import List, Dict, Optional
from datetime import datetime

class ResearchArticleDatabase:
    """
    Simple in-memory database for research articles
    """
    
    def __init__(self):
        self.articles = []
    
    def add_article(self, article: Dict) -> bool:
        """Add an article to the database"""
        article['added_at'] = datetime.now().isoformat()
        self.articles.append(article)
        return True
    
    def get_articles(self, limit: int = 100) -> List[Dict]:
        """Get recent articles"""
        return self.articles[-limit:]
    
    def search_articles(self, query: str) -> List[Dict]:
        """Search articles by title or keywords"""
        results = []
        query_lower = query.lower()
        for article in self.articles:
            if query_lower in article.get('title', '').lower():
                results.append(article)
        return results
