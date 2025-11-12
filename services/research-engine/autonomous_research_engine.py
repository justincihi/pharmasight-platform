"""
Autonomous Research Engine for PharmaSight™
Performs daily research activities with rate limiting and cost management
Screens for novel drug targets using chemical similarity and patent landscape
"""

import json
import time
import requests
from datetime import datetime, timedelta
from typing import List, Dict, Optional
import sys
sys.path.append('/home/ubuntu/pharmasight-latest/src')

from research_article_database import ResearchArticleDatabase
from rdkit_analog_generator import generate_novel_analogs

class AutonomousResearchEngine:
    """
    Autonomous research engine that performs daily literature searches,
    screens for novel drug targets, and manages IP opportunities
    """
    
    def __init__(self, 
                 max_api_calls_per_day=50,
                 max_articles_per_search=10):
        self.article_db = ResearchArticleDatabase()
        self.max_api_calls_per_day = max_api_calls_per_day
        self.max_articles_per_search = max_articles_per_search
        self.api_call_count = 0
        self.session_log = {
            "session_start": datetime.now().isoformat(),
            "api_calls_made": 0,
            "articles_found": 0,
            "analogs_generated": 0,
            "errors": []
        }
    
    def run_daily_research_cycle(self, research_goals: List[str]) -> Dict:
        """
        Run a complete daily research cycle
        
        Args:
            research_goals: List of research topics to investigate
        
        Returns:
            Summary of research activities
        """
        print("=" * 80)
        print("PHARMASIGHT™ AUTONOMOUS RESEARCH ENGINE - DAILY CYCLE")
        print("=" * 80)
        print(f"Session Start: {self.session_log['session_start']}")
        print(f"Max API Calls: {self.max_api_calls_per_day}")
        print(f"Research Goals: {len(research_goals)}")
        print("=" * 80)
        
        for goal in research_goals:
            if self.api_call_count >= self.max_api_calls_per_day:
                print(f"\n⚠️  API call limit reached ({self.max_api_calls_per_day})")
                break
            
            print(f"\n🔬 Researching: {goal}")
            self._research_goal(goal)
            
            # Rate limiting: 2 seconds between searches
            time.sleep(2)
        
        # Generate summary
        summary = self._generate_session_summary()
        self._save_session_log(summary)
        
        return summary
    
    def _research_goal(self, goal: str):
        """Research a specific goal"""
        try:
            # Search PubMed for relevant articles
            articles = self._search_pubmed(goal)
            
            if articles:
                print(f"   ✅ Found {len(articles)} articles")
                
                # Add to article database
                for article in articles:
                    self.article_db.add_article(
                        title=article['title'],
                        authors=article.get('authors', []),
                        doi=article.get('doi'),
                        pmid=article.get('pmid'),
                        journal=article.get('journal'),
                        year=article.get('year'),
                        abstract=article.get('abstract'),
                        keywords=article.get('keywords', []),
                        relevance_score=article.get('relevance_score', 80.0),
                        source_database="PubMed",
                        research_goal=goal
                    )
                    self.session_log['articles_found'] += 1
                
                # Save article database
                self.article_db.save()
                print(f"   📁 Articles added to database")
            else:
                print(f"   ℹ️  No new articles found")
        
        except Exception as e:
            print(f"   ❌ Error: {str(e)}")
            self.session_log['errors'].append({
                "goal": goal,
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            })
    
    def _search_pubmed(self, query: str, max_results: int = None) -> List[Dict]:
        """
        Search PubMed for articles (rate-limited)
        
        Args:
            query: Search query
            max_results: Maximum number of results (default: self.max_articles_per_search)
        
        Returns:
            List of article dictionaries
        """
        if max_results is None:
            max_results = self.max_articles_per_search
        
        # Check API call limit
        if self.api_call_count >= self.max_api_calls_per_day:
            return []
        
        try:
            # Use NCBI E-utilities API (free, no API key required for basic usage)
            # Search for article IDs
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "retmode": "json",
                "sort": "relevance"
            }
            
            self.api_call_count += 1
            self.session_log['api_calls_made'] += 1
            
            response = requests.get(search_url, params=search_params, timeout=10)
            response.raise_for_status()
            
            search_data = response.json()
            pmids = search_data.get('esearchresult', {}).get('idlist', [])
            
            if not pmids:
                return []
            
            # Rate limiting: 3 requests per second for NCBI
            time.sleep(0.4)
            
            # Fetch article details
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(pmids),
                "retmode": "json"
            }
            
            self.api_call_count += 1
            self.session_log['api_calls_made'] += 1
            
            response = requests.get(fetch_url, params=fetch_params, timeout=10)
            response.raise_for_status()
            
            fetch_data = response.json()
            
            articles = []
            for pmid in pmids:
                if pmid in fetch_data.get('result', {}):
                    article_data = fetch_data['result'][pmid]
                    
                    # Extract authors
                    authors = []
                    if 'authors' in article_data:
                        authors = [author.get('name', '') for author in article_data['authors'][:5]]
                    
                    # Extract DOI
                    doi = None
                    if 'articleids' in article_data:
                        for article_id in article_data['articleids']:
                            if article_id.get('idtype') == 'doi':
                                doi = article_id.get('value')
                                break
                    
                    article = {
                        "title": article_data.get('title', 'Unknown Title'),
                        "authors": authors,
                        "pmid": pmid,
                        "doi": doi,
                        "journal": article_data.get('fulljournalname', ''),
                        "year": int(article_data.get('pubdate', '').split()[0]) if article_data.get('pubdate') else None,
                        "abstract": None,  # Would require additional API call
                        "keywords": [],
                        "relevance_score": 85.0,  # Default relevance
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    }
                    articles.append(article)
            
            return articles
        
        except requests.exceptions.RequestException as e:
            print(f"   ⚠️  PubMed API error: {str(e)}")
            return []
        except Exception as e:
            print(f"   ⚠️  Error parsing PubMed results: {str(e)}")
            return []
    
    def screen_for_novel_targets(self, compound_smiles: str, compound_name: str) -> Dict:
        """
        Screen for novel drug targets using chemical similarity
        Focus on patent-free opportunities
        
        Args:
            compound_smiles: SMILES string of parent compound
            compound_name: Name of parent compound
        
        Returns:
            Dictionary with screening results
        """
        print(f"\n🎯 Screening for novel targets: {compound_name}")
        
        results = {
            "compound_name": compound_name,
            "compound_smiles": compound_smiles,
            "analogs_generated": 0,
            "patent_free_analogs": 0,
            "high_ip_opportunity": 0,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            # Generate analogs
            analogs = generate_novel_analogs(
                parent_smiles=compound_smiles,
                parent_name=compound_name,
                num_analogs=10
            )
            
            results['analogs_generated'] = len(analogs)
            self.session_log['analogs_generated'] += len(analogs)
            
            # Count patent-free and high IP opportunity analogs
            for analog in analogs:
                if analog.get('patent_status') == 'Patent-Free (Novel)':
                    results['patent_free_analogs'] += 1
                if analog.get('patent_opportunity_score', 0) >= 90:
                    results['high_ip_opportunity'] += 1
            
            print(f"   ✅ Generated {results['analogs_generated']} analogs")
            print(f"   📊 Patent-Free: {results['patent_free_analogs']}")
            print(f"   💎 High IP Opportunity: {results['high_ip_opportunity']}")
        
        except Exception as e:
            print(f"   ❌ Error: {str(e)}")
            results['error'] = str(e)
        
        return results
    
    def _generate_session_summary(self) -> Dict:
        """Generate summary of research session"""
        self.session_log['session_end'] = datetime.now().isoformat()
        
        # Calculate duration
        start = datetime.fromisoformat(self.session_log['session_start'])
        end = datetime.fromisoformat(self.session_log['session_end'])
        duration = (end - start).total_seconds()
        
        self.session_log['duration_seconds'] = duration
        self.session_log['api_calls_remaining'] = self.max_api_calls_per_day - self.api_call_count
        
        return self.session_log
    
    def _save_session_log(self, summary: Dict):
        """Save session log to file"""
        log_path = f"/home/ubuntu/pharmasight-latest/research_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        
        with open(log_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"\n📁 Session log saved: {log_path}")
    
    def print_session_summary(self):
        """Print session summary"""
        print("\n" + "=" * 80)
        print("SESSION SUMMARY")
        print("=" * 80)
        print(f"Duration: {self.session_log.get('duration_seconds', 0):.1f} seconds")
        print(f"API Calls Made: {self.session_log['api_calls_made']}/{self.max_api_calls_per_day}")
        print(f"Articles Found: {self.session_log['articles_found']}")
        print(f"Analogs Generated: {self.session_log['analogs_generated']}")
        print(f"Errors: {len(self.session_log['errors'])}")
        print("=" * 80)

# Default research goals
DEFAULT_RESEARCH_GOALS = [
    "psilocybin 5-HT2A receptor depression treatment",
    "ketamine NMDA antagonist rapid antidepressant",
    "MDMA PTSD therapy serotonin oxytocin",
    "kava lactones anxiolytic GABA modulation",
    "muscimol GABA-A agonist sleep anxiety"
]

def run_daily_cycle():
    """Run a daily research cycle"""
    engine = AutonomousResearchEngine(
        max_api_calls_per_day=50,  # Conservative limit
        max_articles_per_search=5   # 5 articles per goal = 25 total
    )
    
    summary = engine.run_daily_research_cycle(DEFAULT_RESEARCH_GOALS)
    engine.print_session_summary()
    
    return summary

if __name__ == "__main__":
    print("🚀 Starting Autonomous Research Engine...")
    summary = run_daily_cycle()
    print("\n✅ Daily research cycle complete!")

