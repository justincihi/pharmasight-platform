"""
PharmaSight™ Research Engine Microservice
Autonomous research, article database, and RDKit integration
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import logging
from datetime import datetime

from autonomous_research_engine import AutonomousResearchEngine
from research_article_database import ResearchArticleDatabase
from research_rdkit_integration import ResearchRDKitIntegration

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="PharmaSight Research Engine",
    description="Autonomous research, literature scanning, and analog discovery",
    version="1.0.0"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize services
research_engine = AutonomousResearchEngine(
    max_api_calls_per_day=50,
    max_articles_per_search=5
)
article_db = ResearchArticleDatabase()
rdkit_integration = ResearchRDKitIntegration()

# Request/Response Models
class ResearchGoal(BaseModel):
    goal: str

class ResearchGoalsList(BaseModel):
    goals: List[str]

class ArticleSearchRequest(BaseModel):
    keyword: str
    limit: Optional[int] = 10

class RDKitSyncRequest(BaseModel):
    max_compounds: Optional[int] = 3

# Health check
@app.get("/@app.get("/health")
async def health():
    """Health check endpoint"""
    return {"status": "healthy", "service": "research-engine"}

# PubChem Integration
@app.get("/pubchem/compound/{compound_name}")
async def get_pubchem_compound(compound_name: str):
    """Get compound information from PubChem"""
    try:
        from api_integrations import PubChemAPI
        api = PubChemAPI()
        data = api.get_compound_by_name(compound_name)
        
        if data:
            return {"success": True, "compound": data}
        return {"success": False, "error": "Compound not found"}
    except Exception as e:
        return {"success": False, "error": str(e)}

@app.get("/pubchem/search/{query}")
async def search_pubchem(query: str, limit: int = 10):
    """Search compounds in PubChem"""
    try:
        from api_integrations import PubChemAPI
        api = PubChemAPI()
        results = api.search_compounds(query, limit=limit)
        
        return {"success": True, "results": results, "count": len(results)}
    except Exception as e:
        return {"success": False, "error": str(e)}rch-engine",
        "timestamp": datetime.now().isoformat()
    }

# Research Engine Endpoints
@app.post("/research/run-cycle")
async def run_research_cycle(request: ResearchGoalsList):
    """Run a daily research cycle with specified goals"""
    try:
        logger.info(f"Starting research cycle with {len(request.goals)} goals")
        summary = research_engine.run_daily_research_cycle(request.goals)
        return {
            "success": True,
            "summary": summary,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Research cycle failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/research/session-summary")
async def get_session_summary():
    """Get the current research session summary"""
    try:
        summary = research_engine.get_session_summary()
        return {
            "success": True,
            "summary": summary
        }
    except Exception as e:
        logger.error(f"Failed to get session summary: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# Article Database Endpoints
@app.get("/articles/all")
async def get_all_articles():
    """Get all articles from the database"""
    try:
        articles = article_db.get_all_articles()
        return {
            "success": True,
            "count": len(articles),
            "articles": articles
        }
    except Exception as e:
        logger.error(f"Failed to get articles: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/articles/search")
async def search_articles(request: ArticleSearchRequest):
    """Search articles by keyword"""
    try:
        articles = article_db.search_by_keyword(request.keyword, request.limit)
        return {
            "success": True,
            "count": len(articles),
            "articles": articles
        }
    except Exception as e:
        logger.error(f"Article search failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/articles/statistics")
async def get_article_statistics():
    """Get database statistics"""
    try:
        stats = article_db.get_statistics()
        return {
            "success": True,
            "statistics": stats
        }
    except Exception as e:
        logger.error(f"Failed to get statistics: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/articles/by-year/{year}")
async def get_articles_by_year(year: int):
    """Get articles published in a specific year"""
    try:
        articles = article_db.get_articles_by_year(year)
        return {
            "success": True,
            "year": year,
            "count": len(articles),
            "articles": articles
        }
    except Exception as e:
        logger.error(f"Failed to get articles by year: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# RDKit Integration Endpoints
@app.post("/rdkit/sync")
async def sync_research_with_rdkit(request: RDKitSyncRequest):
    """Sync research goals with RDKit for analog generation"""
    try:
        logger.info(f"Starting RDKit sync with max_compounds={request.max_compounds}")
        summary = rdkit_integration.sync_research_goals_with_rdkit(
            max_compounds=request.max_compounds
        )
        return {
            "success": True,
            "summary": summary,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"RDKit sync failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/rdkit/sync-summary")
async def get_rdkit_sync_summary():
    """Get the last RDKit sync summary"""
    try:
        # This would need to be implemented in the integration class
        return {
            "success": True,
            "message": "Check sync logs for details"
        }
    except Exception as e:
        logger.error(f"Failed to get sync summary: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# Analog Discovery Endpoints
@app.get("/discoveries/all")
async def get_all_discoveries():
    """Get all analog discoveries"""
    try:
        discoveries = rdkit_integration.get_all_discoveries()
        return {
            "success": True,
            "total_sessions": discoveries.get("total_sessions", 0),
            "total_analogs": discoveries.get("total_analogs", 0),
            "discoveries": discoveries.get("discovery_sessions", [])
        }
    except Exception as e:
        logger.error(f"Failed to get discoveries: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/discoveries/statistics")
async def get_discovery_statistics():
    """Get analog discovery statistics"""
    try:
        stats = rdkit_integration.get_discovery_statistics()
        return {
            "success": True,
            "statistics": stats
        }
    except Exception as e:
        logger.error(f"Failed to get discovery statistics: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8006)

