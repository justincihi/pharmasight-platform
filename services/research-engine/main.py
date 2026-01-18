"""
PharmaSight™ Research Engine Microservice - Enhanced
Autonomous research, article database, RDKit integration, and PubChem API
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
    title="PharmaSight Research Engine - Enhanced",
    description="Autonomous research with PubChem integration, literature scanning, and analog discovery",
    version="2.0.0"
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
class ResearchGoalsList(BaseModel):
    goals: List[str]

class ArticleSearchRequest(BaseModel):
    keyword: str
    limit: Optional[int] = 10

class RDKitSyncRequest(BaseModel):
    max_compounds: Optional[int] = 3

# Health check
@app.get("/health")
async def health():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "research-engine-enhanced",
        "version": "2.0.0",
        "features": [
            "autonomous_research",
            "article_database",
            "rdkit_integration",
            "pubchem_api"
        ],
        "timestamp": datetime.now().isoformat()
    }

# PubChem Integration Endpoints
@app.get("/pubchem/compound/{compound_name}")
async def get_pubchem_compound(compound_name: str):
    """Get compound information from PubChem by name"""
    try:
        from api_integrations import PubChemAPI
        api = PubChemAPI()
        data = api.get_compound_by_name(compound_name)
        
        if data:
            return {
                "success": True,
                "compound": data,
                "source": "PubChem",
                "timestamp": datetime.now().isoformat()
            }
        return {"success": False, "error": "Compound not found in PubChem"}
    except Exception as e:
        logger.error(f"PubChem API error: {str(e)}")
        return {"success": False, "error": str(e)}

@app.get("/pubchem/cid/{cid}")
async def get_pubchem_by_cid(cid: int):
    """Get compound information from PubChem by CID"""
    try:
        from api_integrations import PubChemAPI
        api = PubChemAPI()
        data = api.get_compound_by_cid(cid)
        
        if data:
            return {
                "success": True,
                "compound": data,
                "source": "PubChem",
                "cid": cid,
                "timestamp": datetime.now().isoformat()
            }
        return {"success": False, "error": f"CID {cid} not found in PubChem"}
    except Exception as e:
        logger.error(f"PubChem API error: {str(e)}")
        return {"success": False, "error": str(e)}

@app.get("/pubchem/properties/{compound_name}")
async def get_pubchem_properties(compound_name: str):
    """Get molecular properties from PubChem"""
    try:
        from api_integrations import PubChemAPI
        api = PubChemAPI()
        properties = api.get_compound_properties(compound_name)
        
        if properties:
            return {
                "success": True,
                "properties": properties,
                "compound": compound_name,
                "timestamp": datetime.now().isoformat()
            }
        return {"success": False, "error": "Properties not found"}
    except Exception as e:
        logger.error(f"PubChem properties error: {str(e)}")
        return {"success": False, "error": str(e)}

# Research Engine Endpoints
@app.post("/research/run-cycle")
async def run_research_cycle(request: ResearchGoalsList):
    """Run a daily research cycle with specified goals"""
    try:
        logger.info(f"Starting research cycle with {len(request.goals)} goals")
        results = research_engine.run_daily_cycle(request.goals)
        
        return {
            "success": True,
            "summary": results,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Research cycle error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/research/session-summary")
async def get_session_summary():
    """Get summary of the last research session"""
    try:
        summary = research_engine.get_session_summary()
        return {
            "success": True,
            "summary": summary,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Session summary error: {str(e)}")
        return {"success": False, "error": str(e)}

# Article Database Endpoints
@app.get("/articles/all")
async def get_all_articles():
    """Get all research articles from the database"""
    try:
        articles = article_db.get_all_articles()
        return {
            "success": True,
            "articles": articles,
            "count": len(articles),
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Get articles error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/articles/search")
async def search_articles(request: ArticleSearchRequest):
    """Search articles by keyword"""
    try:
        articles = article_db.search_articles(request.keyword, limit=request.limit)
        return {
            "success": True,
            "articles": articles,
            "count": len(articles),
            "keyword": request.keyword,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Search articles error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/articles/statistics")
async def get_article_statistics():
    """Get statistics about the article database"""
    try:
        stats = article_db.get_statistics()
        return {
            "success": True,
            "statistics": stats,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Statistics error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/articles/by-year/{year}")
async def get_articles_by_year(year: int):
    """Get articles published in a specific year"""
    try:
        articles = article_db.get_articles_by_year(year)
        return {
            "success": True,
            "articles": articles,
            "count": len(articles),
            "year": year,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Get by year error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# RDKit Integration Endpoints
@app.post("/rdkit/sync")
async def sync_with_rdkit(request: RDKitSyncRequest):
    """Synchronize research findings with RDKit analog generation"""
    try:
        logger.info(f"Starting RDKit sync with max_compounds={request.max_compounds}")
        results = rdkit_integration.sync_research_with_rdkit(max_compounds=request.max_compounds)
        
        return {
            "success": True,
            "results": results,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"RDKit sync error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/rdkit/sync-summary")
async def get_rdkit_sync_summary():
    """Get summary of the last RDKit synchronization"""
    try:
        summary = rdkit_integration.get_sync_summary()
        return {
            "success": True,
            "summary": summary,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"RDKit summary error: {str(e)}")
        return {"success": False, "error": str(e)}

# Analog Discoveries Endpoints
@app.get("/discoveries/all")
async def get_all_discoveries():
    """Get all analog discoveries"""
    try:
        discoveries = rdkit_integration.get_all_discoveries()
        return {
            "success": True,
            "discoveries": discoveries,
            "total_analogs": sum(len(d.get("analogs", [])) for d in discoveries),
            "total_sessions": len(discoveries),
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Get discoveries error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/discoveries/statistics")
async def get_discovery_statistics():
    """Get statistics about analog discoveries"""
    try:
        stats = rdkit_integration.get_discovery_statistics()
        return {
            "success": True,
            "statistics": stats,
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Discovery statistics error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# Combined Research + PubChem Endpoint
@app.post("/research/enhanced-cycle")
async def run_enhanced_research_cycle(request: ResearchGoalsList):
    """Run enhanced research cycle with PubChem validation"""
    try:
        logger.info(f"Starting enhanced research cycle with PubChem integration")
        
        # Run standard research cycle
        research_results = research_engine.run_daily_cycle(request.goals)
        
        # Validate compounds with PubChem
        from api_integrations import PubChemAPI
        pubchem = PubChemAPI()
        
        validated_compounds = []
        for goal in request.goals:
            # Extract compound name from goal (simple extraction)
            words = goal.split()
            if words:
                compound_name = words[0]
                pubchem_data = pubchem.get_compound_by_name(compound_name)
                if pubchem_data:
                    validated_compounds.append({
                        "name": compound_name,
                        "pubchem_data": pubchem_data,
                        "validated": True
                    })
        
        return {
            "success": True,
            "research_results": research_results,
            "pubchem_validation": {
                "validated_compounds": validated_compounds,
                "count": len(validated_compounds)
            },
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Enhanced research cycle error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8006)

