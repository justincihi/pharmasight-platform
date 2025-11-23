from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware
import os
import json
from datetime import datetime

app = FastAPI(title="PharmaSight Frontend", version="1.0.0")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")

# Templates
templates = Jinja2Templates(directory="templates")

# Mock data for demo
MOCK_ARTICLES = [
    {
        "title": "Novel Ketamine Analogs for Treatment-Resistant Depression",
        "authors": "Smith J, Johnson M, Williams K",
        "journal": "Journal of Psychopharmacology",
        "year": 2024,
        "doi": "10.1177/0269881124001234",
        "pmid": "38123456",
        "abstract": "This study investigates novel ketamine analogs with improved safety profiles...",
        "relevance_score": 95
    },
    {
        "title": "NMDA Receptor Modulation in Psychiatric Disorders",
        "authors": "Brown A, Davis R, Miller S",
        "journal": "Nature Neuroscience",
        "year": 2023,
        "doi": "10.1038/s41593-023-01234-5",
        "pmid": "37234567",
        "abstract": "Comprehensive review of NMDA receptor mechanisms in mental health...",
        "relevance_score": 88
    },
    {
        "title": "Rapid-Acting Antidepressants: Mechanisms and Clinical Applications",
        "authors": "Garcia M, Rodriguez L, Martinez P",
        "journal": "Molecular Psychiatry",
        "year": 2023,
        "doi": "10.1038/s41380-023-02345-6",
        "pmid": "36345678",
        "abstract": "Analysis of rapid-acting antidepressant mechanisms including ketamine derivatives...",
        "relevance_score": 92
    }
]

MOCK_DISCOVERIES = [
    {
        "analog_id": "ANALOG_001",
        "name": "4-Methoxy-Ketamine",
        "smiles": "CNC1(CCCCC1=O)c1ccc(OC)cc1Cl",
        "parent_compound": "Ketamine",
        "transformation": "Methoxy substitution at para position",
        "patent_status": "Patent-Free",
        "ip_opportunity_score": 95,
        "discovery_date": "2024-11-15",
        "predicted_properties": {
            "molecular_weight": 273.76,
            "logP": 2.8,
            "hbd": 1,
            "hba": 3
        }
    },
    {
        "analog_id": "ANALOG_002",
        "name": "3-Fluoro-Ketamine",
        "smiles": "CNC1(CCCCC1=O)c1cc(F)ccc1Cl",
        "parent_compound": "Ketamine",
        "transformation": "Fluorine substitution at meta position",
        "patent_status": "Patent-Free",
        "ip_opportunity_score": 88,
        "discovery_date": "2024-11-14",
        "predicted_properties": {
            "molecular_weight": 261.72,
            "logP": 2.9,
            "hbd": 1,
            "hba": 2
        }
    },
    {
        "analog_id": "ANALOG_003",
        "name": "N-Ethyl-Ketamine",
        "smiles": "CCNC1(CCCCC1=O)c1ccccc1Cl",
        "parent_compound": "Ketamine",
        "transformation": "N-methylation to N-ethylation",
        "patent_status": "Patent-Free",
        "ip_opportunity_score": 82,
        "discovery_date": "2024-11-13",
        "predicted_properties": {
            "molecular_weight": 257.76,
            "logP": 3.1,
            "hbd": 1,
            "hba": 2
        }
    }
]

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Serve the main frontend page"""
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/health")
async def health():
    """Health check endpoint"""
    return {"status": "healthy", "service": "frontend", "timestamp": datetime.now().isoformat()}

# Mock API endpoints for demo
@app.get("/api/research/articles")
async def get_articles():
    """Get all research articles (mock data)"""
    return {
        "status": "success",
        "count": len(MOCK_ARTICLES),
        "articles": MOCK_ARTICLES
    }

@app.get("/api/research/discoveries")
async def get_discoveries():
    """Get all analog discoveries (mock data)"""
    return {
        "status": "success",
        "count": len(MOCK_DISCOVERIES),
        "total_analogs": 119,
        "patent_free": 104,
        "discoveries": MOCK_DISCOVERIES
    }

@app.get("/api/stats")
async def get_stats():
    """Get platform statistics (mock data)"""
    return {
        "compounds": 523,
        "analogs": 119,
        "articles": 26,
        "receptors": 72,
        "patent_free_percentage": 87,
        "last_updated": datetime.now().isoformat()
    }

@app.post("/api/research/run-cycle")
async def run_research_cycle():
    """Simulate running a research cycle"""
    return {
        "status": "success",
        "message": "Research cycle completed",
        "articles_found": 3,
        "new_discoveries": 2,
        "timestamp": datetime.now().isoformat()
    }

@app.post("/api/compound/analyze")
async def analyze_compound(request: Request):
    """Mock compound analysis"""
    data = await request.json()
    smiles = data.get("smiles", "")
    
    return {
        "status": "success",
        "smiles": smiles,
        "properties": {
            "molecular_weight": 243.73,
            "logP": 2.7,
            "hbd": 1,
            "hba": 2,
            "tpsa": 29.1,
            "rotatable_bonds": 2
        },
        "predictions": {
            "blood_brain_barrier": "High",
            "absorption": "Good",
            "toxicity_risk": "Low"
        }
    }

@app.post("/api/analog/generate")
async def generate_analogs(request: Request):
    """Mock analog generation"""
    data = await request.json()
    smiles = data.get("smiles", "")
    count = data.get("count", 5)
    
    return {
        "status": "success",
        "parent_smiles": smiles,
        "analogs_generated": count,
        "analogs": [
            {
                "id": f"ANALOG_{i:03d}",
                "smiles": f"Modified_{smiles}_{i}",
                "similarity": 0.95 - (i * 0.05),
                "patent_free": True
            }
            for i in range(1, count + 1)
        ]
    }

@app.get("/PUBLIC_ANALOG_DISCOVERY_REGISTRY.md")
async def registry():
    """Serve the public IP registry"""
    registry_path = "../PUBLIC_ANALOG_DISCOVERY_REGISTRY.md"
    if os.path.exists(registry_path):
        return FileResponse(registry_path, media_type="text/markdown")
    return {"error": "Registry not found"}

@app.get("/docs-page", response_class=HTMLResponse)
async def docs_page(request: Request):
    """API documentation page"""
    return """
    <html>
        <head>
            <title>PharmaSight API Documentation</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; background: #0a0e27; color: #fff; }
                h1 { color: #00a8e8; }
                .endpoint { background: #1a1f3a; padding: 20px; margin: 20px 0; border-radius: 8px; }
                code { background: #0a0e27; padding: 4px 8px; border-radius: 4px; color: #00d9ff; }
            </style>
        </head>
        <body>
            <h1>PharmaSight™ API Documentation</h1>
            <div class="endpoint">
                <h2>Research API</h2>
                <p><code>GET /api/research/articles</code> - Get all research articles</p>
                <p><code>GET /api/research/discoveries</code> - Get all analog discoveries</p>
                <p><code>POST /api/research/run-cycle</code> - Run research cycle</p>
            </div>
            <div class="endpoint">
                <h2>Compound Analysis API</h2>
                <p><code>POST /api/compound/analyze</code> - Analyze compound</p>
            </div>
            <div class="endpoint">
                <h2>Analog Generation API</h2>
                <p><code>POST /api/analog/generate</code> - Generate analogs</p>
            </div>
            <div class="endpoint">
                <h2>Statistics API</h2>
                <p><code>GET /api/stats</code> - Get platform statistics</p>
            </div>
            <p><a href="/" style="color: #00a8e8;">← Back to Platform</a></p>
        </body>
    </html>
    """

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=3000)

