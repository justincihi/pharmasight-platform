from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict
from rdkit import Chem
from rdkit_analog_generator import RDKitAnalogGenerator
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="Analog Generation Service")

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize RDKit generator
generator = RDKitAnalogGenerator()

class AnalogRequest(BaseModel):
    parent_smiles: str
    parent_name: Optional[str] = "Unknown"
    num_analogs: Optional[int] = 10
    target_properties: Optional[str] = "all"

class AnalogResponse(BaseModel):
    parent_compound: str
    parent_smiles: str
    analogs: List[Dict]
    count: int
    generation_method: str

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "service": "analog-generation",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "generator_ready": True
    }

@app.post("/generate", response_model=AnalogResponse)
async def generate_analogs(request: AnalogRequest):
    """
    Generate novel molecular analogs using RDKit transformations.
    
    Args:
        parent_smiles: SMILES string of parent compound
        parent_name: Name of parent compound (optional)
        num_analogs: Number of analogs to generate (default: 10)
        target_properties: Filter criteria (all, patent-free, high-similarity)
    
    Returns:
        List of generated analogs with properties and predictions
    """
    try:
        # Validate SMILES
        mol = Chem.MolFromSmiles(request.parent_smiles)
        if mol is None:
            raise HTTPException(
                status_code=400, 
                detail=f"Invalid SMILES string: {request.parent_smiles}"
            )
        
        logger.info(f"Generating {request.num_analogs} analogs for {request.parent_name}")
        
        # Generate analogs using RDKit
        analogs = generator.generate_analogs(
            parent_smiles=request.parent_smiles,
            parent_name=request.parent_name,
            num_analogs=request.num_analogs
        )
        
        # Filter based on target properties
        filtered_analogs = []
        for analog in analogs:
            passes_filter = False
            
            if request.target_properties == "all":
                passes_filter = True
            elif request.target_properties == "patent-free" and "Patent-Free" in analog.get("patent_status", ""):
                passes_filter = True
            elif request.target_properties == "high-similarity" and analog.get("similarity", 0) > 0.8:
                passes_filter = True
            elif request.target_properties == "high-value" and analog.get("patent_opportunity_score", 0) >= 90:
                passes_filter = True
            
            if passes_filter:
                filtered_analogs.append(analog)
        
        logger.info(f"Generated {len(filtered_analogs)} analogs (filtered from {len(analogs)})")
        
        return {
            "parent_compound": request.parent_name,
            "parent_smiles": request.parent_smiles,
            "analogs": filtered_analogs,
            "count": len(filtered_analogs),
            "generation_method": "RDKit Structural Transformation"
        }
        
    except Exception as e:
        logger.error(f"Error generating analogs: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/transformations")
async def list_transformations():
    """List available transformation strategies."""
    return {
        "strategies": [
            "Add methyl group",
            "Add fluorine",
            "Add hydroxyl",
            "Replace hydrogen with halogen",
            "Extend carbon chain",
            "Add methoxy group",
            "Ring expansion",
            "Add amino group"
        ],
        "count": 8
    }

@app.post("/batch-generate")
async def batch_generate_analogs(compounds: List[AnalogRequest]):
    """
    Generate analogs for multiple parent compounds in batch.
    
    Args:
        compounds: List of analog generation requests
    
    Returns:
        List of analog generation results
    """
    results = []
    
    for request in compounds:
        try:
            result = await generate_analogs(request)
            results.append({
                "parent": request.parent_name,
                "status": "success",
                "data": result
            })
        except Exception as e:
            results.append({
                "parent": request.parent_name,
                "status": "error",
                "error": str(e)
            })
    
    return {
        "results": results,
        "total": len(results),
        "successful": len([r for r in results if r["status"] == "success"]),
        "failed": len([r for r in results if r["status"] == "error"])
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8002)

