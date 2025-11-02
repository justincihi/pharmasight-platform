from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import Optional, Dict, Any
import redis
import pickle
from functools import wraps

# RDKit imports will be available in the container
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

app = FastAPI()

# Connect to Redis
redis_client = redis.Redis(host='redis', port=6379, db=0)

def cache_result(expiration=3600):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            cache_key = f"{func.__name__}:{pickle.dumps((args, kwargs))}"
            cached = redis_client.get(cache_key)
            if cached:
                return pickle.loads(cached)
            
            result = func(*args, **kwargs)
            redis_client.setex(cache_key, expiration, pickle.dumps(result))
            return result
        return wrapper
    return decorator

class CompoundAnalysisRequest(BaseModel):
    compound: str
    analysis_type: Optional[str] = "full"

# This would be populated from a database in a real application
COMPOUND_DATABASE = {
    "psilocybin": {
        "smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "name": "Psilocybin",
        "molecular_weight": 284.25,
        "therapeutic_area": "Psychedelic Therapy",
        "status": "Phase II Clinical Trials",
    },
    "lsd": {
        "smiles": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C",
        "name": "LSD",
        "molecular_weight": 323.43,
        "therapeutic_area": "Psychedelic Research",
        "status": "Research Phase",
    },
    "mdma": {
        "smiles": "CC(CC1=CC2=C(C=C1)OCO2)NC",
        "name": "MDMA",
        "molecular_weight": 193.25,
        "therapeutic_area": "PTSD Therapy",
        "status": "Phase III Clinical Trials",
    },
    # ... Add more compounds from the original database
}

def get_compound_data(identifier: str) -> Optional[Dict[str, Any]]:
    """Retrieve compound data by name or SMILES."""
    identifier_lower = identifier.lower()
    if identifier_lower in COMPOUND_DATABASE:
        return COMPOUND_DATABASE[identifier_lower]
    for data in COMPOUND_DATABASE.values():
        if data["smiles"] == identifier:
            return data
    return None

@cache_result(expiration=1800)
def perform_rdkit_analysis(smiles: str) -> Dict[str, Any]:
    """Perform molecular analysis using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES string"}

        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "num_h_donors": Descriptors.NumHDonors(mol),
            "num_h_acceptors": Descriptors.NumHAcceptors(mol),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
        }
    except Exception as e:
        return {"error": f"RDKit analysis failed: {str(e)}"}

def generate_svg_structure(smiles: str) -> str:
    """Generate SVG representation of chemical structure using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "<svg>Invalid SMILES</svg>"
    
    # Generate a 2D depiction of the molecule
    Draw.PrepareMolForDrawing(mol)
    
    # Use the rdkit MolDraw2DSVG drawer
    drawer = Draw.MolDraw2DSVG(300, 200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    return drawer.GetDrawingText()

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    try:
        # Test Redis connection
        redis_client.ping()
        redis_status = "connected"
    except Exception as e:
        redis_status = f"error: {str(e)}"
    
    return {
        "status": "healthy",
        "service": "compound-analysis",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "redis": redis_status
    }

@app.post("/analyze", response_model=Dict[str, Any])
async def analyze_compound(request: CompoundAnalysisRequest):
    """Analyzes a chemical compound and returns its properties."""
    compound_data = get_compound_data(request.compound)
    
    if not compound_data:
        # If not in DB, assume it's a SMILES string
        smiles = request.compound
    else:
        smiles = compound_data.get("smiles")

    if not smiles:
        raise HTTPException(status_code=404, detail="Compound not found and not a valid SMILES string.")

    analysis_results = perform_rdkit_analysis(smiles)
    if "error" in analysis_results:
        raise HTTPException(status_code=400, detail=analysis_results["error"])

    svg_image = generate_svg_structure(smiles)

    response = {
        "base_properties": compound_data or {"name": "Custom Compound", "smiles": smiles},
        "rdkit_analysis": analysis_results,
        "svg_image": svg_image,
    }
    
    return response
