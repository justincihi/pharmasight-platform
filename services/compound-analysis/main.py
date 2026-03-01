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


# ============================================================================
# ADDITIONAL ENDPOINTS FOR DASHBOARD INTEGRATION
# ============================================================================

class ADMETRequest(BaseModel):
    smiles: str

class ToxicityRequest(BaseModel):
    smiles: str

class DockingRequest(BaseModel):
    ligand_smiles: str
    receptor_pdb: str


@app.post("/admet/predict")
async def predict_admet(request: ADMETRequest):
    """Predict ADMET properties for a compound."""
    try:
        mol = Chem.MolFromSmiles(request.smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        # Basic ADMET prediction using RDKit descriptors
        admet_prediction = {
            "absorption": {
                "caco2_permeability": "high" if Descriptors.MolLogP(mol) > 0 else "low",
                "intestinal_absorption": "good" if Descriptors.TPSA(mol) < 140 else "poor",
                "bioavailability_score": min(100, max(0, 100 - (Descriptors.TPSA(mol) / 2)))
            },
            "distribution": {
                "logp": round(Descriptors.MolLogP(mol), 2),
                "tpsa": round(Descriptors.TPSA(mol), 2),
                "bbb_penetration": "likely" if Descriptors.TPSA(mol) < 90 and Descriptors.MolLogP(mol) > 0 else "unlikely"
            },
            "metabolism": {
                "cyp450_substrate": "potential" if Descriptors.MolWt(mol) < 500 else "unlikely",
                "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol)
            },
            "excretion": {
                "renal_clearance": "likely" if Descriptors.MolWt(mol) < 300 else "hepatic"
            },
            "toxicity": {
                "ames_toxicity": "non-mutagenic",  # Placeholder
                "hepatotoxicity": "low risk",  # Placeholder
                "skin_sensitization": "non-sensitizer"  # Placeholder
            }
        }

        return {"success": True, "data": admet_prediction}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/toxicity/predict")
async def predict_toxicity(request: ToxicityRequest):
    """Predict toxicity for a compound."""
    try:
        mol = Chem.MolFromSmiles(request.smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        # Basic toxicity prediction
        toxicity_prediction = {
            "ames_mutagenicity": "non-mutagenic",  # Placeholder
            "carcinogenicity": "non-carcinogenic",  # Placeholder
            "acute_toxicity": {
                "ld50_oral_rat": "low toxicity",
                "classification": "Category 4"
            },
            "hepatotoxicity": "low risk" if Descriptors.MolLogP(mol) < 5 else "moderate risk",
            "cardiotoxicity": "low risk",
            "cytotoxicity": "low",
            "alerts": []
        }

        # Add alerts for concerning features
        if Descriptors.NumAromaticRings(mol) > 3:
            toxicity_prediction["alerts"].append("Multiple aromatic rings - check for mutagenicity")

        if Descriptors.MolWt(mol) > 500:
            toxicity_prediction["alerts"].append("High molecular weight - potential absorption issues")

        return {"success": True, "data": toxicity_prediction}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/docking/simulate")
async def simulate_docking(request: DockingRequest):
    """Simulate molecular docking."""
    try:
        mol = Chem.MolFromSmiles(request.ligand_smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        # Placeholder docking results
        # In production, this would call AutoDock Vina
        docking_results = {
            "receptor": request.receptor_pdb,
            "ligand": request.ligand_smiles,
            "binding_affinity": round(-7.5 + (Descriptors.MolWt(mol) / 100), 2),  # Simulated
            "best_pose": {
                "score": -7.5,
                "rmsd": 0.0
            },
            "poses": [
                {"score": -7.5, "rmsd": 0.0},
                {"score": -7.2, "rmsd": 1.8},
                {"score": -6.9, "rmsd": 2.1}
            ],
            "notes": "Docking simulation - production version would use AutoDock Vina"
        }

        return {"success": True, "data": docking_results}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
