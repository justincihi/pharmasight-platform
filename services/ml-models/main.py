from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any
import pandas as pd
import numpy as np

# Try to import RDKit for feature generation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

app = FastAPI(title="PharmaSight ML Models Service", version="1.0.0")

# Load pre-trained models (paths would be to mounted volumes in a real setup)
admet_model = None # joblib.load('models/admet_predictor.pkl')
property_model = None # joblib.load('models/property_predictor.pkl')

class ADMETRequest(BaseModel):
    smiles_list: List[str]

class PropertyRequest(BaseModel):
    smiles: str

def featurize_smiles(smiles_list: List[str]) -> pd.DataFrame:
    """Generate molecular descriptors from SMILES strings using RDKit."""
    if not RDKIT_AVAILABLE:
        # Return dummy features if RDKit is not available
        return pd.DataFrame({'smiles': smiles_list, 'feature1': range(len(smiles_list))})
    
    features = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Invalid SMILES, use zeros
            features.append({
                'molecular_weight': 0,
                'logp': 0,
                'tpsa': 0,
                'hbd': 0,
                'hba': 0,
                'rotatable_bonds': 0
            })
        else:
            features.append({
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
            })
    
    return pd.DataFrame(features)

def predict_admet_heuristic(smiles: str) -> Dict[str, Any]:
    """Heuristic ADMET prediction based on molecular properties."""
    if not RDKIT_AVAILABLE:
        return {"toxicity": 0.1, "absorption": 0.9, "warning": "RDKit not available"}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    # Simple heuristic rules
    lipinski_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])
    
    # Estimate absorption (high if Lipinski compliant and good TPSA)
    absorption = max(0, min(1, (1 - lipinski_violations/4) * (1 - abs(tpsa - 90)/200)))
    
    # Estimate toxicity (higher with more violations and extreme logP)
    toxicity = min(1, (lipinski_violations/4 + abs(logp - 2)/10) / 2)
    
    return {
        "absorption": round(absorption, 3),
        "toxicity": round(toxicity, 3),
        "lipinski_violations": lipinski_violations,
        "bioavailability_score": round((1 - toxicity) * absorption, 3)
    }

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "service": "ml-models",
        "rdkit_available": RDKIT_AVAILABLE,
        "models_loaded": {
            "admet": admet_model is not None,
            "property": property_model is not None
        }
    }

@app.post("/predict/admet", response_model=Dict[str, Any])
async def predict_admet(request: ADMETRequest):
    """Predicts ADMET properties for a list of SMILES strings."""
    if not admet_model:
        # Use heuristic predictions if no model is loaded
        predictions = {}
        for smiles in request.smiles_list:
            predictions[smiles] = predict_admet_heuristic(smiles)
        
        return {
            "warning": "ADMET model not loaded. Using heuristic predictions.",
            "predictions": predictions
        }
    
    features = featurize_smiles(request.smiles_list)
    predictions = admet_model.predict(features)
    
    results = {smiles: pred for smiles, pred in zip(request.smiles_list, predictions)}
    return {"predictions": results}

@app.post("/predict/properties", response_model=Dict[str, Any])
async def predict_properties(request: PropertyRequest):
    """Predicts molecular properties for a single SMILES string."""
    if not property_model:
        # Use RDKit descriptors if no model is loaded
        if not RDKIT_AVAILABLE:
            return {
                "warning": "Property model not loaded and RDKit not available.",
                "logp": 2.5,
                "solubility": -3.0
            }
        
        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        return {
            "warning": "Property model not loaded. Using RDKit descriptors.",
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol)
        }
        
    features = featurize_smiles([request.smiles])
    prediction = property_model.predict(features)[0]
    
    return {"prediction": prediction}
