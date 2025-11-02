from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Optional, Dict
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

app = FastAPI()

# This would be populated from a database in a real application
COMPOUND_DATABASE = {
    "psilocybin": {"smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12", "name": "Psilocybin"},
    "ketamine": {"smiles": "CNC1(CCCCC1=O)c2ccccc2Cl", "name": "Ketamine"},
    "mdma": {"smiles": "CC(CC1=CC2=C(C=C1)OCO2)NC", "name": "MDMA"},
    # ... Add more compounds
}

ANALOG_GENERATION_DATABASE = {
    "psilocybin": [
        {"name": "4-HO-MET", "smiles": "CCN(C)CCc1c[nH]c2ccc(O)cc12", "similarity": 0.92, "patent_status": "Expired"},
        {"name": "4-AcO-DMT", "smiles": "CC(=O)Oc1ccc2c(c1)c(CCN(C)C)c[nH]2", "similarity": 0.88, "patent_status": "Expired"},
        {"name": "CYP-001", "smiles": "CN(C)CCc1c[nH]c2ccc(OC(F)(F)F)cc12", "similarity": 0.85, "patent_status": "Patented (US11,234,567)"},
    ],
    # ... Add more analogs
}

class AnalogRequest(BaseModel):
    parent_compound: str
    target_properties: Optional[str] = "all"

def get_smiles_from_identifier(identifier: str) -> Optional[str]:
    """Get SMILES from compound name or return identifier if it's already a SMILES."""
    if identifier.lower() in COMPOUND_DATABASE:
        return COMPOUND_DATABASE[identifier.lower()]["smiles"]
    # Basic check if it's a SMILES string
    if Chem.MolFromSmiles(identifier):
        return identifier
    return None

def calculate_similarity(smiles1: str, smiles2: str) -> float:
    """Calculate Tanimoto similarity between two SMILES strings."""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if not mol1 or not mol2:
        return 0.0
    
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    
    return DataStructs.TanimotoSimilarity(fp1, fp2)

@app.post("/generate", response_model=Dict)
async def generate_analogs(request: AnalogRequest):
    """Generates analogs for a given parent compound."""
    parent_smiles = get_smiles_from_identifier(request.parent_compound)
    if not parent_smiles:
        raise HTTPException(status_code=404, detail="Parent compound not found or invalid SMILES.")

    analogs = ANALOG_GENERATION_DATABASE.get(request.parent_compound.lower(), [])

    # In a real application, you would generate analogs on the fly
    # For now, we just filter the pre-defined list
    
    filtered_analogs = []
    for analog in analogs:
        passes_filter = False
        if request.target_properties == "all":
            passes_filter = True
        elif request.target_properties == "patent-free" and "Expired" in analog["patent_status"]:
            passes_filter = True
        elif request.target_properties == "high-similarity" and analog["similarity"] > 0.9:
            passes_filter = True
        
        if passes_filter:
            # Recalculate similarity for demonstration
            analog["similarity"] = calculate_similarity(parent_smiles, analog["smiles"])
            filtered_analogs.append(analog)

    return {
        "parent_compound": request.parent_compound,
        "analogs": filtered_analogs,
        "count": len(filtered_analogs)
    }
