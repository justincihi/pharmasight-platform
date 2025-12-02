from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import Optional, Dict, Any, List
import redis
import pickle
from functools import wraps
import sys
import os

# RDKit imports
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, QED, AllChem
from rdkit import RDConfig

# Add SA_Score to path
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
try:
    import sascorer
    SA_SCORE_AVAILABLE = True
except ImportError:
    SA_SCORE_AVAILABLE = False
    print("Warning: SA_Score not available")

# Add NP_Score to path
sys.path.append(os.path.join(RDConfig.RDContribDir, 'NP_Score'))
try:
    import npscorer
    NP_SCORE_AVAILABLE = True
except ImportError:
    NP_SCORE_AVAILABLE = False
    print("Warning: NP_Score not available")

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
    include_viability: Optional[bool] = True

class ViabilityAnalysisRequest(BaseModel):
    smiles: str
    include_3d: Optional[bool] = False

# Compound database
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

def calculate_synthetic_accessibility(mol) -> Dict[str, Any]:
    """Calculate synthetic accessibility score (1=easy, 10=difficult)."""
    if not SA_SCORE_AVAILABLE:
        return {
            "sa_score": None,
            "sa_interpretation": "SA Score module not available",
            "synthetic_feasibility": "Unknown"
        }
    
    try:
        sa_score = sascorer.calculateScore(mol)
        
        # Interpret the score
        if sa_score <= 3:
            feasibility = "Very Easy"
        elif sa_score <= 5:
            feasibility = "Easy"
        elif sa_score <= 7:
            feasibility = "Moderate"
        elif sa_score <= 9:
            feasibility = "Difficult"
        else:
            feasibility = "Very Difficult"
        
        return {
            "sa_score": round(sa_score, 2),
            "sa_interpretation": f"Score of {sa_score:.2f} on scale of 1-10",
            "synthetic_feasibility": feasibility
        }
    except Exception as e:
        return {
            "sa_score": None,
            "sa_interpretation": f"Error: {str(e)}",
            "synthetic_feasibility": "Unknown"
        }

def calculate_drug_likeness(mol) -> Dict[str, Any]:
    """Calculate QED (Quantitative Estimate of Drug-likeness)."""
    try:
        qed_score = QED.qed(mol)
        qed_properties = QED.properties(mol)
        
        # Interpret QED score
        if qed_score >= 0.7:
            interpretation = "Excellent drug-like properties"
        elif qed_score >= 0.5:
            interpretation = "Good drug-like properties"
        elif qed_score >= 0.3:
            interpretation = "Moderate drug-like properties"
        else:
            interpretation = "Poor drug-like properties"
        
        return {
            "qed_score": round(qed_score, 3),
            "qed_interpretation": interpretation,
            "qed_properties": {
                "MW": round(qed_properties.MW, 2),
                "ALOGP": round(qed_properties.ALOGP, 2),
                "HBA": qed_properties.HBA,
                "HBD": qed_properties.HBD,
                "PSA": round(qed_properties.PSA, 2),
                "ROTB": qed_properties.ROTB,
                "AROM": qed_properties.AROM,
                "ALERTS": qed_properties.ALERTS
            }
        }
    except Exception as e:
        return {
            "qed_score": None,
            "qed_interpretation": f"Error: {str(e)}",
            "qed_properties": {}
        }

def calculate_natural_product_likeness(mol) -> Dict[str, Any]:
    """Calculate natural product-likeness score."""
    if not NP_SCORE_AVAILABLE:
        return {
            "np_score": None,
            "np_interpretation": "NP Score module not available",
            "natural_product_like": "Unknown"
        }
    
    try:
        np_score = npscorer.scoreMol(mol)
        
        # Interpret the score (typically ranges from -5 to +5)
        if np_score >= 1:
            interpretation = "Highly natural product-like"
        elif np_score >= 0:
            interpretation = "Moderately natural product-like"
        elif np_score >= -1:
            interpretation = "Slightly natural product-like"
        else:
            interpretation = "Not natural product-like"
        
        return {
            "np_score": round(np_score, 2),
            "np_interpretation": interpretation,
            "natural_product_like": "Yes" if np_score >= 0 else "No"
        }
    except Exception as e:
        return {
            "np_score": None,
            "np_interpretation": f"Error: {str(e)}",
            "natural_product_like": "Unknown"
        }

def calculate_lipinski_rule_of_five(mol) -> Dict[str, Any]:
    """Calculate Lipinski's Rule of Five parameters."""
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    violations = 0
    rules = {}
    
    rules["molecular_weight"] = {
        "value": round(mw, 2),
        "limit": "≤ 500 Da",
        "passes": mw <= 500
    }
    if mw > 500:
        violations += 1
    
    rules["logp"] = {
        "value": round(logp, 2),
        "limit": "≤ 5",
        "passes": logp <= 5
    }
    if logp > 5:
        violations += 1
    
    rules["h_bond_donors"] = {
        "value": hbd,
        "limit": "≤ 5",
        "passes": hbd <= 5
    }
    if hbd > 5:
        violations += 1
    
    rules["h_bond_acceptors"] = {
        "value": hba,
        "limit": "≤ 10",
        "passes": hba <= 10
    }
    if hba > 10:
        violations += 1
    
    return {
        "violations": violations,
        "passes_ro5": violations <= 1,  # Allow 1 violation
        "rules": rules,
        "interpretation": f"{violations} violation(s) - " + 
                         ("Good oral bioavailability expected" if violations <= 1 else "Poor oral bioavailability likely")
    }

def calculate_viability_score(mol) -> Dict[str, Any]:
    """Calculate overall viability score (0-100)."""
    scores = {}
    weights = {}
    
    # Synthetic Accessibility (20% weight)
    sa_data = calculate_synthetic_accessibility(mol)
    if sa_data["sa_score"] is not None:
        # Invert SA score (lower is better) and normalize to 0-100
        sa_normalized = max(0, 100 - (sa_data["sa_score"] * 10))
        scores["synthetic_accessibility"] = sa_normalized
        weights["synthetic_accessibility"] = 0.20
    
    # Drug-likeness QED (15% weight)
    qed_data = calculate_drug_likeness(mol)
    if qed_data["qed_score"] is not None:
        scores["drug_likeness"] = qed_data["qed_score"] * 100
        weights["drug_likeness"] = 0.15
    
    # Lipinski's Rule of Five (15% weight)
    ro5_data = calculate_lipinski_rule_of_five(mol)
    ro5_score = 100 if ro5_data["passes_ro5"] else max(0, 100 - (ro5_data["violations"] * 25))
    scores["lipinski_compliance"] = ro5_score
    weights["lipinski_compliance"] = 0.15
    
    # Placeholder scores for future integration (50% weight total)
    # These will be replaced with actual predictions
    scores["safety_profile"] = 75  # 25% weight - placeholder
    weights["safety_profile"] = 0.25
    
    scores["efficacy_prediction"] = 70  # 20% weight - placeholder
    weights["efficacy_prediction"] = 0.20
    
    scores["patent_freedom"] = 80  # 20% weight - placeholder (will be calculated from similarity search)
    weights["patent_freedom"] = 0.05  # Reduced weight for now
    
    # Calculate weighted average
    total_weight = sum(weights.values())
    if total_weight > 0:
        viability_score = sum(scores[k] * weights[k] for k in scores) / total_weight
    else:
        viability_score = 0
    
    # Determine priority level
    if viability_score >= 80:
        priority = "Very High"
        recommendation = "Excellent candidate - proceed to synthesis"
    elif viability_score >= 65:
        priority = "High"
        recommendation = "Good candidate - further analysis recommended"
    elif viability_score >= 50:
        priority = "Medium"
        recommendation = "Moderate candidate - consider optimization"
    else:
        priority = "Low"
        recommendation = "Poor candidate - significant improvements needed"
    
    return {
        "overall_score": round(viability_score, 1),
        "priority": priority,
        "recommendation": recommendation,
        "component_scores": {k: round(v, 1) for k, v in scores.items()},
        "weights": weights
    }

@cache_result(expiration=1800)
def perform_rdkit_analysis(smiles: str) -> Dict[str, Any]:
    """Perform molecular analysis using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES string"}

        return {
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "num_h_donors": Descriptors.NumHDonors(mol),
            "num_h_acceptors": Descriptors.NumHAcceptors(mol),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
            "num_saturated_rings": Descriptors.NumSaturatedRings(mol),
            "fraction_csp3": round(Descriptors.FractionCSP3(mol), 3),
        }
    except Exception as e:
        return {"error": f"RDKit analysis failed: {str(e)}"}

@cache_result(expiration=3600)
def perform_viability_analysis(smiles: str, include_3d: bool = False) -> Dict[str, Any]:
    """Perform comprehensive viability analysis."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES string"}
        
        # Add hydrogens if 3D analysis requested
        if include_3d:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        
        analysis = {
            "synthetic_accessibility": calculate_synthetic_accessibility(mol),
            "drug_likeness": calculate_drug_likeness(mol),
            "natural_product_likeness": calculate_natural_product_likeness(mol),
            "lipinski_rule_of_five": calculate_lipinski_rule_of_five(mol),
            "viability_score": calculate_viability_score(mol)
        }
        
        return analysis
    except Exception as e:
        return {"error": f"Viability analysis failed: {str(e)}"}

def generate_svg_structure(smiles: str) -> str:
    """Generate SVG representation of chemical structure using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "<svg>Invalid SMILES</svg>"
    
    Draw.PrepareMolForDrawing(mol)
    drawer = Draw.MolDraw2DSVG(300, 200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    return drawer.GetDrawingText()

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    try:
        redis_client.ping()
        redis_status = "connected"
    except Exception as e:
        redis_status = f"error: {str(e)}"
    
    return {
        "status": "healthy",
        "service": "compound-analysis-enhanced",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "redis": redis_status,
        "features": {
            "sa_score": SA_SCORE_AVAILABLE,
            "np_score": NP_SCORE_AVAILABLE,
            "qed": True,
            "viability_scoring": True
        }
    }

@app.post("/analyze", response_model=Dict[str, Any])
async def analyze_compound(request: CompoundAnalysisRequest):
    """Analyzes a chemical compound and returns its properties."""
    compound_data = get_compound_data(request.compound)
    
    if not compound_data:
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
    
    # Add viability analysis if requested
    if request.include_viability:
        viability_results = perform_viability_analysis(smiles)
        if "error" not in viability_results:
            response["viability_analysis"] = viability_results
    
    return response

@app.post("/viability", response_model=Dict[str, Any])
async def analyze_viability(request: ViabilityAnalysisRequest):
    """Perform comprehensive viability analysis on a compound."""
    viability_results = perform_viability_analysis(request.smiles, request.include_3d)
    
    if "error" in viability_results:
        raise HTTPException(status_code=400, detail=viability_results["error"])
    
    # Add basic RDKit analysis
    basic_analysis = perform_rdkit_analysis(request.smiles)
    
    return {
        "smiles": request.smiles,
        "basic_properties": basic_analysis,
        "viability_analysis": viability_results,
        "svg_image": generate_svg_structure(request.smiles)
    }

@app.post("/batch_viability", response_model=List[Dict[str, Any]])
async def batch_analyze_viability(smiles_list: List[str]):
    """Perform viability analysis on multiple compounds."""
    results = []
    
    for smiles in smiles_list:
        try:
            viability_results = perform_viability_analysis(smiles)
            if "error" not in viability_results:
                results.append({
                    "smiles": smiles,
                    "viability_score": viability_results["viability_score"]["overall_score"],
                    "priority": viability_results["viability_score"]["priority"],
                    "full_analysis": viability_results
                })
            else:
                results.append({
                    "smiles": smiles,
                    "error": viability_results["error"]
                })
        except Exception as e:
            results.append({
                "smiles": smiles,
                "error": str(e)
            })
    
    # Sort by viability score (highest first)
    results.sort(key=lambda x: x.get("viability_score", 0), reverse=True)
    
    return results
