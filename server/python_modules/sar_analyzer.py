#!/usr/bin/env python3
"""
Structure-Activity Relationship (SAR) Analysis Module
Analyzes molecular scaffolds, R-groups, and predicts activity cliffs
"""

import json
import sys
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Draw
    from rdkit.Chem.Scaffolds import MurckoScaffold
    import numpy as np
except ImportError as e:
    print(f"Error: Required package not installed: {e}", file=sys.stderr)
    sys.exit(1)


@dataclass
class RGroup:
    """Represents an R-group substitution"""
    position: int
    smiles: str
    name: str
    frequency: int
    avg_activity: float
    activity_std: float


@dataclass
class ActivityCliff:
    """Represents an activity cliff (large activity change from small structural change)"""
    compound1_smiles: str
    compound2_smiles: str
    compound1_activity: float
    compound2_activity: float
    activity_difference: float
    structural_similarity: float
    cliff_steepness: float


@dataclass
class SARResult:
    """Complete SAR analysis result"""
    parent_smiles: str
    scaffold: str
    r_groups: Dict[int, List[RGroup]]
    activity_cliffs: List[ActivityCliff]
    pharmacophore_features: List[str]
    key_substitutions: List[Dict[str, Any]]
    execution_time: float


def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False


def extract_scaffold(smiles: str) -> Optional[str]:
    """Extract Murcko scaffold from molecule"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def identify_r_groups(smiles: str, scaffold_smiles: str) -> Dict[int, List[str]]:
    """
    Identify R-group positions and variations
    Returns dictionary of position -> list of R-group SMILES
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        scaffold = Chem.MolFromSmiles(scaffold_smiles)
        
        if mol is None or scaffold is None:
            return {}
        
        # Use BRICS to identify R-group positions
        from rdkit.Chem import BRICS
        
        brics_bonds = BRICS.BRICSDecompose(mol)
        r_groups = {}
        
        for i, bond_smiles in enumerate(brics_bonds):
            r_groups[i] = [bond_smiles]
        
        return r_groups
    except Exception:
        return {}


def identify_pharmacophore_features(mol: Chem.Mol) -> List[str]:
    """Identify key pharmacophore features"""
    features = []
    
    # Hydrogen bond donors
    hbd = Descriptors.NumHDonors(mol)
    if hbd > 0:
        features.append(f"HBD ({hbd})")
    
    # Hydrogen bond acceptors
    hba = Descriptors.NumHAcceptors(mol)
    if hba > 0:
        features.append(f"HBA ({hba})")
    
    # Aromatic rings
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    if aromatic_rings > 0:
        features.append(f"Aromatic ({aromatic_rings})")
    
    # Rotatable bonds
    rotatable = Descriptors.NumRotatableBonds(mol)
    if rotatable > 0:
        features.append(f"Rotatable ({rotatable})")
    
    # Charged atoms
    charged = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0)
    if charged > 0:
        features.append(f"Charged ({charged})")
    
    # Halogen atoms
    halogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I'])
    if halogens > 0:
        features.append(f"Halogen ({halogens})")
    
    return features


def calculate_tanimoto_similarity(smiles1: str, smiles2: str) -> float:
    """Calculate Tanimoto similarity between two molecules"""
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return 0.0
        
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        
        return float(DataStructs.TanimotoSimilarity(fp1, fp2))
    except Exception:
        return 0.0


def identify_activity_cliffs(
    compounds: List[Dict[str, Any]],
    similarity_threshold: float = 0.85,
    activity_threshold: float = 2.0
) -> List[ActivityCliff]:
    """
    Identify activity cliffs in a series of compounds
    Activity cliff: similar compounds with large activity differences
    """
    cliffs = []
    
    try:
        from rdkit.DataStructs import TanimotoSimilarity
        
        for i in range(len(compounds)):
            for j in range(i + 1, len(compounds)):
                comp1 = compounds[i]
                comp2 = compounds[j]
                
                mol1 = Chem.MolFromSmiles(comp1["smiles"])
                mol2 = Chem.MolFromSmiles(comp2["smiles"])
                
                if mol1 is None or mol2 is None:
                    continue
                
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
                
                similarity = TanimotoSimilarity(fp1, fp2)
                
                # Check if compounds are similar but have different activities
                if similarity >= similarity_threshold:
                    activity_diff = abs(comp1.get("activity", 0) - comp2.get("activity", 0))
                    
                    if activity_diff >= activity_threshold:
                        cliff_steepness = activity_diff / (1 - similarity + 0.001)
                        
                        cliff = ActivityCliff(
                            compound1_smiles=comp1["smiles"],
                            compound2_smiles=comp2["smiles"],
                            compound1_activity=comp1.get("activity", 0),
                            compound2_activity=comp2.get("activity", 0),
                            activity_difference=activity_diff,
                            structural_similarity=similarity,
                            cliff_steepness=cliff_steepness
                        )
                        cliffs.append(cliff)
    except Exception as e:
        print(f"Error identifying activity cliffs: {e}", file=sys.stderr)
    
    return cliffs


def predict_r_group_improvements(
    r_groups: Dict[int, List[RGroup]],
    current_activity: float
) -> List[Dict[str, Any]]:
    """
    Predict which R-group substitutions might improve activity
    """
    improvements = []
    
    for position, groups in r_groups.items():
        for group in groups:
            # Predict improvement based on activity statistics
            if group.avg_activity > current_activity:
                predicted_improvement = group.avg_activity - current_activity
                confidence = 1.0 - (group.activity_std / (group.avg_activity + 0.001))
                
                improvements.append({
                    "position": position,
                    "r_group": group.smiles,
                    "predicted_improvement": round(predicted_improvement, 2),
                    "confidence": round(max(0, min(1, confidence)), 2),
                    "frequency": group.frequency,
                    "avg_activity": round(group.avg_activity, 2)
                })
    
    # Sort by predicted improvement
    improvements.sort(key=lambda x: x["predicted_improvement"], reverse=True)
    
    return improvements


def analyze_sar(
    parent_smiles: str,
    compound_series: Optional[List[Dict[str, Any]]] = None
) -> Optional[SARResult]:
    """
    Perform comprehensive SAR analysis on a compound
    
    Args:
        parent_smiles: SMILES string of parent compound
        compound_series: Optional list of related compounds with activities
    
    Returns:
        SARResult object with complete analysis
    """
    import time
    start_time = time.time()
    
    # Validate SMILES
    if not validate_smiles(parent_smiles):
        print(f"Error: Invalid SMILES string: {parent_smiles}", file=sys.stderr)
        return None
    
    try:
        mol = Chem.MolFromSmiles(parent_smiles)
        
        # Extract scaffold
        scaffold = extract_scaffold(parent_smiles)
        
        # Identify R-groups
        r_groups_dict = {}
        if scaffold:
            r_groups_dict = identify_r_groups(parent_smiles, scaffold)
        
        # Identify pharmacophore features
        pharmacophore = identify_pharmacophore_features(mol)
        
        # Identify activity cliffs if compound series provided
        cliffs = []
        key_substitutions = []
        if compound_series:
            cliffs = identify_activity_cliffs(compound_series)
            # Predict R-group improvements
            key_substitutions = predict_r_group_improvements(r_groups_dict, 0.5)
        
        execution_time = time.time() - start_time
        
        result = SARResult(
            parent_smiles=parent_smiles,
            scaffold=scaffold or parent_smiles,
            r_groups=r_groups_dict,
            activity_cliffs=cliffs,
            pharmacophore_features=pharmacophore,
            key_substitutions=key_substitutions,
            execution_time=execution_time
        )
        
        return result
    except Exception as e:
        print(f"Error analyzing SAR: {e}", file=sys.stderr)
        return None


def cliff_to_dict(cliff: ActivityCliff) -> Dict[str, Any]:
    """Convert ActivityCliff to dictionary"""
    return {
        "compound1_smiles": cliff.compound1_smiles,
        "compound2_smiles": cliff.compound2_smiles,
        "compound1_activity": round(cliff.compound1_activity, 2),
        "compound2_activity": round(cliff.compound2_activity, 2),
        "activity_difference": round(cliff.activity_difference, 2),
        "structural_similarity": round(cliff.structural_similarity, 3),
        "cliff_steepness": round(cliff.cliff_steepness, 2)
    }


def result_to_dict(result: SARResult) -> Dict[str, Any]:
    """Convert SARResult to dictionary"""
    return {
        "parent_smiles": result.parent_smiles,
        "scaffold": result.scaffold,
        "r_groups": {str(k): [asdict(rg) for rg in v] for k, v in result.r_groups.items()},
        "activity_cliffs": [cliff_to_dict(c) for c in result.activity_cliffs],
        "pharmacophore_features": result.pharmacophore_features,
        "key_substitutions": result.key_substitutions,
        "execution_time": round(result.execution_time, 3)
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python sar_analyzer.py <smiles>", file=sys.stderr)
        sys.exit(1)
    
    smiles = sys.argv[1]
    
    result = analyze_sar(smiles)
    
    if result:
        output = result_to_dict(result)
        print(json.dumps(output, indent=2))
    else:
        print(json.dumps({"error": "Failed to analyze SAR"}, indent=2), file=sys.stderr)
        sys.exit(1)
