#!/usr/bin/env python3
"""
Real-time Structure-Activity Relationship (SAR) Explorer
Visualizes how molecular changes affect receptor binding
"""

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem import rdMolDescriptors, rdFMCS
from typing import Dict, List, Tuple
import base64
from io import BytesIO

class SARExplorer:
    """Interactive SAR visualization and analysis"""
    
    def __init__(self):
        self.sar_data = []
        self.activity_cliff_threshold = 0.3
        
    def analyze_sar(self, compound_series: List[Dict]) -> Dict:
        """Analyze SAR for a series of compounds"""
        analysis = {
            "series_size": len(compound_series),
            "sar_matrix": [],
            "activity_cliffs": [],
            "key_features": [],
            "modification_effects": [],
            "interactive_visualization": {}
        }
        
        # Build SAR matrix
        for i, compound in enumerate(compound_series):
            mol = Chem.MolFromSmiles(compound['smiles'])
            if not mol:
                continue
            
            sar_entry = {
                "id": compound.get('id', f"CPD_{i+1}"),
                "smiles": compound['smiles'],
                "activity": compound.get('activity', 0),
                "descriptors": self.calculate_sar_descriptors(mol),
                "fragments": self.identify_key_fragments(mol),
                "modifications": []
            }
            
            # Compare with other compounds
            if i > 0:
                modifications = self.identify_modifications(
                    compound_series[0]['smiles'], 
                    compound['smiles']
                )
                sar_entry["modifications"] = modifications
            
            analysis["sar_matrix"].append(sar_entry)
        
        # Identify activity cliffs
        analysis["activity_cliffs"] = self.identify_activity_cliffs(analysis["sar_matrix"])
        
        # Extract key features
        analysis["key_features"] = self.extract_key_features(analysis["sar_matrix"])
        
        # Analyze modification effects
        analysis["modification_effects"] = self.analyze_modification_effects(analysis["sar_matrix"])
        
        # Generate visualization data
        analysis["interactive_visualization"] = self.generate_visualization_data(analysis["sar_matrix"])
        
        return analysis
    
    def calculate_sar_descriptors(self, mol) -> Dict:
        """Calculate SAR-relevant descriptors"""
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "formal_charge": Chem.rdmolops.GetFormalCharge(mol),
            "num_het_atoms": Descriptors.NumHeteroatoms(mol),
            "ring_count": Descriptors.RingCount(mol),
            "molar_refractivity": Descriptors.MolMR(mol)
        }
    
    def identify_key_fragments(self, mol) -> List[str]:
        """Identify important molecular fragments"""
        fragments = []
        
        # Common pharmacophore patterns
        patterns = {
            "phenyl": "c1ccccc1",
            "pyridine": "c1ccncc1",
            "morpholine": "C1COCCN1",
            "piperidine": "C1CCNCC1",
            "amide": "C(=O)N",
            "ester": "C(=O)O",
            "sulfonamide": "S(=O)(=O)N",
            "amine": "[NX3;H2,H1,H0]",
            "hydroxyl": "[OH]",
            "methoxy": "OC",
            "halogen": "[F,Cl,Br,I]",
            "cyano": "C#N",
            "nitro": "[N+](=O)[O-]",
            "carboxylic": "C(=O)O",
            "ketone": "C(=O)C"
        }
        
        for name, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                fragments.append(name)
        
        return fragments
    
    def identify_modifications(self, ref_smiles: str, query_smiles: str) -> List[Dict]:
        """Identify structural modifications between two molecules"""
        modifications = []
        
        ref_mol = Chem.MolFromSmiles(ref_smiles)
        query_mol = Chem.MolFromSmiles(query_smiles)
        
        if not ref_mol or not query_mol:
            return modifications
        
        # Find maximum common substructure
        mcs = rdFMCS.FindMCS([ref_mol, query_mol])
        
        if mcs:
            common_smarts = mcs.smartsString
            common_mol = Chem.MolFromSmarts(common_smarts)
            
            # Identify what's different
            ref_atoms = ref_mol.GetNumAtoms()
            query_atoms = query_mol.GetNumAtoms()
            common_atoms = common_mol.GetNumAtoms() if common_mol else 0
            
            if ref_atoms != query_atoms:
                modifications.append({
                    "type": "size_change",
                    "description": f"Atom count: {ref_atoms} → {query_atoms}",
                    "impact": "molecular_size"
                })
            
            # Check for specific modifications
            ref_frags = set(self.identify_key_fragments(ref_mol))
            query_frags = set(self.identify_key_fragments(query_mol))
            
            added = query_frags - ref_frags
            removed = ref_frags - query_frags
            
            for frag in added:
                modifications.append({
                    "type": "addition",
                    "fragment": frag,
                    "description": f"Added {frag} group",
                    "impact": "unknown"
                })
            
            for frag in removed:
                modifications.append({
                    "type": "removal",
                    "fragment": frag,
                    "description": f"Removed {frag} group",
                    "impact": "unknown"
                })
        
        return modifications
    
    def identify_activity_cliffs(self, sar_matrix: List[Dict]) -> List[Dict]:
        """Identify activity cliffs - similar structures with different activities"""
        cliffs = []
        
        for i in range(len(sar_matrix)):
            for j in range(i + 1, len(sar_matrix)):
                mol1 = Chem.MolFromSmiles(sar_matrix[i]['smiles'])
                mol2 = Chem.MolFromSmiles(sar_matrix[j]['smiles'])
                
                if mol1 and mol2:
                    # Calculate similarity
                    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
                    
                    from rdkit import DataStructs
                    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                    
                    # Check activity difference
                    activity_diff = abs(sar_matrix[i]['activity'] - sar_matrix[j]['activity'])
                    
                    # Activity cliff: high similarity, large activity difference
                    if similarity > 0.8 and activity_diff > self.activity_cliff_threshold:
                        cliffs.append({
                            "compound1": sar_matrix[i]['id'],
                            "compound2": sar_matrix[j]['id'],
                            "similarity": round(similarity, 3),
                            "activity_difference": round(activity_diff, 3),
                            "cliff_steepness": round(activity_diff / (1 - similarity + 0.001), 2),
                            "modifications": self.identify_modifications(
                                sar_matrix[i]['smiles'],
                                sar_matrix[j]['smiles']
                            )
                        })
        
        # Sort by cliff steepness
        cliffs.sort(key=lambda x: x['cliff_steepness'], reverse=True)
        return cliffs
    
    def extract_key_features(self, sar_matrix: List[Dict]) -> List[Dict]:
        """Extract molecular features that correlate with activity"""
        features = []
        
        if not sar_matrix:
            return features
        
        # Collect all fragments and their activities
        fragment_activities = {}
        
        for entry in sar_matrix:
            for fragment in entry['fragments']:
                if fragment not in fragment_activities:
                    fragment_activities[fragment] = []
                fragment_activities[fragment].append(entry['activity'])
        
        # Analyze fragment contributions
        for fragment, activities in fragment_activities.items():
            if len(activities) > 1:
                avg_activity = np.mean(activities)
                std_activity = np.std(activities)
                
                # Compare with compounds without this fragment
                without_fragment = []
                for entry in sar_matrix:
                    if fragment not in entry['fragments']:
                        without_fragment.append(entry['activity'])
                
                if without_fragment:
                    avg_without = np.mean(without_fragment)
                    contribution = avg_activity - avg_without
                    
                    features.append({
                        "fragment": fragment,
                        "occurrence": len(activities),
                        "avg_activity": round(avg_activity, 3),
                        "contribution": round(contribution, 3),
                        "consistency": round(1 / (std_activity + 0.01), 2),
                        "importance": "positive" if contribution > 0.1 else "negative" if contribution < -0.1 else "neutral"
                    })
        
        # Sort by absolute contribution
        features.sort(key=lambda x: abs(x['contribution']), reverse=True)
        return features[:10]  # Top 10 features
    
    def analyze_modification_effects(self, sar_matrix: List[Dict]) -> List[Dict]:
        """Analyze how specific modifications affect activity"""
        effects = []
        
        # Track modification impacts
        modification_impacts = {}
        
        for i in range(1, len(sar_matrix)):
            mods = sar_matrix[i].get('modifications', [])
            activity_change = sar_matrix[i]['activity'] - sar_matrix[0]['activity']
            
            for mod in mods:
                key = f"{mod['type']}_{mod.get('fragment', 'unknown')}"
                if key not in modification_impacts:
                    modification_impacts[key] = []
                modification_impacts[key].append(activity_change)
        
        # Summarize impacts
        for mod_key, impacts in modification_impacts.items():
            if impacts:
                effects.append({
                    "modification": mod_key,
                    "occurrences": len(impacts),
                    "avg_impact": round(np.mean(impacts), 3),
                    "impact_range": [round(min(impacts), 3), round(max(impacts), 3)],
                    "consistency": round(1 / (np.std(impacts) + 0.01), 2),
                    "recommendation": "beneficial" if np.mean(impacts) > 0.1 else "detrimental" if np.mean(impacts) < -0.1 else "neutral"
                })
        
        # Sort by average impact
        effects.sort(key=lambda x: abs(x['avg_impact']), reverse=True)
        return effects
    
    def generate_visualization_data(self, sar_matrix: List[Dict]) -> Dict:
        """Generate data for interactive SAR visualization"""
        viz_data = {
            "compounds": [],
            "activity_range": [0, 0],
            "descriptor_ranges": {},
            "heatmap_data": [],
            "structure_images": []
        }
        
        activities = []
        descriptors_all = {key: [] for key in ['mw', 'logp', 'hbd', 'hba', 'tpsa']}
        
        for entry in sar_matrix:
            activities.append(entry['activity'])
            
            # Collect descriptor values
            for key in descriptors_all:
                if key in entry['descriptors']:
                    descriptors_all[key].append(entry['descriptors'][key])
            
            # Generate structure image (base64)
            mol = Chem.MolFromSmiles(entry['smiles'])
            if mol:
                img = Draw.MolToImage(mol, size=(200, 200))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                
                viz_data["structure_images"].append({
                    "id": entry['id'],
                    "image": f"data:image/png;base64,{img_str}"
                })
            
            viz_data["compounds"].append({
                "id": entry['id'],
                "smiles": entry['smiles'],
                "activity": entry['activity'],
                "descriptors": entry['descriptors']
            })
        
        # Set ranges
        if activities:
            viz_data["activity_range"] = [min(activities), max(activities)]
        
        for key, values in descriptors_all.items():
            if values:
                viz_data["descriptor_ranges"][key] = [min(values), max(values)]
        
        # Generate heatmap data for descriptor-activity correlation
        for desc_key in ['mw', 'logp', 'tpsa']:
            if desc_key in viz_data["descriptor_ranges"]:
                correlation_data = []
                for entry in sar_matrix:
                    if desc_key in entry['descriptors']:
                        correlation_data.append({
                            "descriptor": desc_key,
                            "value": entry['descriptors'][desc_key],
                            "activity": entry['activity']
                        })
                viz_data["heatmap_data"].append(correlation_data)
        
        return viz_data
    
    def predict_activity(self, smiles: str, sar_model: Dict = None) -> Dict:
        """Predict activity based on SAR analysis"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        prediction = {
            "smiles": smiles,
            "predicted_activity": 0.5,  # Default
            "confidence": 0.0,
            "similar_compounds": [],
            "key_features_present": [],
            "missing_features": []
        }
        
        # Extract features
        fragments = self.identify_key_fragments(mol)
        descriptors = self.calculate_sar_descriptors(mol)
        
        # Simple prediction based on fragments (in production, use ML model)
        activity_score = 0.5
        confidence = 0.0
        
        # Check for positive features
        positive_features = ['phenyl', 'pyridine', 'amide', 'morpholine']
        negative_features = ['nitro', 'halogen']
        
        for feat in positive_features:
            if feat in fragments:
                activity_score += 0.1
                confidence += 0.15
                prediction["key_features_present"].append(feat)
        
        for feat in negative_features:
            if feat in fragments:
                activity_score -= 0.05
                confidence += 0.1
        
        # Apply descriptor-based rules
        if descriptors['mw'] < 500 and descriptors['logp'] < 5:
            activity_score += 0.1
            confidence += 0.1
        
        prediction["predicted_activity"] = round(min(1.0, max(0.0, activity_score)), 3)
        prediction["confidence"] = round(min(1.0, confidence), 2) * 100
        
        # Identify missing beneficial features
        for feat in positive_features:
            if feat not in fragments:
                prediction["missing_features"].append(feat)
        
        return prediction