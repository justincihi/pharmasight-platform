#!/usr/bin/env python3
"""
Pharmacophore Modeling Module
Identifies key molecular features required for receptor binding
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from typing import Dict, List, Tuple
import json

class PharmacophoreModeler:
    """3D Pharmacophore generation and screening"""
    
    def __init__(self):
        self.pharmacophore_definitions = self.load_pharmacophore_definitions()
        self.feature_definitions = self.define_pharmacophoric_features()
        
    def load_pharmacophore_definitions(self) -> Dict:
        """Load predefined pharmacophore models for known targets"""
        return {
            "5-HT2A": {
                "features": [
                    {"type": "Aromatic", "required": True, "position": [0, 0, 0]},
                    {"type": "HBondAcceptor", "required": True, "position": [3.5, 1.2, 0]},
                    {"type": "Hydrophobic", "required": False, "position": [5.0, 0, 0]},
                    {"type": "PositiveIonizable", "required": True, "position": [7.0, 2.0, 0]}
                ],
                "distances": {
                    "Aromatic-HBondAcceptor": [3.0, 4.5],
                    "Aromatic-PositiveIonizable": [6.5, 8.0],
                    "HBondAcceptor-PositiveIonizable": [3.5, 5.0]
                }
            },
            "D2": {
                "features": [
                    {"type": "Aromatic", "required": True, "position": [0, 0, 0]},
                    {"type": "PositiveIonizable", "required": True, "position": [5.5, 1.0, 0]},
                    {"type": "HBondDonor", "required": False, "position": [3.0, 3.0, 0]}
                ],
                "distances": {
                    "Aromatic-PositiveIonizable": [5.0, 6.5],
                    "Aromatic-HBondDonor": [2.5, 4.0]
                }
            },
            "MOR": {
                "features": [
                    {"type": "Aromatic", "required": True, "position": [0, 0, 0]},
                    {"type": "PositiveIonizable", "required": True, "position": [4.5, 2.0, 0]},
                    {"type": "HBondAcceptor", "required": True, "position": [2.0, 3.5, 0]},
                    {"type": "Hydrophobic", "required": True, "position": [6.0, 0, 1.0]}
                ],
                "distances": {
                    "Aromatic-PositiveIonizable": [4.0, 5.5],
                    "PositiveIonizable-HBondAcceptor": [2.5, 4.0],
                    "Aromatic-Hydrophobic": [5.5, 7.0]
                }
            },
            "GABA-A": {
                "features": [
                    {"type": "Aromatic", "required": True, "position": [0, 0, 0]},
                    {"type": "HBondAcceptor", "required": True, "position": [3.0, 1.5, 0]},
                    {"type": "Hydrophobic", "required": True, "position": [5.5, 0, 0]},
                    {"type": "Aromatic", "required": False, "position": [3.0, -3.0, 0]}
                ],
                "distances": {
                    "Aromatic-HBondAcceptor": [2.5, 3.5],
                    "Aromatic-Hydrophobic": [5.0, 6.5]
                }
            },
            "Beta2": {
                "features": [
                    {"type": "HBondDonor", "required": True, "position": [0, 0, 0]},
                    {"type": "PositiveIonizable", "required": True, "position": [3.0, 1.0, 0]},
                    {"type": "Aromatic", "required": False, "position": [5.0, 0, 0]}
                ],
                "distances": {
                    "HBondDonor-PositiveIonizable": [2.5, 3.5],
                    "PositiveIonizable-Aromatic": [4.0, 6.0]
                }
            }
        }
    
    def define_pharmacophoric_features(self) -> Dict:
        """Define SMARTS patterns for pharmacophoric features"""
        return {
            "Aromatic": ["a", "c1ccccc1", "c1ccncc1", "c1ccccn1"],
            "HBondDonor": ["[OH]", "[NH]", "[NH2]", "[nH]"],
            "HBondAcceptor": ["[O]", "[N]", "C=O", "S=O", "[nX2]"],
            "PositiveIonizable": ["[NH3+]", "[NH2+]", "[NH+]", "[N+]", "C(N)N"],
            "NegativeIonizable": ["C(=O)[O-]", "S(=O)(=O)[O-]", "P(=O)([O-])[O-]"],
            "Hydrophobic": ["[CH3]", "[CH2][CH3]", "C(C)(C)C", "c1ccccc1"],
            "Halogen": ["[F]", "[Cl]", "[Br]", "[I]"],
            "MetalCoordination": ["[O-]", "[S-]", "[N]", "C=O", "C=S"]
        }
    
    def extract_pharmacophore(self, active_compounds: List[Dict]) -> Dict:
        """Extract common pharmacophore from active compounds"""
        pharmacophore = {
            "num_compounds": len(active_compounds),
            "common_features": [],
            "distance_constraints": [],
            "excluded_volumes": [],
            "confidence": 0.0
        }
        
        if len(active_compounds) < 2:
            return {"error": "Need at least 2 active compounds"}
        
        # Extract features from each compound
        all_features = []
        for compound in active_compounds:
            mol = Chem.MolFromSmiles(compound['smiles'])
            if mol:
                features = self.extract_molecular_features(mol)
                all_features.append(features)
        
        # Find common features
        common_features = self.find_common_features(all_features)
        pharmacophore["common_features"] = common_features
        
        # Calculate distance constraints
        if len(common_features) >= 2:
            pharmacophore["distance_constraints"] = self.calculate_distance_constraints(
                active_compounds, common_features
            )
        
        # Estimate confidence
        pharmacophore["confidence"] = self.calculate_pharmacophore_confidence(
            all_features, common_features
        )
        
        # Generate excluded volumes (regions where bulk is not tolerated)
        pharmacophore["excluded_volumes"] = self.identify_excluded_volumes(active_compounds)
        
        return pharmacophore
    
    def extract_molecular_features(self, mol) -> List[Dict]:
        """Extract pharmacophoric features from a molecule"""
        features = []
        
        for feature_type, patterns in self.feature_definitions.items():
            for pattern_smarts in patterns:
                pattern = Chem.MolFromSmarts(pattern_smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    matches = mol.GetSubstructMatches(pattern)
                    for match in matches:
                        # Get 3D coordinates if available
                        if mol.GetNumConformers() > 0:
                            conf = mol.GetConformer()
                            pos = conf.GetAtomPosition(match[0])
                            coords = [pos.x, pos.y, pos.z]
                        else:
                            # Generate 2D coordinates as fallback
                            AllChem.Compute2DCoords(mol)
                            coords = [0, 0, 0]  # Simplified
                        
                        features.append({
                            "type": feature_type,
                            "atoms": match,
                            "position": coords,
                            "pattern": pattern_smarts
                        })
                    break  # Found this feature type, move to next
        
        return features
    
    def find_common_features(self, all_features: List[List[Dict]]) -> List[Dict]:
        """Find features common to all or most compounds"""
        if not all_features:
            return []
        
        # Count feature occurrences
        feature_counts = {}
        total_compounds = len(all_features)
        
        for compound_features in all_features:
            seen_types = set()
            for feature in compound_features:
                feat_type = feature["type"]
                if feat_type not in seen_types:
                    seen_types.add(feat_type)
                    feature_counts[feat_type] = feature_counts.get(feat_type, 0) + 1
        
        # Features present in at least 75% of compounds
        threshold = 0.75 * total_compounds
        common_features = []
        
        for feat_type, count in feature_counts.items():
            if count >= threshold:
                common_features.append({
                    "type": feat_type,
                    "occurrence": count / total_compounds,
                    "required": count == total_compounds
                })
        
        return common_features
    
    def calculate_distance_constraints(self, compounds: List[Dict], 
                                      features: List[Dict]) -> List[Dict]:
        """Calculate distance ranges between pharmacophoric features"""
        constraints = []
        
        if len(features) < 2:
            return constraints
        
        # For each pair of features
        for i in range(len(features)):
            for j in range(i + 1, len(features)):
                feat1 = features[i]["type"]
                feat2 = features[j]["type"]
                
                distances = []
                # Calculate distances in each compound
                # Simplified - in production would use 3D conformers
                
                # Generate mock distances based on feature types
                if feat1 == "Aromatic" and feat2 == "PositiveIonizable":
                    distances = [5.5, 6.0, 5.8]  # Mock data
                elif feat1 == "HBondAcceptor" and feat2 == "HBondDonor":
                    distances = [2.8, 3.2, 3.0]
                else:
                    distances = [4.0, 4.5, 4.2]
                
                if distances:
                    constraints.append({
                        "feature1": feat1,
                        "feature2": feat2,
                        "min_distance": round(min(distances), 1),
                        "max_distance": round(max(distances), 1),
                        "avg_distance": round(np.mean(distances), 1),
                        "tolerance": round(np.std(distances), 2)
                    })
        
        return constraints
    
    def identify_excluded_volumes(self, compounds: List[Dict]) -> List[Dict]:
        """Identify regions where bulk is not tolerated"""
        excluded_volumes = []
        
        # Simplified excluded volume identification
        # In production, would analyze inactive compounds vs active
        
        excluded_volumes.append({
            "position": [8.0, 3.0, 0],
            "radius": 2.0,
            "description": "Steric clash region"
        })
        
        return excluded_volumes
    
    def calculate_pharmacophore_confidence(self, all_features: List[List[Dict]], 
                                          common_features: List[Dict]) -> float:
        """Calculate confidence score for pharmacophore model"""
        if not all_features or not common_features:
            return 0.0
        
        # Base confidence on feature consistency
        total_possible = len(all_features) * len(common_features)
        total_found = sum(len(f) for f in all_features)
        
        consistency = min(1.0, total_found / (total_possible + 1))
        
        # Adjust for number of compounds
        compound_factor = min(1.0, len(all_features) / 5)  # Max confidence at 5+ compounds
        
        # Adjust for feature diversity
        diversity_factor = min(1.0, len(common_features) / 3)  # Good if 3+ features
        
        confidence = (consistency * 0.5 + compound_factor * 0.3 + diversity_factor * 0.2) * 100
        
        return round(confidence, 1)
    
    def screen_compound(self, smiles: str, pharmacophore: Dict) -> Dict:
        """Screen a compound against a pharmacophore model"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        screening_result = {
            "smiles": smiles,
            "matches_pharmacophore": False,
            "matched_features": [],
            "missing_features": [],
            "distance_violations": [],
            "fit_score": 0.0
        }
        
        # Extract features
        mol_features = self.extract_molecular_features(mol)
        mol_feature_types = set(f["type"] for f in mol_features)
        
        # Check required features
        for pharm_feature in pharmacophore.get("common_features", []):
            if pharm_feature["required"]:
                if pharm_feature["type"] in mol_feature_types:
                    screening_result["matched_features"].append(pharm_feature["type"])
                else:
                    screening_result["missing_features"].append(pharm_feature["type"])
        
        # Check distance constraints (simplified)
        if len(screening_result["matched_features"]) >= 2:
            # In production, would calculate actual 3D distances
            # For now, assume some violations for demonstration
            if np.random.random() > 0.7:
                screening_result["distance_violations"].append({
                    "constraint": "Aromatic-PositiveIonizable",
                    "expected": "5.5-6.5 Å",
                    "found": "7.2 Å"
                })
        
        # Calculate fit score
        if not screening_result["missing_features"]:
            screening_result["matches_pharmacophore"] = True
            base_score = 1.0
        else:
            base_score = len(screening_result["matched_features"]) / (
                len(screening_result["matched_features"]) + 
                len(screening_result["missing_features"])
            )
        
        # Penalize for distance violations
        violation_penalty = 0.1 * len(screening_result["distance_violations"])
        screening_result["fit_score"] = round(max(0, base_score - violation_penalty), 2)
        
        return screening_result
    
    def generate_pharmacophore_query(self, receptor: str) -> Dict:
        """Generate a pharmacophore query for a specific receptor"""
        if receptor in self.pharmacophore_definitions:
            pharm = self.pharmacophore_definitions[receptor]
            return {
                "receptor": receptor,
                "pharmacophore": pharm,
                "description": f"Pharmacophore model for {receptor} receptor",
                "num_features": len(pharm["features"]),
                "num_constraints": len(pharm["distances"])
            }
        else:
            # Generate generic pharmacophore
            return {
                "receptor": receptor,
                "pharmacophore": {
                    "features": [
                        {"type": "Aromatic", "required": True},
                        {"type": "HBondAcceptor", "required": False},
                        {"type": "Hydrophobic", "required": False}
                    ],
                    "distances": {}
                },
                "description": f"Generic pharmacophore for {receptor}",
                "num_features": 3,
                "num_constraints": 0
            }
    
    def optimize_pharmacophore(self, pharmacophore: Dict, 
                               validation_compounds: List[Dict]) -> Dict:
        """Refine pharmacophore model based on validation set"""
        optimized = pharmacophore.copy()
        
        # Test each compound against pharmacophore
        results = []
        for compound in validation_compounds:
            result = self.screen_compound(compound['smiles'], pharmacophore)
            result['activity'] = compound.get('activity', 0)
            results.append(result)
        
        # Analyze results to optimize
        active_results = [r for r in results if r.get('activity', 0) > 0.7]
        inactive_results = [r for r in results if r.get('activity', 0) < 0.3]
        
        # Identify features that discriminate actives from inactives
        if active_results and inactive_results:
            # Features more common in actives
            active_features = set()
            for r in active_results:
                active_features.update(r['matched_features'])
            
            inactive_features = set()
            for r in inactive_results:
                inactive_features.update(r['matched_features'])
            
            discriminating = active_features - inactive_features
            
            # Update pharmacophore
            for feature in optimized.get('common_features', []):
                if feature['type'] in discriminating:
                    feature['importance'] = 'high'
                    feature['required'] = True
        
        optimized['optimized'] = True
        optimized['validation_accuracy'] = self.calculate_validation_accuracy(results)
        
        return optimized
    
    def calculate_validation_accuracy(self, results: List[Dict]) -> float:
        """Calculate accuracy of pharmacophore predictions"""
        correct = 0
        total = len(results)
        
        for result in results:
            predicted_active = result['fit_score'] > 0.6
            actual_active = result.get('activity', 0) > 0.5
            
            if predicted_active == actual_active:
                correct += 1
        
        return round(correct / total * 100, 1) if total > 0 else 0.0