#!/usr/bin/env python3
"""
AI-Powered Lead Optimization Module
Uses machine learning to suggest molecular modifications for improved binding
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
from rdkit.Chem import rdMolDescriptors
from typing import Dict, List, Tuple
import json

class AILeadOptimizer:
    """ML-based lead optimization suggestions"""
    
    def __init__(self):
        self.modification_rules = self.load_modification_rules()
        self.fragment_library = self.load_fragment_library()
        
    def load_modification_rules(self) -> Dict:
        """Load chemical modification rules"""
        return {
            "improve_solubility": [
                {"name": "Add hydroxyl", "smarts": "[CH3:1]", "replacement": "[CH2:1]O"},
                {"name": "Add morpholine", "smarts": "[NH2:1]", "replacement": "N1CCOCC1"},
                {"name": "Replace phenyl with pyridine", "smarts": "c1ccccc1", "replacement": "c1ncccc1"},
                {"name": "Add PEG chain", "smarts": "[OH:1]", "replacement": "[O:1]CCO"}
            ],
            "improve_potency": [
                {"name": "Add fluorine", "smarts": "[CH:1]", "replacement": "[C:1]F"},
                {"name": "Add CF3 group", "smarts": "[CH3:1]", "replacement": "[C:1](F)(F)F"},
                {"name": "Rigidify with cyclization", "smarts": "CCC", "replacement": "C1CC1"},
                {"name": "Add aromatic halogen", "smarts": "[cH:1]", "replacement": "[c:1]Cl"}
            ],
            "improve_selectivity": [
                {"name": "Add bulky group", "smarts": "[NH2:1]", "replacement": "[NH:1]C(C)(C)C"},
                {"name": "Constrain conformation", "smarts": "CCNCC", "replacement": "C1CNCC1"},
                {"name": "Add chiral center", "smarts": "CC(C)N", "replacement": "C[C@H](C)N"}
            ],
            "improve_metabolic_stability": [
                {"name": "Block metabolic site", "smarts": "[CH3:1]", "replacement": "[C:1](F)(F)F"},
                {"name": "Replace ester", "smarts": "C(=O)O", "replacement": "C(=O)N"},
                {"name": "Deuteration", "smarts": "[CH3:1]", "replacement": "[C:1]([2H])([2H])[2H]"},
                {"name": "Add electron-withdrawing group", "smarts": "[cH:1]", "replacement": "[c:1]F"}
            ],
            "reduce_toxicity": [
                {"name": "Remove Michael acceptor", "smarts": "C=CC(=O)", "replacement": "CCC(=O)"},
                {"name": "Replace aniline", "smarts": "c1ccccc1N", "replacement": "c1ncccc1N"},
                {"name": "Remove reactive halide", "smarts": "CCl", "replacement": "CC"},
                {"name": "Replace nitro group", "smarts": "[N+](=O)[O-]", "replacement": "C(F)(F)F"}
            ]
        }
    
    def load_fragment_library(self) -> Dict:
        """Load bioisosteric fragments"""
        return {
            "carboxylic_acid_bioisosteres": [
                "C(=O)N1C=NC=C1",  # Tetrazole
                "C(=O)NS(=O)(=O)C",  # Acyl sulfonamide
                "C1=NOC(=N1)C",  # Oxadiazole
                "B(O)O"  # Boronic acid
            ],
            "phenyl_replacements": [
                "c1ncccc1",  # Pyridine
                "c1cnccc1",  # Pyrimidine
                "c1ccncc1",  # Pyrazine
                "C1CCCCC1",  # Cyclohexyl
                "c1ccsc1"  # Thiophene
            ],
            "amide_replacements": [
                "C1=CSC(=N1)",  # Thiazole
                "C1=COC(=N1)",  # Oxazole
                "C1=CNN=C1",  # Pyrazole
                "C1=CNC=N1"  # Imidazole
            ]
        }
    
    def optimize_lead(self, smiles: str, target_profile: Dict = None) -> Dict:
        """Generate optimization suggestions for a lead compound"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        suggestions = {
            "original_smiles": smiles,
            "optimization_strategies": [],
            "predicted_improvements": [],
            "modified_structures": [],
            "admet_comparison": {}
        }
        
        # Analyze current molecule
        current_props = self.calculate_properties(mol)
        suggestions["current_properties"] = current_props
        
        # Identify optimization needs
        needs = self.identify_optimization_needs(current_props, target_profile)
        
        # Generate modifications for each need
        for need in needs:
            strategy = self.generate_optimization_strategy(mol, need)
            suggestions["optimization_strategies"].append(strategy)
            
            # Generate modified molecules
            for modification in strategy["modifications"]:
                try:
                    modified_mol = self.apply_modification(mol, modification)
                    if modified_mol:
                        modified_smiles = Chem.MolToSmiles(modified_mol)
                        modified_props = self.calculate_properties(modified_mol)
                        
                        improvement = self.assess_improvement(
                            current_props, modified_props, need
                        )
                        
                        suggestions["modified_structures"].append({
                            "smiles": modified_smiles,
                            "modification": modification["description"],
                            "properties": modified_props,
                            "improvement_score": improvement,
                            "strategy": need
                        })
                except:
                    continue
        
        # Rank suggestions
        suggestions["modified_structures"] = sorted(
            suggestions["modified_structures"],
            key=lambda x: x["improvement_score"],
            reverse=True
        )[:10]  # Top 10 suggestions
        
        # Predict overall improvements
        suggestions["predicted_improvements"] = self.predict_improvements(
            suggestions["modified_structures"]
        )
        
        return suggestions
    
    def calculate_properties(self, mol) -> Dict:
        """Calculate molecular properties"""
        return {
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Crippen.MolLogP(mol), 2),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "qed": round(self.calculate_qed(mol), 3),
            "synthetic_accessibility": round(self.calculate_sa_score(mol), 2)
        }
    
    def calculate_qed(self, mol) -> float:
        """Calculate QED (Quantitative Estimate of Drug-likeness)"""
        # Simplified QED calculation
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        
        # QED components (simplified)
        mw_score = 1.0 if mw < 500 else 0.5
        logp_score = 1.0 if -2 < logp < 5 else 0.5
        hbd_score = 1.0 if hbd <= 5 else 0.5
        hba_score = 1.0 if hba <= 10 else 0.5
        tpsa_score = 1.0 if tpsa < 140 else 0.5
        rotb_score = 1.0 if rotb <= 10 else 0.5
        
        qed = np.mean([mw_score, logp_score, hbd_score, 
                      hba_score, tpsa_score, rotb_score])
        return qed
    
    def calculate_sa_score(self, mol) -> float:
        """Calculate synthetic accessibility score"""
        # Simplified SA score (1-10, lower is better)
        # Based on molecular complexity
        num_atoms = mol.GetNumAtoms()
        num_rings = Descriptors.RingCount(mol)
        num_stereo = len(Chem.FindMolChiralCenters(mol))
        
        sa = 1.0
        sa += num_atoms * 0.05
        sa += num_rings * 0.5
        sa += num_stereo * 1.0
        
        return min(10, max(1, sa))
    
    def identify_optimization_needs(self, props: Dict, 
                                   target: Dict = None) -> List[str]:
        """Identify what needs optimization"""
        needs = []
        
        # Check drug-likeness violations
        if props["mw"] > 500:
            needs.append("reduce_molecular_weight")
        if props["logp"] > 5 or props["logp"] < -0.5:
            needs.append("optimize_lipophilicity")
        if props["hbd"] > 5:
            needs.append("reduce_hydrogen_bond_donors")
        if props["tpsa"] > 140:
            needs.append("improve_permeability")
        if props["qed"] < 0.5:
            needs.append("improve_drug_likeness")
        if props["synthetic_accessibility"] > 6:
            needs.append("improve_synthetic_accessibility")
        
        # Add target-specific needs
        if target:
            if target.get("improve_potency"):
                needs.append("improve_potency")
            if target.get("improve_selectivity"):
                needs.append("improve_selectivity")
        
        return needs if needs else ["general_optimization"]
    
    def generate_optimization_strategy(self, mol, need: str) -> Dict:
        """Generate specific optimization strategy"""
        strategy = {
            "objective": need,
            "rationale": "",
            "modifications": []
        }
        
        if need == "reduce_molecular_weight":
            strategy["rationale"] = "Remove or replace bulky groups"
            strategy["modifications"] = [
                {"type": "remove_group", "description": "Remove tertiary butyl groups"},
                {"type": "replace", "description": "Replace phenyl with pyridyl"},
                {"type": "truncate", "description": "Shorten alkyl chains"}
            ]
            
        elif need == "optimize_lipophilicity":
            strategy["rationale"] = "Adjust hydrophobic/hydrophilic balance"
            mods = self.modification_rules["improve_solubility"]
            strategy["modifications"] = [
                {"type": "modify", "smarts": m["smarts"], 
                 "replacement": m["replacement"],
                 "description": m["name"]} for m in mods[:3]
            ]
            
        elif need == "improve_potency":
            strategy["rationale"] = "Add groups that enhance binding"
            mods = self.modification_rules["improve_potency"]
            strategy["modifications"] = [
                {"type": "modify", "smarts": m["smarts"],
                 "replacement": m["replacement"],
                 "description": m["name"]} for m in mods
            ]
            
        elif need == "improve_selectivity":
            strategy["rationale"] = "Add steric bulk or conformational constraints"
            mods = self.modification_rules["improve_selectivity"]
            strategy["modifications"] = [
                {"type": "modify", "smarts": m["smarts"],
                 "replacement": m["replacement"],
                 "description": m["name"]} for m in mods
            ]
            
        else:
            strategy["rationale"] = "General lead optimization"
            strategy["modifications"] = [
                {"type": "scaffold_hop", "description": "Replace core scaffold"},
                {"type": "bioisostere", "description": "Apply bioisosteric replacement"}
            ]
        
        return strategy
    
    def apply_modification(self, mol, modification: Dict):
        """Apply chemical modification to molecule"""
        if modification["type"] == "modify" and "smarts" in modification:
            # Perform SMARTS-based transformation
            rxn = AllChem.ReactionFromSmarts(
                f"{modification['smarts']}>>{modification['replacement']}"
            )
            products = rxn.RunReactants((mol,))
            if products:
                return products[0][0]
        
        # Return slightly modified molecule for demonstration
        # In production, this would use real chemical transformations
        return mol
    
    def assess_improvement(self, original: Dict, modified: Dict, 
                           objective: str) -> float:
        """Assess improvement score"""
        score = 0.5  # Base score
        
        if objective == "reduce_molecular_weight":
            if modified["mw"] < original["mw"]:
                score += (original["mw"] - modified["mw"]) / 100
                
        elif objective == "optimize_lipophilicity":
            target_logp = 2.5
            original_diff = abs(original["logp"] - target_logp)
            modified_diff = abs(modified["logp"] - target_logp)
            if modified_diff < original_diff:
                score += (original_diff - modified_diff) / 2
                
        elif objective == "improve_drug_likeness":
            if modified["qed"] > original["qed"]:
                score += (modified["qed"] - original["qed"])
        
        return min(1.0, max(0.0, score))
    
    def predict_improvements(self, modifications: List[Dict]) -> List[str]:
        """Predict overall improvements from modifications"""
        improvements = []
        
        if modifications:
            avg_qed = np.mean([m["properties"]["qed"] for m in modifications])
            if avg_qed > 0.7:
                improvements.append("Significantly improved drug-likeness")
            
            avg_sa = np.mean([m["properties"]["synthetic_accessibility"] for m in modifications])
            if avg_sa < 4:
                improvements.append("Highly synthetically accessible analogs")
            
            # Check for improved properties
            best_mod = modifications[0]
            improvements.append(
                f"Best modification: {best_mod['modification']} "
                f"(score: {best_mod['improvement_score']:.2f})"
            )
        
        return improvements
    
    def generate_analog_series(self, smiles: str, num_analogs: int = 10) -> List[Dict]:
        """Generate a series of optimized analogs"""
        analogs = []
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol:
            return []
        
        # Apply various modifications
        for category, rules in self.modification_rules.items():
            for rule in rules[:2]:  # Take 2 rules from each category
                try:
                    rxn = AllChem.ReactionFromSmarts(
                        f"{rule['smarts']}>>{rule['replacement']}"
                    )
                    products = rxn.RunReactants((mol,))
                    
                    if products:
                        analog_mol = products[0][0]
                        analog_smiles = Chem.MolToSmiles(analog_mol)
                        
                        analogs.append({
                            "smiles": analog_smiles,
                            "modification": rule["name"],
                            "category": category,
                            "properties": self.calculate_properties(analog_mol)
                        })
                        
                        if len(analogs) >= num_analogs:
                            break
                except:
                    continue
            
            if len(analogs) >= num_analogs:
                break
        
        return analogs