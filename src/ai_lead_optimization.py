#!/usr/bin/env python3
"""
AI-Powered Lead Optimization Module
Uses machine learning to suggest molecular modifications for improved binding
Enhanced with specific molecular modifications and 2D structure visualization
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Draw
from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator
from rdkit import DataStructs
from typing import Dict, List, Tuple
import json
import base64
import io

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
    
    def mol_to_base64_image(self, mol, size=(300, 200), highlight_atoms=None) -> str:
        """Convert molecule to base64 encoded PNG image"""
        try:
            if highlight_atoms:
                img = Draw.MolToImage(mol, size=size, highlightAtoms=highlight_atoms)
            else:
                img = Draw.MolToImage(mol, size=size)
            
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            buffer.seek(0)
            img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            return f"data:image/png;base64,{img_base64}"
        except:
            return ""
    
    def get_compound_specific_modifications(self, mol, smiles: str) -> List[Dict]:
        """Analyze molecule structure and suggest specific, actionable modifications"""
        modifications = []
        
        # Identify structural features
        features = self._identify_structural_features(mol)
        
        # Arylcyclohexylamine modifications (ketamine-like compounds)
        arylcyclohexyl_pattern = Chem.MolFromSmarts('c1ccccc1[C;!$(C=O)]N')
        if arylcyclohexyl_pattern and mol.HasSubstructMatch(arylcyclohexyl_pattern):
            modifications.extend(self._get_arylcyclohexylamine_mods(mol, smiles))
        
        # Aromatic ring modifications
        if features.get('aromatic_rings', 0) >= 1:
            modifications.extend(self._get_aromatic_mods(mol, smiles, features))
        
        # Amine modifications
        if features.get('has_amine', False):
            modifications.extend(self._get_amine_mods(mol, smiles, features))
        
        # Ether/alkoxy modifications
        if features.get('has_ether', False):
            modifications.extend(self._get_ether_mods(mol, smiles, features))
        
        # Alkyl chain modifications
        if len(features.get('alkyl_chains', [])) >= 1:
            modifications.extend(self._get_alkyl_mods(mol, smiles, features))
        
        return modifications
    
    def _identify_structural_features(self, mol) -> Dict:
        """Identify key structural features of the molecule"""
        from rdkit.Chem import Fragments
        
        return {
            'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
            'aliphatic_rings': rdMolDescriptors.CalcNumAliphaticRings(mol),
            'has_amine': Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol) + Fragments.fr_NH0(mol) > 0,
            'has_ether': Fragments.fr_ether(mol) > 0,
            'has_phenol': Fragments.fr_phenol(mol) > 0,
            'has_halogen': Fragments.fr_halogen(mol) > 0,
            'alkyl_chains': mol.GetSubstructMatches(Chem.MolFromSmarts('CC')) if Chem.MolFromSmarts('CC') else [],
            'mw': Descriptors.MolWt(mol),
            'logp': Crippen.MolLogP(mol)
        }
    
    def _get_arylcyclohexylamine_mods(self, mol, smiles: str) -> List[Dict]:
        """Modifications specific to arylcyclohexylamine scaffolds (ketamine, MXE, etc.)"""
        mods = []
        
        # Ortho-methoxy substitution variations
        mods.append({
            "modification_type": "aromatic_substitution",
            "name": "Move methoxy to meta position",
            "description": "Shifting the methoxy group from ortho to meta position can alter receptor binding selectivity and metabolic stability",
            "expected_effect": "May increase sigma receptor affinity while reducing NMDA potency",
            "smarts_from": "COc1ccccc1",
            "smarts_to": "c1cc(OC)ccc1",
            "difficulty": "moderate"
        })
        
        mods.append({
            "modification_type": "aromatic_substitution", 
            "name": "Add para-fluorine",
            "description": "Adding fluorine at para position blocks CYP-mediated metabolism and may increase binding affinity",
            "expected_effect": "Improved metabolic stability, potentially increased potency",
            "smarts_from": "[cH:1]1[cH:2][c:3]([OC])[cH:4][cH:5][c:6]1",
            "smarts_to": "[c:1]1[cH:2][c:3]([OC])[cH:4][c:5](F)[c:6]1",
            "difficulty": "easy"
        })
        
        mods.append({
            "modification_type": "amine_modification",
            "name": "N-ethyl to N-methyl",
            "description": "Reducing N-alkyl chain length typically increases potency but may reduce duration of action",
            "expected_effect": "Higher NMDA affinity, shorter duration",
            "smarts_from": "NCC",
            "smarts_to": "NC",
            "difficulty": "easy"
        })
        
        mods.append({
            "modification_type": "amine_modification",
            "name": "N-cyclopropyl substitution",
            "description": "Cyclopropyl groups can improve metabolic stability while maintaining or increasing receptor affinity",
            "expected_effect": "Improved metabolic stability, altered receptor profile",
            "smarts_from": "[NH1:1]C",
            "smarts_to": "[N:1]C1CC1",
            "difficulty": "moderate"
        })
        
        mods.append({
            "modification_type": "ring_modification",
            "name": "Replace phenyl with 2-thienyl",
            "description": "Thienyl bioisostere can alter receptor selectivity and improve CNS penetration",
            "expected_effect": "Different receptor binding profile, may reduce off-target effects",
            "smarts_from": "c1ccccc1",
            "smarts_to": "c1ccsc1",
            "difficulty": "moderate"
        })
        
        return mods
    
    def _get_aromatic_mods(self, mol, smiles: str, features: Dict) -> List[Dict]:
        """Modifications for aromatic rings"""
        mods = []
        
        if features.get('has_halogen', False):
            mods.append({
                "modification_type": "halogen_exchange",
                "name": "Halogen exchange (Cl → F)",
                "description": "Fluorine is smaller and more electronegative, often improving binding while reducing molecular weight",
                "expected_effect": "Improved binding affinity, better metabolic stability",
                "difficulty": "easy"
            })
        else:
            mods.append({
                "modification_type": "halogenation",
                "name": "Add ortho-fluorine to aromatic ring",
                "description": "Ortho-fluorine blocks metabolic hydroxylation and can improve potency through electronic effects",
                "expected_effect": "Metabolic stability, potential potency increase",
                "difficulty": "easy"
            })
        
        mods.append({
            "modification_type": "ring_replacement",
            "name": "Phenyl to pyridyl replacement",
            "description": "Pyridine can improve water solubility and alter receptor binding through hydrogen bond acceptance",
            "expected_effect": "Improved solubility, altered receptor selectivity",
            "difficulty": "moderate"
        })
        
        return mods
    
    def _get_amine_mods(self, mol, smiles: str, features: Dict) -> List[Dict]:
        """Modifications for amine groups"""
        mods = []
        
        secondary_amine = Chem.MolFromSmarts('[NH1;!$(NC=O)]')
        if secondary_amine and mol.HasSubstructMatch(secondary_amine):
            mods.append({
                "modification_type": "amine_modification",
                "name": "Secondary amine methylation",
                "description": "N-methylation typically increases lipophilicity and CNS penetration, may also affect receptor binding",
                "expected_effect": "Improved CNS penetration, potentially altered receptor profile",
                "difficulty": "easy"
            })
        
        primary_amine = Chem.MolFromSmarts('[NH2]')
        if primary_amine and mol.HasSubstructMatch(primary_amine):
            mods.append({
                "modification_type": "amine_modification",
                "name": "Convert primary amine to dimethylamine",
                "description": "Tertiary amines are more lipophilic and resistant to MAO metabolism",
                "expected_effect": "Metabolic stability, altered receptor binding",
                "difficulty": "easy"
            })
        
        return mods
    
    def _get_ether_mods(self, mol, smiles: str, features: Dict) -> List[Dict]:
        """Modifications for ether groups"""
        return [
            {
                "modification_type": "ether_modification",
                "name": "Methoxy to ethoxy",
                "description": "Ethoxy groups are more lipophilic and may alter binding kinetics",
                "expected_effect": "Increased lipophilicity, potentially slower onset",
                "difficulty": "easy"
            },
            {
                "modification_type": "ether_modification",
                "name": "Methoxy to difluoromethoxy",
                "description": "Difluoromethoxy blocks O-demethylation metabolism while maintaining similar electronic properties",
                "expected_effect": "Significantly improved metabolic stability",
                "difficulty": "moderate"
            }
        ]
    
    def _get_alkyl_mods(self, mol, smiles: str, features: Dict) -> List[Dict]:
        """Modifications for alkyl chains"""
        return [
            {
                "modification_type": "chain_modification",
                "name": "Alkyl chain rigidification",
                "description": "Converting flexible alkyl chains to cyclopropyl or cyclobutyl reduces conformational entropy, potentially improving binding",
                "expected_effect": "Improved binding affinity, metabolic stability",
                "difficulty": "moderate"
            },
            {
                "modification_type": "chain_modification",
                "name": "Alpha-deuteration",
                "description": "Deuterium at alpha-carbon slows oxidative metabolism without affecting pharmacodynamics",
                "expected_effect": "Extended half-life, reduced dosing frequency",
                "difficulty": "moderate"
            }
        ]
    
    def optimize_lead(self, smiles: str, target_profile: Dict = None) -> Dict:
        """Generate optimization suggestions for a lead compound with 2D visualization"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        suggestions = {
            "original_smiles": smiles,
            "original_structure_image": self.mol_to_base64_image(mol, size=(350, 250)),
            "optimization_strategies": [],
            "predicted_improvements": [],
            "modified_structures": [],
            "specific_modifications": [],
            "admet_comparison": {}
        }
        
        # Analyze current molecule
        current_props = self.calculate_properties(mol)
        suggestions["current_properties"] = current_props
        
        # Get compound-specific modifications (the new, detailed suggestions)
        specific_mods = self.get_compound_specific_modifications(mol, smiles)
        suggestions["specific_modifications"] = specific_mods
        
        # Apply specific modifications to generate actual structures
        for mod in specific_mods:
            if "smarts_from" in mod and "smarts_to" in mod:
                try:
                    rxn_smarts = f"{mod['smarts_from']}>>{mod['smarts_to']}"
                    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                    if rxn:
                        products = rxn.RunReactants((mol,))
                        if products:
                            product_mol = products[0][0]
                            Chem.SanitizeMol(product_mol)
                            product_smiles = Chem.MolToSmiles(product_mol)
                            product_props = self.calculate_properties(product_mol)
                            
                            suggestions["modified_structures"].append({
                                "smiles": product_smiles,
                                "structure_image": self.mol_to_base64_image(product_mol, size=(300, 200)),
                                "modification": mod["name"],
                                "modification_type": mod["modification_type"],
                                "description": mod["description"],
                                "expected_effect": mod.get("expected_effect", ""),
                                "difficulty": mod.get("difficulty", "moderate"),
                                "properties": product_props,
                                "property_changes": {
                                    "mw_change": round(product_props["mw"] - current_props["mw"], 2),
                                    "logp_change": round(product_props["logp"] - current_props["logp"], 2),
                                    "tpsa_change": round(product_props["tpsa"] - current_props["tpsa"], 2)
                                }
                            })
                except Exception as e:
                    # If reaction fails, still include the suggestion without structure
                    suggestions["modified_structures"].append({
                        "smiles": None,
                        "structure_image": None,
                        "modification": mod["name"],
                        "modification_type": mod["modification_type"],
                        "description": mod["description"],
                        "expected_effect": mod.get("expected_effect", ""),
                        "difficulty": mod.get("difficulty", "moderate"),
                        "properties": None,
                        "requires_synthesis": True
                    })
        
        # Also identify optimization needs based on properties
        needs = self.identify_optimization_needs(current_props, target_profile)
        
        # Generate strategies for each need
        for need in needs:
            strategy = self.generate_optimization_strategy(mol, need)
            suggestions["optimization_strategies"].append(strategy)
            
            # Try to apply rule-based modifications
            for modification in strategy["modifications"]:
                try:
                    modified_mol = self.apply_modification(mol, modification)
                    if modified_mol and Chem.MolToSmiles(modified_mol) != smiles:
                        Chem.SanitizeMol(modified_mol)
                        modified_smiles = Chem.MolToSmiles(modified_mol)
                        modified_props = self.calculate_properties(modified_mol)
                        
                        improvement = self.assess_improvement(
                            current_props, modified_props, need
                        )
                        
                        # Avoid duplicates
                        existing_smiles = [m.get("smiles") for m in suggestions["modified_structures"]]
                        if modified_smiles not in existing_smiles:
                            suggestions["modified_structures"].append({
                                "smiles": modified_smiles,
                                "structure_image": self.mol_to_base64_image(modified_mol, size=(300, 200)),
                                "modification": modification["description"],
                                "modification_type": modification.get("type", "general"),
                                "description": modification["description"],
                                "expected_effect": f"Addresses: {need}",
                                "difficulty": "easy",
                                "properties": modified_props,
                                "improvement_score": improvement,
                                "strategy": need
                            })
                except:
                    continue
        
        # Sort by those with actual structures first, then by improvement score
        suggestions["modified_structures"] = sorted(
            [m for m in suggestions["modified_structures"] if m.get("smiles")],
            key=lambda x: x.get("improvement_score", 0.5),
            reverse=True
        )[:15]  # Top 15 suggestions
        
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
            # Filter to only include modifications with properties
            mods_with_props = [m for m in modifications if m.get("properties")]
            
            if mods_with_props:
                avg_qed = np.mean([m["properties"]["qed"] for m in mods_with_props])
                if avg_qed > 0.7:
                    improvements.append("Significantly improved drug-likeness")
                
                avg_sa = np.mean([m["properties"]["synthetic_accessibility"] for m in mods_with_props])
                if avg_sa < 4:
                    improvements.append("Highly synthetically accessible analogs")
            
            # Check for best modification
            best_mod = modifications[0]
            mod_name = best_mod.get('modification', 'Unknown')
            
            if 'improvement_score' in best_mod:
                improvements.append(
                    f"Best modification: {mod_name} "
                    f"(score: {best_mod['improvement_score']:.2f})"
                )
            else:
                improvements.append(f"Top suggestion: {mod_name}")
                if best_mod.get('expected_effect'):
                    improvements.append(f"Expected effect: {best_mod['expected_effect']}")
        
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