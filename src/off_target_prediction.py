#!/usr/bin/env python3
"""
Off-Target Prediction System
Predicts potential side effects based on unintended receptor binding
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from typing import Dict, List, Tuple
import json

class OffTargetPredictor:
    """Predict off-target effects and safety issues"""
    
    def __init__(self):
        self.safety_panel = self.load_safety_panel()
        self.side_effect_map = self.load_side_effect_map()
        self.risk_thresholds = self.define_risk_thresholds()
        
    def load_safety_panel(self) -> Dict:
        """Load standard safety pharmacology panel"""
        return {
            "cardiac": {
                "hERG": {
                    "full_name": "Human Ether-à-go-go-Related Gene",
                    "concern": "QT prolongation, arrhythmia",
                    "ic50_threshold": 10.0,  # µM
                    "severity": "critical"
                },
                "Nav1.5": {
                    "full_name": "Cardiac sodium channel",
                    "concern": "Cardiac conduction abnormalities",
                    "ic50_threshold": 10.0,
                    "severity": "high"
                },
                "Cav1.2": {
                    "full_name": "L-type calcium channel",
                    "concern": "Negative inotropic effects",
                    "ic50_threshold": 10.0,
                    "severity": "high"
                },
                "5-HT2B": {
                    "full_name": "Serotonin 2B receptor",
                    "concern": "Cardiac valvulopathy",
                    "ic50_threshold": 1.0,
                    "severity": "critical"
                }
            },
            "cns": {
                "D2": {
                    "full_name": "Dopamine D2 receptor",
                    "concern": "Extrapyramidal symptoms",
                    "ic50_threshold": 1.0,
                    "severity": "high"
                },
                "H1": {
                    "full_name": "Histamine H1 receptor",
                    "concern": "Sedation, weight gain",
                    "ic50_threshold": 10.0,
                    "severity": "moderate"
                },
                "M1": {
                    "full_name": "Muscarinic M1 receptor",
                    "concern": "Cognitive impairment, confusion",
                    "ic50_threshold": 10.0,
                    "severity": "moderate"
                },
                "GABA-A": {
                    "full_name": "GABA-A receptor",
                    "concern": "Sedation, dependence",
                    "ic50_threshold": 10.0,
                    "severity": "moderate"
                }
            },
            "peripheral": {
                "Alpha1A": {
                    "full_name": "α1A-adrenergic receptor",
                    "concern": "Orthostatic hypotension",
                    "ic50_threshold": 1.0,
                    "severity": "moderate"
                },
                "M3": {
                    "full_name": "Muscarinic M3 receptor",
                    "concern": "Dry mouth, constipation, urinary retention",
                    "ic50_threshold": 10.0,
                    "severity": "low"
                },
                "Beta2": {
                    "full_name": "β2-adrenergic receptor",
                    "concern": "Tremor, tachycardia",
                    "ic50_threshold": 10.0,
                    "severity": "low"
                }
            },
            "metabolic": {
                "PPAR-gamma": {
                    "full_name": "Peroxisome proliferator-activated receptor γ",
                    "concern": "Weight gain, edema",
                    "ic50_threshold": 10.0,
                    "severity": "moderate"
                },
                "CYP3A4": {
                    "full_name": "Cytochrome P450 3A4",
                    "concern": "Drug-drug interactions",
                    "ic50_threshold": 1.0,
                    "severity": "high"
                },
                "CYP2D6": {
                    "full_name": "Cytochrome P450 2D6",
                    "concern": "Drug-drug interactions",
                    "ic50_threshold": 1.0,
                    "severity": "high"
                }
            },
            "hepatic": {
                "BSEP": {
                    "full_name": "Bile Salt Export Pump",
                    "concern": "Cholestasis, hepatotoxicity",
                    "ic50_threshold": 50.0,
                    "severity": "high"
                },
                "MRP2": {
                    "full_name": "Multidrug Resistance Protein 2",
                    "concern": "Hyperbilirubinemia",
                    "ic50_threshold": 50.0,
                    "severity": "moderate"
                }
            }
        }
    
    def load_side_effect_map(self) -> Dict:
        """Map receptor interactions to clinical side effects"""
        return {
            "hERG": ["QT prolongation", "Torsades de pointes", "Sudden cardiac death"],
            "5-HT2B": ["Valvular heart disease", "Pulmonary hypertension"],
            "D2": ["Parkinsonism", "Akathisia", "Tardive dyskinesia", "Hyperprolactinemia"],
            "H1": ["Sedation", "Weight gain", "Dry mouth", "Dizziness"],
            "M1": ["Memory impairment", "Confusion", "Hallucinations"],
            "M3": ["Dry mouth", "Constipation", "Urinary retention", "Blurred vision"],
            "Alpha1A": ["Orthostatic hypotension", "Dizziness", "Syncope"],
            "GABA-A": ["Sedation", "Ataxia", "Amnesia", "Dependence"],
            "CYP3A4": ["Drug interactions", "Altered drug metabolism"],
            "CYP2D6": ["Drug interactions", "Poor/ultra-rapid metabolism"],
            "BSEP": ["Cholestasis", "Jaundice", "Liver injury"]
        }
    
    def define_risk_thresholds(self) -> Dict:
        """Define risk assessment thresholds"""
        return {
            "binding_affinity": {
                "high_risk": 0.1,  # IC50 < 0.1 µM
                "moderate_risk": 1.0,  # IC50 < 1 µM
                "low_risk": 10.0  # IC50 < 10 µM
            },
            "selectivity": {
                "poor": 10,  # <10-fold selectivity
                "moderate": 100,  # <100-fold
                "good": 1000  # >1000-fold
            }
        }
    
    def predict_off_targets(self, smiles: str, primary_target: str = None) -> Dict:
        """Predict off-target interactions for a compound"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        prediction = {
            "smiles": smiles,
            "primary_target": primary_target,
            "off_target_hits": [],
            "safety_alerts": [],
            "predicted_side_effects": [],
            "risk_score": 0,
            "safety_classification": "",
            "drug_drug_interaction_risk": ""
        }
        
        # Generate molecular fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        
        # Calculate molecular properties
        props = self.calculate_safety_properties(mol)
        
        # Check each safety panel target
        total_risk = 0
        critical_alerts = []
        high_alerts = []
        moderate_alerts = []
        
        for category, targets in self.safety_panel.items():
            for target_name, target_info in targets.items():
                # Predict binding (simplified - in production use ML model)
                binding_prob = self.predict_target_binding(mol, fp, target_name)
                predicted_ic50 = self.probability_to_ic50(binding_prob)
                
                if predicted_ic50 < target_info["ic50_threshold"]:
                    hit = {
                        "target": target_name,
                        "category": category,
                        "predicted_ic50": f"{predicted_ic50:.2f} µM",
                        "binding_probability": f"{binding_prob*100:.1f}%",
                        "concern": target_info["concern"],
                        "severity": target_info["severity"]
                    }
                    prediction["off_target_hits"].append(hit)
                    
                    # Add to risk score
                    if target_info["severity"] == "critical":
                        total_risk += 10
                        critical_alerts.append(target_name)
                    elif target_info["severity"] == "high":
                        total_risk += 5
                        high_alerts.append(target_name)
                    elif target_info["severity"] == "moderate":
                        total_risk += 2
                        moderate_alerts.append(target_name)
                    else:
                        total_risk += 1
                    
                    # Add predicted side effects
                    if target_name in self.side_effect_map:
                        for effect in self.side_effect_map[target_name]:
                            if effect not in prediction["predicted_side_effects"]:
                                prediction["predicted_side_effects"].append(effect)
        
        # Generate safety alerts
        if critical_alerts:
            prediction["safety_alerts"].append({
                "level": "CRITICAL",
                "message": f"Critical safety liability: {', '.join(critical_alerts)}",
                "recommendation": "Consider structural modification or abandon compound"
            })
        
        if high_alerts:
            prediction["safety_alerts"].append({
                "level": "HIGH",
                "message": f"Significant off-targets: {', '.join(high_alerts)}",
                "recommendation": "Requires careful monitoring and selectivity optimization"
            })
        
        if moderate_alerts:
            prediction["safety_alerts"].append({
                "level": "MODERATE",
                "message": f"Moderate concerns: {', '.join(moderate_alerts)}",
                "recommendation": "Monitor in preclinical studies"
            })
        
        # Calculate overall risk score and classification
        prediction["risk_score"] = min(100, total_risk * 5)
        
        if prediction["risk_score"] > 70:
            prediction["safety_classification"] = "High Risk"
        elif prediction["risk_score"] > 40:
            prediction["safety_classification"] = "Moderate Risk"
        elif prediction["risk_score"] > 20:
            prediction["safety_classification"] = "Low Risk"
        else:
            prediction["safety_classification"] = "Acceptable"
        
        # Assess drug-drug interaction risk
        cyp_hits = [h for h in prediction["off_target_hits"] 
                   if h["target"] in ["CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19"]]
        
        if len(cyp_hits) >= 2:
            prediction["drug_drug_interaction_risk"] = "High"
        elif len(cyp_hits) == 1:
            prediction["drug_drug_interaction_risk"] = "Moderate"
        else:
            prediction["drug_drug_interaction_risk"] = "Low"
        
        # Add structural alerts
        structural_alerts = self.check_structural_alerts(mol)
        if structural_alerts:
            prediction["safety_alerts"].extend(structural_alerts)
        
        return prediction
    
    def calculate_safety_properties(self, mol) -> Dict:
        """Calculate properties relevant to safety"""
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
            "num_heteroatoms": Descriptors.NumHeteroatoms(mol),
            "formal_charge": Chem.rdmolops.GetFormalCharge(mol)
        }
    
    def predict_target_binding(self, mol, fp, target: str) -> float:
        """Predict binding probability to a specific target"""
        # Simplified prediction - in production use trained ML model
        # Using heuristics based on molecular properties
        
        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "aromatic": Descriptors.NumAromaticRings(mol)
        }
        
        base_prob = 0.1
        
        # Target-specific rules
        if target == "hERG":
            # hERG blockers tend to be basic, lipophilic, with aromatic rings
            if props["logp"] > 3 and props["aromatic"] > 0:
                base_prob += 0.3
            if props["mw"] > 350:
                base_prob += 0.2
        
        elif target == "5-HT2B":
            # 5-HT2B agonists often have indole or phenethylamine scaffolds
            if props["aromatic"] > 1 and props["hba"] > 1:
                base_prob += 0.25
        
        elif target == "D2":
            # D2 ligands often have aromatic rings and basic nitrogen
            if props["aromatic"] > 0 and props["hba"] > 1:
                base_prob += 0.2
        
        elif target == "CYP3A4":
            # Large, lipophilic molecules
            if props["mw"] > 400 and props["logp"] > 3:
                base_prob += 0.3
        
        elif target == "BSEP":
            # Bile acid-like structures
            if props["mw"] > 350 and props["hbd"] > 2:
                base_prob += 0.2
        
        # Add randomness for simulation
        import random
        base_prob += random.uniform(-0.1, 0.2)
        
        return max(0, min(1, base_prob))
    
    def probability_to_ic50(self, probability: float) -> float:
        """Convert binding probability to predicted IC50"""
        if probability > 0.8:
            return 0.01  # 10 nM
        elif probability > 0.6:
            return 0.1  # 100 nM
        elif probability > 0.4:
            return 1.0  # 1 µM
        elif probability > 0.2:
            return 10.0  # 10 µM
        else:
            return 100.0  # 100 µM
    
    def check_structural_alerts(self, mol) -> List[Dict]:
        """Check for problematic structural features"""
        alerts = []
        
        # PAINS (Pan-Assay Interference Compounds) patterns
        pains_patterns = {
            "quinone": "[O,N]=C1C=CC(=O)C=C1",
            "rhodanine": "C1(=O)NC(=S)SC1",
            "hydroxyphenyl_hydrazone": "[OH]c1ccccc1/N=N/",
            "michael_acceptor": "C=CC(=O)",
            "alkyl_halide": "[C][Cl,Br,I]",
            "aldehyde": "[CH](=O)",
            "epoxide": "C1OC1",
            "aziridine": "C1NC1"
        }
        
        for name, smarts in pains_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                alerts.append({
                    "level": "WARNING",
                    "message": f"Structural alert: {name.replace('_', ' ').title()}",
                    "recommendation": "May cause assay interference or reactivity"
                })
        
        return alerts
    
    def calculate_selectivity(self, primary_ic50: float, 
                            off_target_ic50s: List[float]) -> Dict:
        """Calculate selectivity indices"""
        if not off_target_ic50s:
            return {
                "status": "No off-targets detected",
                "selectivity_window": ">1000-fold",
                "classification": "Highly selective"
            }
        
        min_off_target = min(off_target_ic50s)
        selectivity_index = min_off_target / primary_ic50
        
        if selectivity_index > 1000:
            classification = "Highly selective"
        elif selectivity_index > 100:
            classification = "Selective"
        elif selectivity_index > 10:
            classification = "Moderately selective"
        else:
            classification = "Non-selective"
        
        return {
            "selectivity_index": round(selectivity_index, 1),
            "selectivity_window": f"{round(selectivity_index, 0)}-fold",
            "classification": classification,
            "closest_off_target_ic50": f"{min_off_target:.2f} µM"
        }
    
    def generate_safety_report(self, prediction: Dict) -> Dict:
        """Generate comprehensive safety assessment report"""
        report = {
            "compound": prediction["smiles"],
            "overall_safety": prediction["safety_classification"],
            "risk_score": prediction["risk_score"],
            "key_findings": [],
            "recommendations": [],
            "preclinical_strategy": []
        }
        
        # Summarize key findings
        if prediction["off_target_hits"]:
            critical_hits = [h for h in prediction["off_target_hits"] 
                           if h["severity"] == "critical"]
            if critical_hits:
                report["key_findings"].append(
                    f"Critical off-targets identified: {len(critical_hits)}"
                )
        
        if prediction["predicted_side_effects"]:
            report["key_findings"].append(
                f"Predicted side effects: {', '.join(prediction['predicted_side_effects'][:5])}"
            )
        
        # Generate recommendations
        if prediction["risk_score"] > 70:
            report["recommendations"].append("Consider major structural redesign")
            report["recommendations"].append("Focus on improving selectivity")
        elif prediction["risk_score"] > 40:
            report["recommendations"].append("Optimize selectivity against key off-targets")
            report["recommendations"].append("Consider biased agonism strategies")
        else:
            report["recommendations"].append("Proceed with standard safety studies")
        
        # Suggest preclinical strategy
        if "hERG" in [h["target"] for h in prediction["off_target_hits"]]:
            report["preclinical_strategy"].append("Early hERG patch clamp assay")
            report["preclinical_strategy"].append("Telemetry studies in conscious animals")
        
        if "5-HT2B" in [h["target"] for h in prediction["off_target_hits"]]:
            report["preclinical_strategy"].append("5-HT2B functional assay")
            report["preclinical_strategy"].append("Chronic valvulopathy studies")
        
        if prediction["drug_drug_interaction_risk"] in ["High", "Moderate"]:
            report["preclinical_strategy"].append("CYP inhibition panel")
            report["preclinical_strategy"].append("Drug-drug interaction studies")
        
        return report