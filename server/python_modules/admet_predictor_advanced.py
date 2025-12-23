#!/usr/bin/env python3
"""
Advanced ADMET Prediction Module
Predicts Absorption, Distribution, Metabolism, Excretion, and Toxicity
Including BBB permeability, SimCyp, and NONMEM-style predictions
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from typing import Dict, List, Optional
import json

class ADMETPredictor:
    """Comprehensive ADMET property prediction"""
    
    def __init__(self):
        self.load_prediction_models()
        
    def load_prediction_models(self):
        """Load or initialize prediction models"""
        # In production, these would be trained ML models
        # For now, using rule-based predictions
        pass
    
    def predict_all_properties(self, smiles: str) -> Dict:
        """Predict all ADMET properties for a compound"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        # Calculate molecular descriptors
        descriptors = self.calculate_descriptors(mol)
        
        # Make predictions
        predictions = {
            "absorption": self.predict_absorption(mol, descriptors),
            "distribution": self.predict_distribution(mol, descriptors),
            "metabolism": self.predict_metabolism(mol, descriptors),
            "excretion": self.predict_excretion(mol, descriptors),
            "toxicity": self.predict_toxicity(mol, descriptors),
            "bbb_permeability": self.predict_bbb_permeability(mol, descriptors),
            "cyp_interactions": self.predict_cyp_interactions(mol),
            "pkpd_parameters": self.predict_pkpd_parameters(mol, descriptors),
            "clinical_trials_readiness": self.assess_clinical_readiness(mol, descriptors)
        }
        
        # Add overall drug-likeness score
        predictions["overall_score"] = self.calculate_overall_score(predictions)
        predictions["recommendations"] = self.generate_recommendations(predictions)
        
        return predictions
    
    def calculate_descriptors(self, mol) -> Dict:
        """Calculate molecular descriptors"""
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "logd": Crippen.MolLogP(mol) - 0.5,  # Simplified LogD estimation
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "heavy_atoms": Descriptors.HeavyAtomCount(mol),
            "molar_refractivity": Crippen.MolMR(mol),
            "formal_charge": Chem.rdmolops.GetFormalCharge(mol),
            "num_rings": Descriptors.RingCount(mol),
            "fraction_sp3": Descriptors.FractionCSP3(mol)
        }
    
    def predict_absorption(self, mol, descriptors: Dict) -> Dict:
        """Predict absorption properties"""
        # Human Intestinal Absorption (HIA)
        hia_score = 100
        if descriptors["mw"] > 500:
            hia_score -= 20
        if descriptors["logp"] < -0.4 or descriptors["logp"] > 5.6:
            hia_score -= 15
        if descriptors["hbd"] > 5:
            hia_score -= 10
        if descriptors["tpsa"] > 140:
            hia_score -= 25
        
        # Caco-2 permeability (intestinal epithelium)
        if descriptors["tpsa"] < 60 and descriptors["logp"] > 1:
            caco2 = "High (>8×10⁻⁶ cm/s)"
        elif descriptors["tpsa"] < 90:
            caco2 = "Moderate (1-8×10⁻⁶ cm/s)"
        else:
            caco2 = "Low (<1×10⁻⁶ cm/s)"
        
        # P-glycoprotein substrate/inhibitor
        pgp_substrate = descriptors["mw"] > 400 and descriptors["hbd"] > 2
        pgp_inhibitor = descriptors["logp"] > 3.5 and descriptors["mw"] > 300
        
        # Bioavailability prediction
        f_score = 100
        violations = 0
        if descriptors["mw"] > 500:
            violations += 1
        if descriptors["logp"] > 5:
            violations += 1
        if descriptors["hbd"] > 5:
            violations += 1
        if descriptors["hba"] > 10:
            violations += 1
        if descriptors["rotatable_bonds"] > 10:
            violations += 1
        if descriptors["tpsa"] > 140:
            violations += 1
        
        f_score = max(20, 100 - violations * 15)
        
        return {
            "hia_probability": f"{max(0, min(100, hia_score))}%",
            "caco2_permeability": caco2,
            "pgp_substrate": "Yes" if pgp_substrate else "No",
            "pgp_inhibitor": "Yes" if pgp_inhibitor else "No",
            "oral_bioavailability": f"{f_score}%",
            "solubility_class": self._predict_solubility_class(descriptors),
            "absorption_rate": "Fast" if hia_score > 70 else "Moderate" if hia_score > 40 else "Slow"
        }
    
    def predict_distribution(self, mol, descriptors: Dict) -> Dict:
        """Predict distribution properties"""
        # Volume of distribution (Vd) prediction
        if descriptors["logp"] > 3 and descriptors["mw"] < 600:
            vd = "High (>0.7 L/kg)"
            vd_value = round(0.7 + descriptors["logp"] * 0.2, 2)
        elif descriptors["logp"] > 1:
            vd = "Moderate (0.2-0.7 L/kg)"
            vd_value = round(0.2 + descriptors["logp"] * 0.15, 2)
        else:
            vd = "Low (<0.2 L/kg)"
            vd_value = round(max(0.04, 0.2 - abs(descriptors["logp"]) * 0.05), 2)
        
        # Plasma protein binding (PPB)
        ppb = min(99, max(20, 50 + descriptors["logp"] * 10))
        if descriptors["formal_charge"] != 0:
            ppb -= 10
        
        # Tissue distribution
        if descriptors["logp"] > 2:
            tissue_dist = "Extensive"
        elif descriptors["logp"] > 0:
            tissue_dist = "Moderate"
        else:
            tissue_dist = "Limited"
        
        return {
            "volume_distribution": vd,
            "vd_value": f"{vd_value} L/kg",
            "plasma_protein_binding": f"{round(ppb)}%",
            "tissue_distribution": tissue_dist,
            "fraction_unbound": f"{round(100-ppb)}%",
            "steady_state_conc": "Reached in 4-5 half-lives"
        }
    
    def predict_metabolism(self, mol, descriptors: Dict) -> Dict:
        """Predict metabolism properties"""
        # CYP substrate prediction (simplified)
        cyp_substrates = {
            "CYP3A4": descriptors["mw"] > 300 and descriptors["logp"] > 2,
            "CYP2D6": descriptors["mw"] < 500 and descriptors["hbd"] < 3,
            "CYP2C9": descriptors["mw"] < 400 and descriptors["formal_charge"] < 0,
            "CYP2C19": descriptors["mw"] < 450 and descriptors["aromatic_rings"] > 0,
            "CYP1A2": descriptors["aromatic_rings"] > 1 and descriptors["mw"] < 350
        }
        
        # Metabolic stability
        if descriptors["rotatable_bonds"] > 7:
            stability = "Low"
            t_half = "< 30 min"
        elif descriptors["rotatable_bonds"] > 4:
            stability = "Moderate"
            t_half = "30-120 min"
        else:
            stability = "High"
            t_half = "> 120 min"
        
        # Phase I/II metabolism
        phase1 = []
        phase2 = []
        
        if descriptors["aromatic_rings"] > 0:
            phase1.append("Hydroxylation")
        if descriptors["logp"] > 3:
            phase1.append("N-dealkylation")
        if "O" in Chem.MolToSmiles(mol):
            phase1.append("O-dealkylation")
        
        if descriptors["hbd"] > 0:
            phase2.append("Glucuronidation")
        if "N" in Chem.MolToSmiles(mol):
            phase2.append("Acetylation")
        if descriptors["aromatic_rings"] > 0:
            phase2.append("Sulfation")
        
        return {
            "cyp_substrates": [k for k, v in cyp_substrates.items() if v],
            "metabolic_stability": stability,
            "hepatic_clearance": self._predict_clearance(descriptors),
            "half_life": t_half,
            "phase1_reactions": phase1 if phase1 else ["None predicted"],
            "phase2_reactions": phase2 if phase2 else ["None predicted"],
            "first_pass_effect": "High" if descriptors["logp"] > 2 else "Low"
        }
    
    def predict_excretion(self, mol, descriptors: Dict) -> Dict:
        """Predict excretion properties"""
        # Renal clearance
        if descriptors["mw"] < 400 and descriptors["formal_charge"] != 0:
            renal_clearance = "High"
            clearance_route = "Primarily renal"
        elif descriptors["mw"] > 600:
            renal_clearance = "Low"
            clearance_route = "Primarily hepatic/biliary"
        else:
            renal_clearance = "Moderate"
            clearance_route = "Mixed renal/hepatic"
        
        # Total clearance estimation
        cl_total = 0.693 * vd_value * 1000 / 240  # Simplified calculation
        
        return {
            "renal_clearance": renal_clearance,
            "clearance_route": clearance_route,
            "total_clearance": f"{round(cl_total, 1)} mL/min/kg",
            "biliary_excretion": "Yes" if descriptors["mw"] > 500 else "No",
            "urine_recovery": "High (>60%)" if renal_clearance == "High" else "Low (<30%)"
        }
    
    def predict_toxicity(self, mol, descriptors: Dict) -> Dict:
        """Predict toxicity properties"""
        # hERG inhibition (cardiotoxicity)
        herg_risk = "Low"
        if descriptors["logp"] > 3.5 and descriptors["mw"] > 300:
            herg_risk = "Medium"
        if descriptors["logp"] > 4.5 and descriptors["hbd"] < 2:
            herg_risk = "High"
        
        # Hepatotoxicity
        hepatotox_risk = "Low"
        if descriptors["logp"] > 3 and descriptors["mw"] > 400:
            hepatotox_risk = "Medium"
        if descriptors["rotatable_bonds"] > 10:
            hepatotox_risk = "High"
        
        # Mutagenicity (simplified Ames test prediction)
        ames_positive = False
        toxic_substructs = ["N=N", "N(=O)O", "[N+](=O)[O-]"]
        smiles = Chem.MolToSmiles(mol)
        for substruct in toxic_substructs:
            if substruct in smiles:
                ames_positive = True
                break
        
        # LD50 prediction (very simplified)
        ld50_estimate = 500 + (5 - descriptors["logp"]) * 200
        ld50_estimate = max(50, min(5000, ld50_estimate))
        
        return {
            "herg_inhibition_risk": herg_risk,
            "hepatotoxicity_risk": hepatotox_risk,
            "ames_mutagenicity": "Positive" if ames_positive else "Negative",
            "ld50_estimate": f"{round(ld50_estimate)} mg/kg",
            "carcinogenicity_risk": "Low" if not ames_positive else "Medium",
            "skin_sensitization": "Low" if descriptors["tpsa"] < 100 else "Medium",
            "safety_margin": self._calculate_safety_margin(ld50_estimate)
        }
    
    def predict_bbb_permeability(self, mol, descriptors: Dict) -> Dict:
        """Predict blood-brain barrier permeability"""
        # BBB penetration score
        bbb_score = 0
        
        # Favorable factors
        if 1 < descriptors["logp"] < 4:
            bbb_score += 30
        if descriptors["mw"] < 450:
            bbb_score += 25
        if descriptors["tpsa"] < 90:
            bbb_score += 25
        if descriptors["hbd"] <= 5:
            bbb_score += 10
        if descriptors["hba"] <= 10:
            bbb_score += 10
        
        # Unfavorable factors
        if descriptors["formal_charge"] != 0:
            bbb_score -= 30
        if descriptors["tpsa"] > 120:
            bbb_score -= 40
        
        bbb_score = max(0, min(100, bbb_score))
        
        # CNS MPO score (CNS multiparameter optimization)
        cns_mpo = self._calculate_cns_mpo(descriptors)
        
        # P-gp efflux at BBB
        pgp_efflux = descriptors["mw"] > 400 and descriptors["logp"] > 2
        
        return {
            "bbb_penetration": "High" if bbb_score > 70 else "Moderate" if bbb_score > 40 else "Low",
            "bbb_score": f"{bbb_score}%",
            "cns_mpo_score": f"{cns_mpo}/6.0",
            "pgp_efflux_risk": "Yes" if pgp_efflux else "No",
            "brain_plasma_ratio": "High (>0.3)" if bbb_score > 70 else "Low (<0.1)",
            "cns_activity": "Likely" if bbb_score > 60 and cns_mpo > 4 else "Unlikely"
        }
    
    def predict_cyp_interactions(self, mol) -> Dict:
        """Predict CYP enzyme interactions"""
        smiles = Chem.MolToSmiles(mol)
        
        interactions = {
            "CYP3A4": {
                "substrate": "Likely" if Descriptors.MolWt(mol) > 300 else "Unlikely",
                "inhibitor": "Possible" if "c1ccccc1" in smiles else "Unlikely",
                "inducer": "Unlikely"
            },
            "CYP2D6": {
                "substrate": "Likely" if "N" in smiles else "Unlikely",
                "inhibitor": "Possible" if Crippen.MolLogP(mol) > 3 else "Unlikely",
                "inducer": "Unlikely"
            },
            "CYP2C9": {
                "substrate": "Possible",
                "inhibitor": "Unlikely",
                "inducer": "Unlikely"
            }
        }
        
        return interactions
    
    def predict_pkpd_parameters(self, mol, descriptors: Dict) -> Dict:
        """Predict PK/PD parameters for SimCyp/NONMEM modeling"""
        # Simplified PK parameters
        ka = 1.0 / (1 + np.exp(-(-2 + descriptors["logp"])))  # Absorption rate constant
        ke = 0.693 / 240  # Elimination rate constant (assume 4h half-life)
        vd_value = 0.2 + descriptors["logp"] * 0.15  # Volume of distribution
        cl = ke * vd_value * 70  # Clearance for 70kg person
        
        # Simplified PD parameters
        if descriptors["mw"] < 500:
            emax = 100
            ec50 = 10 ** (descriptors["logp"] - 2)
        else:
            emax = 80
            ec50 = 10 ** (descriptors["logp"] - 1)
        
        return {
            "absorption_rate_ka": f"{round(ka, 3)} h⁻¹",
            "elimination_rate_ke": f"{round(ke, 4)} h⁻¹",
            "volume_distribution": f"{round(vd_value, 2)} L/kg",
            "clearance": f"{round(cl, 1)} L/h",
            "half_life": f"{round(0.693/ke, 1)} h",
            "time_to_peak": f"{round(np.log(ka/ke)/(ka-ke), 1)} h",
            "steady_state_time": f"{round(5 * 0.693/ke, 0)} h",
            "emax": f"{emax}%",
            "ec50": f"{round(ec50, 2)} µM",
            "therapeutic_index": round(ld50_estimate / ec50, 1)
        }
    
    def assess_clinical_readiness(self, mol, descriptors: Dict) -> Dict:
        """Assess readiness for clinical trials"""
        score = 0
        issues = []
        
        # Check Lipinski's Rule of Five
        lipinski_violations = 0
        if descriptors["mw"] > 500:
            lipinski_violations += 1
            issues.append("MW > 500 Da")
        if descriptors["logp"] > 5:
            lipinski_violations += 1
            issues.append("LogP > 5")
        if descriptors["hbd"] > 5:
            lipinski_violations += 1
            issues.append("H-bond donors > 5")
        if descriptors["hba"] > 10:
            lipinski_violations += 1
            issues.append("H-bond acceptors > 10")
        
        if lipinski_violations == 0:
            score += 40
        elif lipinski_violations == 1:
            score += 25
        elif lipinski_violations == 2:
            score += 10
        
        # Check other drug-like properties
        if descriptors["tpsa"] < 140:
            score += 20
        else:
            issues.append("High TPSA (>140)")
        
        if descriptors["rotatable_bonds"] < 10:
            score += 20
        else:
            issues.append("Too flexible (>10 rotatable bonds)")
        
        if descriptors["aromatic_rings"] > 0:
            score += 10
        
        if descriptors["formal_charge"] == 0:
            score += 10
        else:
            issues.append("Charged molecule")
        
        # Determine clinical phase recommendation
        if score >= 80:
            phase = "Ready for Phase I"
            recommendation = "Proceed to IND-enabling studies"
        elif score >= 60:
            phase = "Lead optimization needed"
            recommendation = "Address identified issues before clinical studies"
        else:
            phase = "Early discovery"
            recommendation = "Significant optimization required"
        
        return {
            "clinical_readiness_score": f"{score}/100",
            "recommended_phase": phase,
            "recommendation": recommendation,
            "issues_to_address": issues if issues else ["None identified"],
            "regulatory_flags": self._check_regulatory_flags(mol)
        }
    
    def _predict_solubility_class(self, descriptors: Dict) -> str:
        """Predict BCS solubility class"""
        if descriptors["logp"] < 2 and descriptors["mw"] < 400:
            return "High (BCS Class I/III)"
        elif descriptors["logp"] > 4:
            return "Low (BCS Class II/IV)"
        else:
            return "Moderate"
    
    def _predict_clearance(self, descriptors: Dict) -> str:
        """Predict hepatic clearance"""
        if descriptors["logp"] > 3 and descriptors["mw"] > 400:
            return "High (>15 mL/min/kg)"
        elif descriptors["logp"] > 1:
            return "Moderate (5-15 mL/min/kg)"
        else:
            return "Low (<5 mL/min/kg)"
    
    def _calculate_safety_margin(self, ld50: float) -> str:
        """Calculate therapeutic safety margin"""
        therapeutic_dose = 10  # mg/kg (assumed)
        margin = ld50 / therapeutic_dose
        
        if margin > 100:
            return f"Excellent ({round(margin)}x)"
        elif margin > 10:
            return f"Good ({round(margin)}x)"
        else:
            return f"Narrow ({round(margin)}x)"
    
    def _calculate_cns_mpo(self, descriptors: Dict) -> float:
        """Calculate CNS multiparameter optimization score"""
        score = 0
        
        # Optimal ranges for CNS drugs
        if 2 <= descriptors["logp"] <= 4:
            score += 1
        if 200 <= descriptors["mw"] <= 360:
            score += 1
        if 40 <= descriptors["tpsa"] <= 90:
            score += 1
        if descriptors["hbd"] <= 2:
            score += 1
        if 2 <= descriptors["hba"] <= 5:
            score += 1
        if 1 <= descriptors["logp"] - descriptors.get("logd", descriptors["logp"]) <= 3:
            score += 1
        
        return score
    
    def _check_regulatory_flags(self, mol) -> List[str]:
        """Check for regulatory concerns"""
        flags = []
        smiles = Chem.MolToSmiles(mol)
        
        # Check for problematic substructures
        if "N=N" in smiles:
            flags.append("Contains azo group")
        if "[N+](=O)[O-]" in smiles:
            flags.append("Contains nitro group")
        if "S(=O)(=O)N" in smiles:
            flags.append("Contains sulfonamide")
        
        if not flags:
            flags.append("No major structural alerts")
        
        return flags
    
    def calculate_overall_score(self, predictions: Dict) -> float:
        """Calculate overall ADMET score"""
        score = 0
        weights = {
            "absorption": 20,
            "distribution": 15,
            "metabolism": 20,
            "excretion": 15,
            "toxicity": 30
        }
        
        # Score each category
        if "High" in str(predictions.get("absorption", {}).get("hia_probability", "")):
            score += weights["absorption"]
        elif "Moderate" in str(predictions.get("absorption", {})):
            score += weights["absorption"] * 0.5
        
        if "Moderate" in str(predictions.get("distribution", {}).get("volume_distribution", "")):
            score += weights["distribution"]
        
        if "Moderate" in str(predictions.get("metabolism", {}).get("metabolic_stability", "")):
            score += weights["metabolism"]
        elif "High" in str(predictions.get("metabolism", {}).get("metabolic_stability", "")):
            score += weights["metabolism"] * 0.7
        
        if "Mixed" in str(predictions.get("excretion", {}).get("clearance_route", "")):
            score += weights["excretion"]
        
        if "Low" in str(predictions.get("toxicity", {}).get("herg_inhibition_risk", "")):
            score += weights["toxicity"] * 0.5
        if "Negative" in str(predictions.get("toxicity", {}).get("ames_mutagenicity", "")):
            score += weights["toxicity"] * 0.5
        
        return round(score, 1)
    
    def generate_recommendations(self, predictions: Dict) -> List[str]:
        """Generate optimization recommendations"""
        recommendations = []
        
        # Check absorption
        if "Low" in str(predictions.get("absorption", {}).get("hia_probability", "")):
            recommendations.append("Improve absorption: reduce MW or TPSA")
        
        # Check distribution
        if "Low" in str(predictions.get("distribution", {}).get("volume_distribution", "")):
            recommendations.append("Increase lipophilicity for better tissue distribution")
        
        # Check metabolism
        if "Low" in str(predictions.get("metabolism", {}).get("metabolic_stability", "")):
            recommendations.append("Reduce metabolic liability: minimize rotatable bonds")
        
        # Check toxicity
        if "High" in str(predictions.get("toxicity", {}).get("herg_inhibition_risk", "")):
            recommendations.append("Reduce hERG risk: decrease LogP or add polar groups")
        
        # Check BBB
        if predictions.get("bbb_permeability", {}).get("cns_activity") == "Likely":
            recommendations.append("CNS penetration detected - consider if intended")
        
        if not recommendations:
            recommendations.append("Compound shows favorable ADMET profile")
            recommendations.append("Consider proceeding to in vitro validation")
        
        return recommendations