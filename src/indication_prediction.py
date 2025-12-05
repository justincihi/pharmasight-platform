#!/usr/bin/env python3
"""
Indication Prediction System
Predicts therapeutic indications based on receptor binding profiles
Phase 6.1 of Advanced Cheminformatics Integration Plan
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdFingerprintGenerator
from typing import Dict, List, Tuple, Optional
import json

class IndicationPredictor:
    """Predict therapeutic indications from receptor binding profiles"""
    
    def __init__(self):
        self.receptor_indication_map = self._load_receptor_indication_map()
        self.indication_weights = self._load_indication_weights()
        self.therapeutic_classes = self._load_therapeutic_classes()
        self.dosage_reference = self._load_dosage_reference()
        
    def _load_receptor_indication_map(self) -> Dict:
        """Map receptor binding patterns to therapeutic indications"""
        return {
            "5-HT2A": {
                "agonist": ["Psychedelic-assisted therapy", "Treatment-resistant depression"],
                "antagonist": ["Antipsychotic", "Sleep disorders", "Migraine prevention"],
                "weight": 0.9
            },
            "5-HT1A": {
                "agonist": ["Anxiolytic", "Antidepressant", "Neuroprotective"],
                "antagonist": ["Cognitive enhancement"],
                "weight": 0.85
            },
            "5-HT2C": {
                "agonist": ["Weight loss", "Anti-addiction"],
                "antagonist": ["Antidepressant augmentation", "Anxiolytic"],
                "weight": 0.7
            },
            "5-HT3": {
                "antagonist": ["Antiemetic", "IBS treatment", "Anxiolytic"],
                "weight": 0.75
            },
            "5-HT7": {
                "antagonist": ["Antidepressant", "Circadian rhythm disorders", "Cognitive enhancement"],
                "weight": 0.6
            },
            "D1": {
                "agonist": ["Parkinson's disease", "Cognitive enhancement", "ADHD"],
                "antagonist": ["Antipsychotic"],
                "weight": 0.8
            },
            "D2": {
                "agonist": ["Parkinson's disease", "Restless leg syndrome", "Prolactinoma"],
                "antagonist": ["Antipsychotic", "Antiemetic", "Tourette syndrome"],
                "weight": 0.95
            },
            "D3": {
                "agonist": ["Parkinson's disease", "Restless leg syndrome"],
                "antagonist": ["Antipsychotic", "Anti-addiction"],
                "weight": 0.7
            },
            "D4": {
                "antagonist": ["ADHD", "Antipsychotic"],
                "weight": 0.55
            },
            "GABA-A": {
                "positive_modulator": ["Anxiolytic", "Sedative/Hypnotic", "Anticonvulsant", "Muscle relaxant"],
                "weight": 0.9
            },
            "GABA-B": {
                "agonist": ["Muscle relaxant", "Anti-spasticity", "Alcohol dependence"],
                "weight": 0.7
            },
            "NMDA": {
                "antagonist": ["Antidepressant", "Analgesic", "Neuroprotective", "Anesthetic"],
                "weight": 0.85
            },
            "AMPA": {
                "positive_modulator": ["Cognitive enhancement", "Antidepressant"],
                "weight": 0.6
            },
            "Mu_Opioid": {
                "agonist": ["Analgesic", "Antitussive", "Antidiarrheal"],
                "partial_agonist": ["Opioid use disorder treatment", "Analgesic"],
                "antagonist": ["Opioid overdose reversal", "Alcohol dependence"],
                "weight": 0.95
            },
            "Kappa_Opioid": {
                "agonist": ["Analgesic", "Antipruritic"],
                "antagonist": ["Antidepressant", "Anti-addiction"],
                "weight": 0.7
            },
            "Delta_Opioid": {
                "agonist": ["Analgesic", "Antidepressant", "Cardioprotective"],
                "weight": 0.6
            },
            "CB1": {
                "agonist": ["Analgesic", "Antiemetic", "Appetite stimulant"],
                "antagonist": ["Weight loss", "Anti-addiction"],
                "weight": 0.75
            },
            "CB2": {
                "agonist": ["Anti-inflammatory", "Analgesic", "Neuroprotective"],
                "weight": 0.65
            },
            "Alpha2_Adrenergic": {
                "agonist": ["Antihypertensive", "Sedative", "ADHD", "Opioid withdrawal"],
                "antagonist": ["Antidepressant augmentation", "Erectile dysfunction"],
                "weight": 0.7
            },
            "Beta1_Adrenergic": {
                "antagonist": ["Antihypertensive", "Antianginal", "Antiarrhythmic", "Heart failure"],
                "weight": 0.85
            },
            "Beta2_Adrenergic": {
                "agonist": ["Bronchodilator", "Asthma", "COPD", "Tocolytic"],
                "weight": 0.8
            },
            "H1": {
                "antagonist": ["Antihistamine", "Antiallergic", "Antiemetic", "Sedative"],
                "weight": 0.7
            },
            "H2": {
                "antagonist": ["Antiulcer", "GERD treatment"],
                "weight": 0.75
            },
            "H3": {
                "antagonist": ["Cognitive enhancement", "Narcolepsy", "ADHD"],
                "weight": 0.6
            },
            "M1": {
                "agonist": ["Cognitive enhancement", "Alzheimer's disease"],
                "weight": 0.7
            },
            "M2": {
                "antagonist": ["Cognitive enhancement", "Antidepressant augmentation"],
                "weight": 0.5
            },
            "M3": {
                "antagonist": ["Bronchodilator", "Overactive bladder", "COPD"],
                "weight": 0.65
            },
            "Sigma1": {
                "agonist": ["Neuroprotective", "Antidepressant", "Analgesic"],
                "weight": 0.55
            },
            "NET": {
                "inhibitor": ["Antidepressant", "ADHD", "Fibromyalgia", "Neuropathic pain"],
                "weight": 0.85
            },
            "DAT": {
                "inhibitor": ["ADHD", "Depression", "Narcolepsy", "Parkinson's"],
                "weight": 0.8
            },
            "SERT": {
                "inhibitor": ["Antidepressant", "Anxiety disorders", "OCD", "PTSD", "Premature ejaculation"],
                "weight": 0.95
            },
            "MAO-A": {
                "inhibitor": ["Antidepressant", "Anxiety disorders"],
                "weight": 0.7
            },
            "MAO-B": {
                "inhibitor": ["Parkinson's disease", "Neuroprotective"],
                "weight": 0.75
            },
            "AChE": {
                "inhibitor": ["Alzheimer's disease", "Myasthenia gravis", "Cognitive enhancement"],
                "weight": 0.85
            },
            "VGSC": {
                "blocker": ["Anticonvulsant", "Neuropathic pain", "Local anesthetic", "Antiarrhythmic"],
                "weight": 0.8
            },
            "VGCC_L": {
                "blocker": ["Antihypertensive", "Antianginal", "Antiarrhythmic"],
                "weight": 0.75
            },
            "VGCC_T": {
                "blocker": ["Anticonvulsant", "Analgesic"],
                "weight": 0.6
            },
            "hERG": {
                "caution": ["Cardiac safety concern - avoid binding"],
                "weight": -1.0
            },
            "COX-1": {
                "inhibitor": ["Antiplatelet", "Anti-inflammatory"],
                "weight": 0.6
            },
            "COX-2": {
                "inhibitor": ["Anti-inflammatory", "Analgesic", "Antipyretic"],
                "weight": 0.8
            },
            "PDE5": {
                "inhibitor": ["Erectile dysfunction", "Pulmonary hypertension"],
                "weight": 0.85
            },
            "PPAR-gamma": {
                "agonist": ["Type 2 diabetes", "Metabolic syndrome"],
                "weight": 0.75
            }
        }
    
    def _load_indication_weights(self) -> Dict:
        """Define weights for different indication categories"""
        return {
            "CNS_Depression": ["Antidepressant", "Anxiolytic", "Antipsychotic"],
            "CNS_Stimulation": ["ADHD", "Cognitive enhancement", "Narcolepsy"],
            "Pain": ["Analgesic", "Neuropathic pain", "Fibromyalgia", "Migraine prevention"],
            "Neurological": ["Parkinson's disease", "Alzheimer's disease", "Anticonvulsant", "Neuroprotective"],
            "Cardiovascular": ["Antihypertensive", "Antianginal", "Antiarrhythmic", "Heart failure"],
            "Metabolic": ["Type 2 diabetes", "Weight loss", "Metabolic syndrome"],
            "Respiratory": ["Asthma", "COPD", "Bronchodilator"],
            "GI": ["Antiemetic", "GERD treatment", "IBS treatment", "Antiulcer"],
            "Addiction": ["Anti-addiction", "Opioid use disorder treatment", "Alcohol dependence"],
            "Sleep": ["Sedative/Hypnotic", "Sleep disorders", "Circadian rhythm disorders"]
        }
    
    def _load_therapeutic_classes(self) -> Dict:
        """Define therapeutic drug classes with typical receptor profiles"""
        return {
            "Typical_Antipsychotic": {
                "primary": ["D2"],
                "secondary": ["5-HT2A", "H1", "Alpha1"],
                "indications": ["Schizophrenia", "Bipolar disorder", "Acute psychosis"]
            },
            "Atypical_Antipsychotic": {
                "primary": ["5-HT2A", "D2"],
                "secondary": ["5-HT1A", "Alpha1", "H1", "M1"],
                "indications": ["Schizophrenia", "Bipolar disorder", "Treatment-resistant depression"]
            },
            "SSRI": {
                "primary": ["SERT"],
                "secondary": [],
                "indications": ["Major depression", "Anxiety disorders", "OCD", "PTSD"]
            },
            "SNRI": {
                "primary": ["SERT", "NET"],
                "secondary": [],
                "indications": ["Major depression", "Anxiety disorders", "Fibromyalgia", "Neuropathic pain"]
            },
            "TCA": {
                "primary": ["SERT", "NET"],
                "secondary": ["H1", "M1", "Alpha1"],
                "indications": ["Major depression", "Neuropathic pain", "Enuresis"]
            },
            "Benzodiazepine": {
                "primary": ["GABA-A"],
                "secondary": [],
                "indications": ["Anxiety disorders", "Insomnia", "Seizures", "Muscle spasm"]
            },
            "Opioid_Analgesic": {
                "primary": ["Mu_Opioid"],
                "secondary": ["Kappa_Opioid", "Delta_Opioid"],
                "indications": ["Moderate-severe pain", "Chronic pain"]
            },
            "Psychedelic": {
                "primary": ["5-HT2A"],
                "secondary": ["5-HT2C", "5-HT1A"],
                "indications": ["Treatment-resistant depression", "PTSD", "End-of-life anxiety"]
            },
            "Dissociative": {
                "primary": ["NMDA"],
                "secondary": ["Sigma1", "Mu_Opioid"],
                "indications": ["Treatment-resistant depression", "Anesthesia", "Chronic pain"]
            },
            "Stimulant": {
                "primary": ["DAT", "NET"],
                "secondary": ["SERT"],
                "indications": ["ADHD", "Narcolepsy"]
            }
        }
    
    def _load_dosage_reference(self) -> Dict:
        """Reference dosage ranges based on receptor affinity patterns"""
        return {
            "ultra_potent": {
                "ki_range": (0.001, 0.1),
                "typical_dose_mg": (0.01, 1),
                "examples": ["Fentanyl", "LSD", "Carfentanil"]
            },
            "high_potency": {
                "ki_range": (0.1, 1),
                "typical_dose_mg": (0.5, 10),
                "examples": ["Haloperidol", "Risperidone", "Clonazepam"]
            },
            "moderate_potency": {
                "ki_range": (1, 10),
                "typical_dose_mg": (5, 50),
                "examples": ["Olanzapine", "Sertraline", "Diazepam"]
            },
            "low_potency": {
                "ki_range": (10, 100),
                "typical_dose_mg": (25, 200),
                "examples": ["Chlorpromazine", "Trazodone", "Quetiapine"]
            },
            "very_low_potency": {
                "ki_range": (100, 1000),
                "typical_dose_mg": (100, 1000),
                "examples": ["Buspirone", "Gabapentin"]
            }
        }
    
    def predict_indications(self, smiles: str, receptor_profile: Optional[Dict] = None,
                            action_modes: Optional[Dict] = None) -> Dict:
        """
        Predict therapeutic indications from SMILES and/or receptor binding profile
        
        Args:
            smiles: SMILES string of compound
            receptor_profile: Optional dict of receptor affinities {receptor: affinity_nM}
            action_modes: Optional dict specifying action at each receptor 
                         {receptor: "agonist"|"antagonist"|"modulator"|"inhibitor"}
        
        Returns:
            Dict with predicted indications, confidence scores, and therapeutic class
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        if receptor_profile is None:
            receptor_profile = self._predict_receptor_profile(mol)
        
        if action_modes is None:
            action_modes = self._infer_action_modes(mol, receptor_profile)
        
        indication_scores = {}
        indication_sources = {}
        therapeutic_class_scores = {}
        
        for receptor, affinity in receptor_profile.items():
            if receptor not in self.receptor_indication_map:
                continue
                
            receptor_data = self.receptor_indication_map[receptor]
            weight = receptor_data.get("weight", 0.5)
            
            if weight < 0:
                continue
            
            affinity_score = self._affinity_to_score(affinity)
            if affinity_score < 0.2:
                continue
            
            action = action_modes.get(receptor, self._guess_primary_action(receptor))
            
            relevant_actions = self._get_relevant_actions(action)
            
            for action_type in relevant_actions:
                if action_type in receptor_data:
                    for indication in receptor_data[action_type]:
                        score = affinity_score * weight * 0.8
                        
                        if indication in indication_scores:
                            if score > indication_scores[indication]:
                                indication_scores[indication] = score
                                indication_sources[indication] = {
                                    "receptor": receptor,
                                    "action": action_type,
                                    "affinity_nM": affinity
                                }
                        else:
                            indication_scores[indication] = score
                            indication_sources[indication] = {
                                "receptor": receptor,
                                "action": action_type,
                                "affinity_nM": affinity
                            }
        
        for class_name, class_data in self.therapeutic_classes.items():
            primary_receptors = class_data["primary"]
            secondary_receptors = class_data["secondary"]
            
            primary_hits = sum(
                1 for r in primary_receptors 
                if r in receptor_profile and receptor_profile[r] < 1000
            )
            primary_score = primary_hits / max(len(primary_receptors), 1)
            
            if primary_score < 0.5:
                therapeutic_class_scores[class_name] = 0.0
                continue
            
            primary_affinity_score = sum(
                self._affinity_to_score(receptor_profile.get(r, 10000))
                for r in primary_receptors
            ) / max(len(primary_receptors), 1)
            
            secondary_score = 0
            if secondary_receptors:
                secondary_hits = sum(
                    1 for r in secondary_receptors 
                    if r in receptor_profile and receptor_profile[r] < 5000
                )
                secondary_score = (secondary_hits / len(secondary_receptors)) * 0.3
            
            total_score = min((primary_affinity_score * 0.7 + secondary_score), 0.85)
            therapeutic_class_scores[class_name] = total_score
        
        sorted_indications = sorted(
            indication_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )[:15]
        
        sorted_classes = sorted(
            therapeutic_class_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )
        valid_classes = [(c, s) for c, s in sorted_classes if s > 0.3][:3]
        
        best_class = valid_classes[0][0] if valid_classes else "Unclassified"
        class_indications = self.therapeutic_classes.get(best_class, {}).get("indications", [])
        
        primary_indications = []
        secondary_indications = []
        
        for indication, score in sorted_indications:
            if score < 0.25:
                continue
                
            confidence = min(score * 85, 85)
            source = indication_sources.get(indication, {})
            
            ind_data = {
                "indication": indication,
                "confidence": round(confidence, 1),
                "evidence_level": self._score_to_evidence(score),
                "primary_receptor": source.get("receptor", "Unknown"),
                "mechanism": source.get("action", "Unknown")
            }
            
            if indication in class_indications and confidence > 40:
                primary_indications.append(ind_data)
            elif confidence > 35:
                secondary_indications.append(ind_data)
        
        dosage_estimate = self._estimate_dosage(receptor_profile, best_class)
        
        return {
            "smiles": smiles,
            "predicted_therapeutic_class": {
                "primary": best_class.replace("_", " ") if valid_classes else "Unclassified",
                "confidence": round(valid_classes[0][1] * 100, 1) if valid_classes else 0,
                "alternatives": [
                    {"class": c[0].replace("_", " "), "confidence": round(c[1] * 100, 1)}
                    for c in valid_classes[1:3]
                ]
            },
            "primary_indications": primary_indications[:5],
            "secondary_indications": secondary_indications[:5],
            "receptor_profile_used": receptor_profile,
            "action_modes_used": action_modes,
            "dosage_estimate": dosage_estimate,
            "confidence_note": "Predictions are computational estimates. Clinical validation required.",
            "development_considerations": self._get_development_considerations(
                sorted_indications, receptor_profile
            )
        }
    
    def _infer_action_modes(self, mol, receptor_profile: Dict) -> Dict:
        """Infer likely action mode (agonist/antagonist) based on structure"""
        action_modes = {}
        
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.ExactMolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        
        for receptor in receptor_profile.keys():
            if receptor in ["SERT", "NET", "DAT", "MAO-A", "MAO-B", "AChE"]:
                action_modes[receptor] = "inhibitor"
            elif receptor in ["VGSC", "VGCC_L", "VGCC_T", "hERG"]:
                action_modes[receptor] = "blocker"
            elif receptor in ["GABA-A", "AMPA"]:
                action_modes[receptor] = "positive_modulator"
            elif receptor in ["COX-1", "COX-2", "PDE5"]:
                action_modes[receptor] = "inhibitor"
            elif receptor.startswith("CYP"):
                action_modes[receptor] = "inhibitor"
            elif receptor in ["5-HT2A", "5-HT2C"]:
                if mw < 350 and hbd >= 2:
                    action_modes[receptor] = "agonist"
                else:
                    action_modes[receptor] = "antagonist"
            elif receptor.startswith("D"):
                if logp > 3.5 and mw > 350:
                    action_modes[receptor] = "antagonist"
                else:
                    action_modes[receptor] = "agonist"
            elif receptor.startswith("5-HT"):
                action_modes[receptor] = "agonist"
            elif receptor.startswith("Alpha") or receptor.startswith("Beta"):
                action_modes[receptor] = "antagonist"
            elif receptor.startswith("M"):
                action_modes[receptor] = "antagonist"
            elif receptor.startswith("H"):
                action_modes[receptor] = "antagonist"
            elif "Opioid" in receptor:
                action_modes[receptor] = "agonist"
            elif receptor.startswith("CB"):
                action_modes[receptor] = "agonist"
            else:
                action_modes[receptor] = "unknown"
        
        return action_modes
    
    def _guess_primary_action(self, receptor: str) -> str:
        """Guess the most common therapeutic action for a receptor"""
        antagonist_dominated = ["D2", "5-HT2A", "H1", "H2", "M1", "M3", "Alpha1A", 
                               "Beta1_Adrenergic", "NMDA"]
        inhibitor_types = ["SERT", "NET", "DAT", "MAO-A", "MAO-B", "AChE", 
                          "COX-1", "COX-2", "PDE5"]
        blocker_types = ["VGSC", "VGCC_L", "VGCC_T", "hERG"]
        
        if receptor in antagonist_dominated:
            return "antagonist"
        elif receptor in inhibitor_types:
            return "inhibitor"
        elif receptor in blocker_types:
            return "blocker"
        elif receptor == "GABA-A":
            return "positive_modulator"
        else:
            return "agonist"
    
    def _get_relevant_actions(self, action: str) -> List[str]:
        """Get list of action types to check based on inferred action"""
        action_mapping = {
            "agonist": ["agonist", "partial_agonist"],
            "antagonist": ["antagonist"],
            "partial_agonist": ["partial_agonist", "agonist"],
            "positive_modulator": ["positive_modulator", "agonist"],
            "inhibitor": ["inhibitor", "antagonist"],
            "blocker": ["blocker", "antagonist", "inhibitor"],
            "unknown": ["agonist", "antagonist"]
        }
        return action_mapping.get(action, ["agonist", "antagonist"])
    
    def _predict_receptor_profile(self, mol) -> Dict:
        """Predict receptor binding profile from molecular structure"""
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = fpgen.GetFingerprint(mol)
        
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        n_atoms = mol.GetNumAtoms()
        n_heteroatoms = Descriptors.NumHeteroatoms(mol)
        
        profile = {}
        
        if 250 < mw < 500 and 2 < logp < 5 and aromatic_rings >= 2:
            profile["5-HT2A"] = max(10, 500 - mw + (logp * 50))
            profile["5-HT2C"] = profile["5-HT2A"] * 1.5
            profile["D2"] = max(10, 600 - mw + (logp * 30))
        
        if hbd <= 3 and hba <= 7 and tpsa < 90:
            profile["GABA-A"] = max(50, 300 - tpsa * 2)
        
        if logp > 3 and aromatic_rings >= 1:
            profile["Mu_Opioid"] = max(20, 100 + logp * 50)
            profile["Kappa_Opioid"] = profile["Mu_Opioid"] * 2
        
        if n_heteroatoms >= 3 and tpsa > 40:
            profile["SERT"] = max(5, 200 - n_heteroatoms * 20)
            profile["NET"] = profile["SERT"] * 1.5
        
        if aromatic_rings >= 2 and logp > 2:
            profile["H1"] = max(100, 500 - aromatic_rings * 100)
            profile["M1"] = max(200, 800 - aromatic_rings * 100)
        
        if mw < 300 and logp < 3:
            profile["NMDA"] = max(100, 500 - mw)
        
        return profile
    
    def _affinity_to_score(self, affinity_nM: float) -> float:
        """Convert binding affinity (nM) to a 0-1 score"""
        if affinity_nM <= 0:
            return 0
        log_affinity = np.log10(affinity_nM)
        score = 1 - (log_affinity / 4)
        return max(0, min(1, score))
    
    def _score_to_evidence(self, score: float) -> str:
        """Convert score to evidence level description"""
        if score > 0.8:
            return "Strong - High affinity for key targets"
        elif score > 0.6:
            return "Moderate - Significant target engagement expected"
        elif score > 0.4:
            return "Suggestive - Some relevant target activity"
        else:
            return "Weak - Limited target engagement"
    
    def _estimate_dosage(self, receptor_profile: Dict, therapeutic_class: Optional[str]) -> Dict:
        """Estimate therapeutic dosage range based on receptor affinity with uncertainty bounds"""
        if not receptor_profile:
            return {
                "error": "No receptor profile available",
                "confidence": "None",
                "uncertainty": "Cannot estimate without receptor data"
            }
        
        relevant_affinities = [v for v in receptor_profile.values() if v < 10000]
        if not relevant_affinities:
            return {
                "potency_class": "Insufficient Binding",
                "estimated_dose_range": "Not estimable",
                "confidence": "Very Low",
                "uncertainty": "No significant receptor binding detected",
                "notes": ["Compound may lack therapeutic activity at typical doses"]
            }
        
        min_affinity = min(relevant_affinities)
        median_affinity = sorted(relevant_affinities)[len(relevant_affinities)//2]
        
        for potency_class, data in self.dosage_reference.items():
            ki_min, ki_max = data["ki_range"]
            if ki_min <= min_affinity < ki_max:
                dose_min, dose_max = data["typical_dose_mg"]
                
                uncertainty_factor = 1 + (len(receptor_profile) * 0.1)
                adjusted_min = max(dose_min * 0.5, dose_min / uncertainty_factor)
                adjusted_max = min(dose_max * 2, dose_max * uncertainty_factor)
                
                confidence_level = "Moderate" if min_affinity < 10 else "Low"
                if therapeutic_class in ["SSRI", "Benzodiazepine", "Opioid_Analgesic"]:
                    confidence_level = "Moderate-High" if min_affinity < 10 else "Moderate"
                
                return {
                    "potency_class": potency_class.replace("_", " ").title(),
                    "estimated_dose_range": f"{dose_min} - {dose_max} mg",
                    "dose_min_mg": dose_min,
                    "dose_max_mg": dose_max,
                    "uncertainty_range": f"{adjusted_min:.1f} - {adjusted_max:.1f} mg",
                    "reference_compounds": data["examples"],
                    "confidence": confidence_level,
                    "based_on_affinity_nM": round(min_affinity, 2),
                    "uncertainty_notes": [
                        "Estimate based on receptor affinity correlation only",
                        "Actual dosing depends on ADME properties, therapeutic index",
                        "Bioavailability and metabolism significantly affect dosing",
                        "Clinical validation required before human use"
                    ],
                    "factors_not_considered": [
                        "Oral bioavailability",
                        "First-pass metabolism",
                        "Protein binding",
                        "Therapeutic window"
                    ]
                }
        
        return {
            "potency_class": "Very Low Potency",
            "estimated_dose_range": "100 - 1000+ mg",
            "dose_min_mg": 100,
            "dose_max_mg": 1000,
            "uncertainty_range": "50 - 2000 mg",
            "reference_compounds": ["High-dose therapeutics"],
            "confidence": "Low",
            "based_on_affinity_nM": round(min_affinity, 2),
            "uncertainty_notes": [
                "Weak receptor binding suggests high doses may be needed",
                "Consider if compound has sufficient potency for practical use",
                "Alternative targets or mechanisms may be more relevant"
            ]
        }
    
    def _get_development_considerations(self, indications: List[Tuple], 
                                        receptor_profile: Dict) -> List[Dict]:
        """Generate development considerations based on predicted profile"""
        considerations = []
        
        indication_names = [ind[0] for ind in indications[:5]]
        
        if any("depression" in i.lower() or "antidepressant" in i.lower() 
               for i in indication_names):
            considerations.append({
                "area": "Regulatory",
                "consideration": "FDA requires suicidality assessment for antidepressants",
                "priority": "High"
            })
        
        if any("psychotic" in i.lower() or "schizophrenia" in i.lower() 
               for i in indication_names):
            considerations.append({
                "area": "Safety",
                "consideration": "Monitor for metabolic syndrome and extrapyramidal symptoms",
                "priority": "High"
            })
        
        if "D2" in receptor_profile and receptor_profile["D2"] < 10:
            considerations.append({
                "area": "Safety",
                "consideration": "Strong D2 binding may cause hyperprolactinemia and EPS",
                "priority": "High"
            })
        
        if "Mu_Opioid" in receptor_profile:
            considerations.append({
                "area": "Regulatory",
                "consideration": "Opioid activity requires DEA scheduling consideration",
                "priority": "Critical"
            })
            considerations.append({
                "area": "Safety",
                "consideration": "Assess abuse liability and respiratory depression risk",
                "priority": "Critical"
            })
        
        if "5-HT2A" in receptor_profile and receptor_profile["5-HT2A"] < 100:
            considerations.append({
                "area": "Clinical",
                "consideration": "Potential for psychedelic effects - consider controlled setting",
                "priority": "High"
            })
        
        if any("ADHD" in i or "Narcolepsy" in i for i in indication_names):
            considerations.append({
                "area": "Regulatory",
                "consideration": "Stimulant activity may require controlled substance scheduling",
                "priority": "High"
            })
        
        if any("pain" in i.lower() or "analgesic" in i.lower() for i in indication_names):
            considerations.append({
                "area": "Clinical",
                "consideration": "Include validated pain assessment scales in trials",
                "priority": "Moderate"
            })
        
        if not considerations:
            considerations.append({
                "area": "General",
                "consideration": "Standard IND-enabling studies recommended",
                "priority": "Moderate"
            })
        
        return considerations
    
    def batch_predict(self, compounds: List[Dict]) -> List[Dict]:
        """
        Predict indications for multiple compounds
        
        Args:
            compounds: List of dicts with 'smiles' and optional 'receptor_profile'
        
        Returns:
            List of prediction results
        """
        results = []
        for compound in compounds[:50]:
            smiles = compound.get("smiles", "")
            receptor_profile = compound.get("receptor_profile")
            
            result = self.predict_indications(smiles, receptor_profile)
            result["compound_id"] = compound.get("id", "")
            result["compound_name"] = compound.get("name", "")
            results.append(result)
        
        return results


predictor = IndicationPredictor()


def predict_indications(smiles: str, receptor_profile: Optional[Dict] = None) -> Dict:
    """Convenience function for indication prediction"""
    return predictor.predict_indications(smiles, receptor_profile)


def batch_predict_indications(compounds: List[Dict]) -> List[Dict]:
    """Convenience function for batch prediction"""
    return predictor.batch_predict(compounds)


def estimate_dosage(receptor_profile: Dict, therapeutic_class: Optional[str] = None) -> Dict:
    """Convenience function for dosage estimation"""
    return predictor._estimate_dosage(receptor_profile, therapeutic_class)
