"""
Enhanced PKPD & DDI Analysis System
Integrates with open-source pharmacokinetic models for comprehensive drug interaction analysis
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Any, Tuple
import datetime

class EnhancedPKPDSystem:
    def __init__(self):
        self.drug_database = self.load_drug_database()
        self.interaction_matrix = self.load_interaction_matrix()
        self.pkpd_models = self.initialize_pkpd_models()
    
    def load_drug_database(self) -> Dict[str, Any]:
        """Load comprehensive drug database with PK/PD parameters"""
        return {
            "buprenorphine": {
                "molecular_weight": 467.64,
                "bioavailability": {"oral": 0.15, "sublingual": 0.55, "iv": 1.0},
                "protein_binding": 0.96,
                "vd": 430,  # L/kg
                "clearance": 1.2,  # L/h/kg
                "half_life": {"alpha": 0.5, "beta": 24},  # hours
                "metabolism": {
                    "primary": "CYP3A4",
                    "secondary": ["CYP2C8", "UGT1A1", "UGT2B7"],
                    "metabolites": ["norbuprenorphine", "buprenorphine-3-glucuronide"]
                },
                "receptors": {
                    "MOR": {"type": "partial_agonist", "ki": 0.2},
                    "KOR": {"type": "antagonist", "ki": 2.5},
                    "DOR": {"type": "antagonist", "ki": 6.1}
                },
                "therapeutic_range": {"min": 0.5, "max": 10},  # ng/mL
                "toxic_level": 50,
                "interactions": ["CYP3A4_inhibitors", "CNS_depressants"]
            },
            "ketamine": {
                "molecular_weight": 237.73,
                "bioavailability": {"oral": 0.17, "nasal": 0.45, "iv": 1.0, "im": 0.93},
                "protein_binding": 0.12,
                "vd": 3.1,  # L/kg
                "clearance": 0.89,  # L/h/kg
                "half_life": {"alpha": 0.5, "beta": 2.5},  # hours
                "metabolism": {
                    "primary": "CYP2B6",
                    "secondary": ["CYP3A4", "CYP2C9"],
                    "metabolites": ["norketamine", "dehydronorketamine"]
                },
                "receptors": {
                    "NMDA": {"type": "antagonist", "ki": 0.53},
                    "AMPA": {"type": "positive_modulator", "ki": 15},
                    "mTOR": {"type": "activator", "ec50": 2.1}
                },
                "therapeutic_range": {"min": 100, "max": 1000},  # ng/mL
                "neuroplasticity_window": 24,  # hours
                "interactions": ["CYP2B6_inhibitors", "NMDA_antagonists"]
            },
            "psilocybin": {
                "molecular_weight": 284.25,
                "bioavailability": {"oral": 0.72},
                "protein_binding": 0.05,
                "vd": 0.84,  # L/kg
                "clearance": 0.15,  # L/h/kg
                "half_life": {"psilocybin": 0.5, "psilocin": 2.5},  # hours
                "metabolism": {
                    "primary": "alkaline_phosphatase",
                    "secondary": ["MAO-A"],
                    "metabolites": ["psilocin", "4-hydroxyindole-3-acetic_acid"]
                },
                "receptors": {
                    "5HT2A": {"type": "agonist", "ki": 0.006},
                    "5HT2C": {"type": "agonist", "ki": 0.054},
                    "5HT1A": {"type": "agonist", "ki": 0.23}
                },
                "therapeutic_range": {"min": 5, "max": 50},  # ng/mL
                "neuroplasticity_window": 72,  # hours
                "interactions": ["MAO_inhibitors", "serotonergic_drugs"]
            },
            "tropacocaine": {
                "molecular_weight": 245.32,
                "bioavailability": {"oral": 0.65, "nasal": 0.85},
                "protein_binding": 0.45,
                "vd": 1.8,  # L/kg
                "clearance": 0.95,  # L/h/kg
                "half_life": {"alpha": 0.3, "beta": 1.2},  # hours
                "metabolism": {
                    "primary": "plasma_esterases",
                    "secondary": ["CYP3A4"],
                    "metabolites": ["benzoylecgonine", "ecgonine"]
                },
                "receptors": {
                    "DAT": {"type": "inhibitor", "ki": 450},
                    "NET": {"type": "inhibitor", "ki": 1200},
                    "SERT": {"type": "inhibitor", "ki": 2800}
                },
                "therapeutic_range": {"min": 50, "max": 300},  # ng/mL
                "cardiotoxicity_threshold": 800,  # ng/mL (vs 400 for cocaine)
                "interactions": ["MAO_inhibitors", "sympathomimetics"]
            }
        }
    
    def load_interaction_matrix(self) -> Dict[str, Dict[str, Any]]:
        """Load drug-drug interaction matrix"""
        return {
            ("buprenorphine", "ketamine"): {
                "severity": "moderate",
                "mechanism": "additive_CNS_depression",
                "effect": "Enhanced sedation and respiratory depression",
                "recommendation": "Monitor closely, consider dose reduction",
                "pk_interaction": False,
                "pd_interaction": True,
                "confidence": 0.85
            },
            ("ketamine", "psilocybin"): {
                "severity": "mild",
                "mechanism": "synergistic_neuroplasticity",
                "effect": "Enhanced neuroplasticity window, potential for improved outcomes",
                "recommendation": "May be beneficial for treatment-resistant depression",
                "pk_interaction": False,
                "pd_interaction": True,
                "confidence": 0.78
            },
            ("buprenorphine", "psilocybin"): {
                "severity": "mild",
                "mechanism": "complementary_pathways",
                "effect": "Buprenorphine may reduce anxiety during psilocybin experience",
                "recommendation": "Potential therapeutic combination for addiction treatment",
                "pk_interaction": False,
                "pd_interaction": True,
                "confidence": 0.72
            }
        }
    
    def initialize_pkpd_models(self) -> Dict[str, Any]:
        """Initialize pharmacokinetic/pharmacodynamic models"""
        return {
            "one_compartment": self.one_compartment_model,
            "two_compartment": self.two_compartment_model,
            "emax_model": self.emax_model,
            "sigmoid_emax": self.sigmoid_emax_model
        }
    
    def one_compartment_model(self, dose: float, ka: float, ke: float, vd: float, t: np.ndarray) -> np.ndarray:
        """One-compartment pharmacokinetic model"""
        if ka == ke:
            # Special case to avoid division by zero
            conc = (dose / vd) * t * ke * np.exp(-ke * t)
        else:
            conc = (dose * ka / (vd * (ka - ke))) * (np.exp(-ke * t) - np.exp(-ka * t))
        return conc
    
    def two_compartment_model(self, dose: float, ka: float, alpha: float, beta: float, 
                             a: float, b: float, t: np.ndarray) -> np.ndarray:
        """Two-compartment pharmacokinetic model"""
        conc = dose * ka * ((a / (ka - alpha)) * np.exp(-alpha * t) + 
                           (b / (ka - beta)) * np.exp(-beta * t) - 
                           ((a + b) / ka) * np.exp(-ka * t))
        return conc
    
    def emax_model(self, conc: np.ndarray, emax: float, ec50: float) -> np.ndarray:
        """Simple Emax pharmacodynamic model"""
        effect = (emax * conc) / (ec50 + conc)
        return effect
    
    def sigmoid_emax_model(self, conc: np.ndarray, emax: float, ec50: float, gamma: float) -> np.ndarray:
        """Sigmoid Emax pharmacodynamic model"""
        effect = (emax * conc**gamma) / (ec50**gamma + conc**gamma)
        return effect
    
    def simulate_pk_profile(self, drug: str, dose: float, route: str = "oral", 
                           duration: float = 24) -> Dict[str, Any]:
        """Simulate pharmacokinetic profile for a drug"""
        if drug not in self.drug_database:
            return {"error": f"Drug {drug} not found in database"}
        
        drug_data = self.drug_database[drug]
        t = np.linspace(0, duration, int(duration * 60))  # 1-minute intervals
        
        # Get PK parameters
        bioavail = drug_data["bioavailability"].get(route, 0.5)
        vd = drug_data["vd"]
        clearance = drug_data["clearance"]
        ke = clearance / vd
        
        # Estimate absorption rate constant
        ka_dict = {"oral": 1.5, "sublingual": 3.0, "nasal": 4.0, "iv": 100, "im": 2.0}
        ka = ka_dict.get(route, 1.5)
        
        # Calculate concentration
        effective_dose = dose * bioavail
        if route == "iv":
            # IV bolus
            conc = (effective_dose / vd) * np.exp(-ke * t)
        else:
            # First-order absorption
            conc = self.one_compartment_model(effective_dose, ka, ke, vd, t)
        
        # Calculate key PK parameters
        cmax = np.max(conc)
        tmax = t[np.argmax(conc)]
        auc = np.trapz(conc, t)
        
        return {
            "drug": drug,
            "dose": dose,
            "route": route,
            "time": t.tolist(),
            "concentration": conc.tolist(),
            "cmax": cmax,
            "tmax": tmax,
            "auc": auc,
            "therapeutic_range": drug_data.get("therapeutic_range", {}),
            "toxic_level": drug_data.get("toxic_level", None)
        }
    
    def analyze_drug_interaction(self, drug1: str, drug2: str, 
                               dose1: float, dose2: float) -> Dict[str, Any]:
        """Analyze drug-drug interaction between two compounds"""
        interaction_key = (drug1, drug2) if (drug1, drug2) in self.interaction_matrix else (drug2, drug1)
        
        if interaction_key not in self.interaction_matrix:
            return {
                "interaction_found": False,
                "message": f"No known interaction between {drug1} and {drug2}",
                "recommendation": "Monitor for unexpected effects"
            }
        
        interaction = self.interaction_matrix[interaction_key]
        
        # Simulate individual PK profiles
        pk1 = self.simulate_pk_profile(drug1, dose1)
        pk2 = self.simulate_pk_profile(drug2, dose2)
        
        # Calculate interaction severity score
        severity_scores = {"mild": 1, "moderate": 2, "severe": 3}
        severity_score = severity_scores.get(interaction["severity"], 1)
        
        # Estimate combined effect
        combined_analysis = {
            "interaction_found": True,
            "severity": interaction["severity"],
            "mechanism": interaction["mechanism"],
            "effect": interaction["effect"],
            "recommendation": interaction["recommendation"],
            "confidence": interaction["confidence"],
            "pk_interaction": interaction["pk_interaction"],
            "pd_interaction": interaction["pd_interaction"],
            "drug1_profile": pk1,
            "drug2_profile": pk2,
            "monitoring_parameters": self.get_monitoring_parameters(drug1, drug2),
            "dose_adjustments": self.suggest_dose_adjustments(drug1, drug2, severity_score)
        }
        
        return combined_analysis
    
    def get_monitoring_parameters(self, drug1: str, drug2: str) -> List[str]:
        """Get monitoring parameters for drug combination"""
        monitoring = []
        
        # Drug-specific monitoring
        drug_monitoring = {
            "buprenorphine": ["respiratory_rate", "sedation_level", "pupil_size"],
            "ketamine": ["blood_pressure", "heart_rate", "dissociation_score"],
            "psilocybin": ["blood_pressure", "anxiety_level", "perceptual_changes"],
            "tropacocaine": ["heart_rate", "blood_pressure", "temperature"]
        }
        
        monitoring.extend(drug_monitoring.get(drug1, []))
        monitoring.extend(drug_monitoring.get(drug2, []))
        
        # Remove duplicates
        return list(set(monitoring))
    
    def suggest_dose_adjustments(self, drug1: str, drug2: str, severity_score: int) -> Dict[str, str]:
        """Suggest dose adjustments based on interaction severity"""
        adjustments = {}
        
        if severity_score >= 2:  # Moderate or severe
            adjustments[drug1] = "Consider 25-50% dose reduction"
            adjustments[drug2] = "Consider 25-50% dose reduction"
        elif severity_score == 1:  # Mild
            adjustments[drug1] = "Monitor closely, no initial adjustment needed"
            adjustments[drug2] = "Monitor closely, no initial adjustment needed"
        
        return adjustments
    
    def neuroplasticity_window_analysis(self, drugs: List[str]) -> Dict[str, Any]:
        """Analyze neuroplasticity windows for combination therapy"""
        windows = {}
        
        for drug in drugs:
            if drug in self.drug_database:
                drug_data = self.drug_database[drug]
                if "neuroplasticity_window" in drug_data:
                    windows[drug] = drug_data["neuroplasticity_window"]
        
        if not windows:
            return {"message": "No neuroplasticity data available for these drugs"}
        
        # Calculate optimal timing
        max_window = max(windows.values())
        overlap_analysis = {
            "individual_windows": windows,
            "maximum_window": max_window,
            "optimal_timing": "Administer drugs simultaneously for maximum overlap",
            "therapeutic_implications": "Extended neuroplasticity window may enhance TMS efficacy",
            "monitoring_duration": f"{max_window + 24} hours post-administration"
        }
        
        return overlap_analysis
    
    def generate_comprehensive_report(self, drugs: List[str], doses: List[float], 
                                    routes: List[str] = None) -> Dict[str, Any]:
        """Generate comprehensive PKPD/DDI analysis report"""
        if routes is None:
            routes = ["oral"] * len(drugs)
        
        report = {
            "analysis_date": datetime.datetime.now().isoformat(),
            "drugs_analyzed": drugs,
            "doses": doses,
            "routes": routes,
            "individual_profiles": [],
            "interactions": [],
            "neuroplasticity_analysis": {},
            "clinical_recommendations": [],
            "monitoring_plan": {}
        }
        
        # Individual drug profiles
        for i, drug in enumerate(drugs):
            profile = self.simulate_pk_profile(drug, doses[i], routes[i])
            report["individual_profiles"].append(profile)
        
        # Pairwise interactions
        for i in range(len(drugs)):
            for j in range(i + 1, len(drugs)):
                interaction = self.analyze_drug_interaction(drugs[i], drugs[j], doses[i], doses[j])
                if interaction.get("interaction_found"):
                    report["interactions"].append(interaction)
        
        # Neuroplasticity analysis
        if any(drug in ["ketamine", "psilocybin", "LSD"] for drug in drugs):
            report["neuroplasticity_analysis"] = self.neuroplasticity_window_analysis(drugs)
        
        # Clinical recommendations
        report["clinical_recommendations"] = self.generate_clinical_recommendations(drugs, report["interactions"])
        
        return report
    
    def generate_clinical_recommendations(self, drugs: List[str], interactions: List[Dict]) -> List[str]:
        """Generate clinical recommendations based on analysis"""
        recommendations = []
        
        # General recommendations
        recommendations.append("Obtain baseline vital signs and mental status assessment")
        recommendations.append("Ensure appropriate monitoring equipment is available")
        
        # Drug-specific recommendations
        if "buprenorphine" in drugs:
            recommendations.append("Monitor respiratory function closely")
            recommendations.append("Have naloxone readily available")
        
        if "ketamine" in drugs:
            recommendations.append("Monitor for dissociative effects")
            recommendations.append("Ensure safe, supervised environment")
        
        if "psilocybin" in drugs:
            recommendations.append("Provide psychological support during experience")
            recommendations.append("Monitor for anxiety or panic reactions")
        
        # Interaction-based recommendations
        for interaction in interactions:
            if interaction["severity"] in ["moderate", "severe"]:
                recommendations.append(f"Special attention to {interaction['mechanism']}")
        
        return recommendations

# Initialize the enhanced PKPD system
def create_enhanced_pkpd_system():
    """Create and return an enhanced PKPD system instance"""
    return EnhancedPKPDSystem()

if __name__ == "__main__":
    # Test the system
    pkpd_system = create_enhanced_pkpd_system()
    
    # Test drug interaction analysis
    print("Analyzing buprenorphine + ketamine interaction...")
    interaction = pkpd_system.analyze_drug_interaction("buprenorphine", "ketamine", 8, 0.5)
    print(f"Interaction severity: {interaction.get('severity', 'Unknown')}")
    print(f"Mechanism: {interaction.get('mechanism', 'Unknown')}")
    
    # Test comprehensive report
    print("\nGenerating comprehensive report...")
    report = pkpd_system.generate_comprehensive_report(
        ["buprenorphine", "ketamine"], 
        [8, 0.5], 
        ["sublingual", "iv"]
    )
    print(f"Found {len(report['interactions'])} interactions")
    print(f"Generated {len(report['clinical_recommendations'])} recommendations")
