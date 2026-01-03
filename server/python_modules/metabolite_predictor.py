"""
Metabolite Prediction Module using RDKit
Implements Phase I and Phase II metabolic transformations based on CYP450 enzymes
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import List, Dict, Tuple
import json

class MetabolitePredictor:
    """
    Predicts metabolites using SMARTS-based reaction rules for common CYP450 transformations
    """
    
    def __init__(self):
        # Phase I CYP450 transformation rules (SMARTS patterns)
        # Format: (name, reactant_SMARTS, product_SMARTS, probability, CYP_enzyme)
        self.phase1_rules = [
            # Aromatic hydroxylation
            ("Aromatic_hydroxylation", "[cH:1]", "[c:1]O", 0.85, "CYP3A4"),
            # Aliphatic hydroxylation
            ("Aliphatic_hydroxylation", "[CH3:1]", "[CH2:1]O", 0.75, "CYP2D6"),
            ("Aliphatic_hydroxylation_2", "[CH2:1]", "[CH:1]O", 0.70, "CYP2D6"),
            # N-dealkylation
            ("N_dealkylation", "[N:1][CH3:2]", "[N:1]", 0.80, "CYP3A4"),
            ("N_dealkylation_2", "[N:1][CH2:2][CH3:3]", "[N:1]", 0.75, "CYP2C19"),
            # O-dealkylation
            ("O_dealkylation", "[O:1][CH3:2]", "[O:1]", 0.78, "CYP2D6"),
            # S-oxidation
            ("S_oxidation", "[S:1]", "[S:1](=O)", 0.70, "CYP3A4"),
            # N-oxidation
            ("N_oxidation", "[N:1]([C:2])([C:3])", "[N+:1]([C:2])([C:3])[O-]", 0.65, "FMO3"),
            # Epoxidation
            ("Epoxidation", "[C:1]=[C:2]", "[C:1]1[O][C:2]1", 0.60, "CYP2E1"),
            # Ketone reduction
            ("Ketone_reduction", "[C:1](=O)[C:2]", "[C:1](O)[C:2]", 0.55, "AKR"),
            # Ester hydrolysis
            ("Ester_hydrolysis", "[C:1](=O)[O][C:2]", "[C:1](=O)O", 0.85, "Esterase"),
            # Amide hydrolysis
            ("Amide_hydrolysis", "[C:1](=O)[N:2]", "[C:1](=O)O", 0.50, "Amidase"),
        ]
        
        # Phase II conjugation rules
        self.phase2_rules = [
            # Glucuronidation (on hydroxyl groups)
            ("Glucuronidation", "[OH:1]", "[O:1]C1OC(C(O)C(O)C1O)C(=O)O", 0.90, "UGT"),
            # Sulfation
            ("Sulfation", "[OH:1]", "[O:1]S(=O)(=O)O", 0.75, "SULT"),
            # Glutathione conjugation
            ("Glutathione_conjugation", "[C:1]=[C:2]", "[C:1][C:2]SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O", 0.60, "GST"),
            # Acetylation
            ("Acetylation", "[NH2:1]", "[NH:1]C(=O)C", 0.70, "NAT"),
            # Methylation
            ("Methylation", "[OH:1]", "[O:1]C", 0.65, "COMT"),
        ]
    
    def predict_metabolites(
        self,
        smiles: str,
        max_metabolites: int = 10,
        include_phase2: bool = True,
        min_probability: float = 0.5
    ) -> List[Dict]:
        """
        Predict metabolites for a given SMILES string
        
        Args:
            smiles: Input molecule SMILES
            max_metabolites: Maximum number of metabolites to return
            include_phase2: Whether to include Phase II metabolites
            min_probability: Minimum probability threshold
            
        Returns:
            List of metabolite dictionaries with structure and metadata
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        metabolites = []
        seen_smiles = set([smiles])  # Track unique metabolites
        
        # Phase I metabolism
        for rule_name, reactant_smarts, product_smarts, probability, enzyme in self.phase1_rules:
            if probability < min_probability:
                continue
                
            rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")
            if rxn is None:
                continue
            
            try:
                products = rxn.RunReactants((mol,))
                for product_set in products:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            metabolite_smiles = Chem.MolToSmiles(product)
                            
                            if metabolite_smiles not in seen_smiles:
                                seen_smiles.add(metabolite_smiles)
                                metabolites.append({
                                    "smiles": metabolite_smiles,
                                    "parent_smiles": smiles,
                                    "transformation": rule_name,
                                    "phase": "Phase I",
                                    "enzyme": enzyme,
                                    "probability": probability,
                                    "molecular_weight": Descriptors.MolWt(product),
                                    "logp": Descriptors.MolLogP(product),
                                })
                        except:
                            continue
            except:
                continue
        
        # Phase II metabolism (on Phase I metabolites)
        if include_phase2:
            phase1_metabolites = metabolites.copy()
            for phase1_met in phase1_metabolites:
                phase1_mol = Chem.MolFromSmiles(phase1_met["smiles"])
                if phase1_mol is None:
                    continue
                
                for rule_name, reactant_smarts, product_smarts, probability, enzyme in self.phase2_rules:
                    if probability < min_probability:
                        continue
                    
                    rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")
                    if rxn is None:
                        continue
                    
                    try:
                        products = rxn.RunReactants((phase1_mol,))
                        for product_set in products:
                            for product in product_set:
                                try:
                                    Chem.SanitizeMol(product)
                                    metabolite_smiles = Chem.MolToSmiles(product)
                                    
                                    if metabolite_smiles not in seen_smiles:
                                        seen_smiles.add(metabolite_smiles)
                                        metabolites.append({
                                            "smiles": metabolite_smiles,
                                            "parent_smiles": phase1_met["smiles"],
                                            "transformation": f"{phase1_met['transformation']} → {rule_name}",
                                            "phase": "Phase II",
                                            "enzyme": enzyme,
                                            "probability": phase1_met["probability"] * probability,
                                            "molecular_weight": Descriptors.MolWt(product),
                                            "logp": Descriptors.MolLogP(product),
                                        })
                                except:
                                    continue
                    except:
                        continue
        
        # Sort by probability and return top N
        metabolites.sort(key=lambda x: x["probability"], reverse=True)
        return metabolites[:max_metabolites]
    
    def predict_metabolic_stability(self, smiles: str) -> Dict:
        """
        Predict metabolic stability based on number and likelihood of metabolites
        
        Returns:
            Dictionary with stability score and analysis
        """
        metabolites = self.predict_metabolites(smiles, max_metabolites=50, min_probability=0.4)
        
        if not metabolites:
            return {
                "stability_score": 95,
                "classification": "High",
                "num_metabolites": 0,
                "avg_probability": 0,
                "analysis": "No major metabolic sites identified - likely stable"
            }
        
        num_metabolites = len(metabolites)
        avg_probability = sum(m["probability"] for m in metabolites) / num_metabolites
        
        # Calculate stability score (0-100, higher = more stable)
        stability_score = max(0, 100 - (num_metabolites * 3) - (avg_probability * 50))
        
        if stability_score >= 70:
            classification = "High"
        elif stability_score >= 40:
            classification = "Moderate"
        else:
            classification = "Low"
        
        return {
            "stability_score": round(stability_score, 1),
            "classification": classification,
            "num_metabolites": num_metabolites,
            "avg_probability": round(avg_probability, 3),
            "analysis": f"Predicted {num_metabolites} metabolites with average probability {avg_probability:.2f}"
        }


def main():
    """Test the metabolite predictor"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python metabolite_predictor.py <SMILES>")
        sys.exit(1)
    
    smiles = sys.argv[1]
    predictor = MetabolitePredictor()
    
    try:
        # Predict metabolites
        metabolites = predictor.predict_metabolites(smiles, max_metabolites=10)
        
        # Predict stability
        stability = predictor.predict_metabolic_stability(smiles)
        
        result = {
            "parent_smiles": smiles,
            "metabolic_stability": stability,
            "metabolites": metabolites
        }
        
        print(json.dumps(result, indent=2))
    except Exception as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
