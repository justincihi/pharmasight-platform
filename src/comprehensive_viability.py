#!/usr/bin/env python3
"""
Comprehensive Viability Analysis Module
Unified scoring system for analog evaluation including:
- Synthetic Accessibility (SA) Score
- Quantitative Estimate of Drug-likeness (QED)
- Natural Product-likeness (NP) Score
- Lipinski's Rule of Five
- FTO (Freedom-to-Operate) Scoring
- Patent Document Generation
"""

import os
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, AllChem, Draw
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import RDConfig
from typing import Dict, Any, List, Optional
from datetime import datetime
import json

sys.path.insert(0, os.path.join(RDConfig.RDContribDir, 'SA_Score'))
try:
    import sascorer
    SA_SCORE_AVAILABLE = True
except ImportError:
    SA_SCORE_AVAILABLE = False
    print("Warning: SA_Score not available - using approximation")

sys.path.insert(0, os.path.join(RDConfig.RDContribDir, 'NP_Score'))
try:
    import npscorer
    fscore = npscorer.readNPModel()
    NP_SCORE_AVAILABLE = True
except ImportError:
    NP_SCORE_AVAILABLE = False
    fscore = None
    print("Warning: NP_Score not available - using approximation")


class ComprehensiveViabilityAnalyzer:
    """Complete viability analysis for drug candidates"""
    
    def __init__(self):
        self.patent_reference_compounds = self._load_patent_references()
        
    def _load_patent_references(self) -> List[Dict]:
        """Load known patented compounds for FTO analysis"""
        return [
            {"name": "Psilocybin", "smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12", "patent_status": "Expired", "expiry": "2000"},
            {"name": "LSD", "smiles": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C", "patent_status": "Expired", "expiry": "1963"},
            {"name": "MDMA", "smiles": "CC(CC1=CC2=C(C=C1)OCO2)NC", "patent_status": "Expired", "expiry": "1970"},
            {"name": "DMT", "smiles": "CN(C)CCc1c[nH]c2ccccc12", "patent_status": "Natural compound", "expiry": "N/A"},
            {"name": "5-MeO-DMT", "smiles": "COc1ccc2[nH]cc(CCN(C)C)c2c1", "patent_status": "Natural compound", "expiry": "N/A"},
            {"name": "Ketamine", "smiles": "CNC1(CCCCC1=O)c2ccccc2Cl", "patent_status": "Expired", "expiry": "1966"},
            {"name": "Mescaline", "smiles": "COc1cc(CCN)cc(OC)c1OC", "patent_status": "Natural compound", "expiry": "N/A"},
            {"name": "Ibogaine", "smiles": "CCC1CC2CC3C(C(=O)OC)C4=Nc5ccc(OC)cc5C4CC3N(CC2C1)C", "patent_status": "Natural compound", "expiry": "N/A"},
        ]
    
    def calculate_sa_score(self, mol) -> Dict[str, Any]:
        """Calculate Synthetic Accessibility Score (1=easy, 10=difficult)"""
        if SA_SCORE_AVAILABLE:
            try:
                sa_score = sascorer.calculateScore(mol)
            except:
                sa_score = self._approximate_sa_score(mol)
        else:
            sa_score = self._approximate_sa_score(mol)
        
        if sa_score <= 3:
            feasibility = "Very Easy"
            interpretation = "Simple synthetic routes available"
        elif sa_score <= 5:
            feasibility = "Easy"
            interpretation = "Standard synthetic procedures"
        elif sa_score <= 7:
            feasibility = "Moderate"
            interpretation = "Requires specialized techniques"
        elif sa_score <= 9:
            feasibility = "Difficult"
            interpretation = "Complex multi-step synthesis"
        else:
            feasibility = "Very Difficult"
            interpretation = "May require novel synthetic methods"
        
        return {
            "sa_score": round(sa_score, 2),
            "feasibility": feasibility,
            "interpretation": interpretation,
            "normalized_score": round(max(0, 100 - (sa_score * 10)), 1)
        }
    
    def _approximate_sa_score(self, mol) -> float:
        """Approximate SA score when module not available"""
        num_atoms = mol.GetNumAtoms()
        num_rings = Descriptors.RingCount(mol)
        num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        num_rotatable = Descriptors.NumRotatableBonds(mol)
        
        score = 1.0
        score += num_atoms * 0.05
        score += num_rings * 0.3
        score += num_stereo * 0.5
        score += num_rotatable * 0.1
        
        return min(10, max(1, score))
    
    def calculate_qed_score(self, mol) -> Dict[str, Any]:
        """Calculate Quantitative Estimate of Drug-likeness (0-1)"""
        try:
            qed_score = QED.qed(mol)
            qed_properties = QED.properties(mol)
            
            if qed_score >= 0.7:
                interpretation = "Excellent drug-like properties"
                priority = "Very High"
            elif qed_score >= 0.5:
                interpretation = "Good drug-like properties"
                priority = "High"
            elif qed_score >= 0.3:
                interpretation = "Moderate drug-like properties"
                priority = "Medium"
            else:
                interpretation = "Poor drug-like properties"
                priority = "Low"
            
            return {
                "qed_score": round(qed_score, 3),
                "interpretation": interpretation,
                "priority": priority,
                "properties": {
                    "MW": round(qed_properties.MW, 2),
                    "ALOGP": round(qed_properties.ALOGP, 2),
                    "HBA": qed_properties.HBA,
                    "HBD": qed_properties.HBD,
                    "PSA": round(qed_properties.PSA, 2),
                    "ROTB": qed_properties.ROTB,
                    "AROM": qed_properties.AROM,
                    "ALERTS": qed_properties.ALERTS
                },
                "normalized_score": round(qed_score * 100, 1)
            }
        except Exception as e:
            return {
                "qed_score": None,
                "interpretation": f"Error: {str(e)}",
                "priority": "Unknown",
                "properties": {},
                "normalized_score": 0
            }
    
    def calculate_np_score(self, mol) -> Dict[str, Any]:
        """Calculate Natural Product-likeness Score (-5 to +5)"""
        if NP_SCORE_AVAILABLE and fscore is not None:
            try:
                np_score = npscorer.scoreMol(mol, fscore)
            except:
                np_score = self._approximate_np_score(mol)
        else:
            np_score = self._approximate_np_score(mol)
        
        if np_score >= 1.5:
            interpretation = "Highly natural product-like"
            natural_like = True
        elif np_score >= 0.5:
            interpretation = "Moderately natural product-like"
            natural_like = True
        elif np_score >= -0.5:
            interpretation = "Slightly natural product-like"
            natural_like = False
        else:
            interpretation = "Synthetic/non-natural product-like"
            natural_like = False
        
        normalized = min(100, max(0, (np_score + 5) * 10))
        
        return {
            "np_score": round(np_score, 2),
            "interpretation": interpretation,
            "natural_product_like": natural_like,
            "normalized_score": round(normalized, 1)
        }
    
    def _approximate_np_score(self, mol) -> float:
        """Approximate NP score when module not available"""
        num_rings = Descriptors.RingCount(mol)
        num_aromatic = Descriptors.NumAromaticRings(mol)
        num_hetero = Descriptors.NumHeteroatoms(mol)
        fsp3 = Descriptors.FractionCSP3(mol)
        
        score = 0.0
        score += (num_rings - num_aromatic) * 0.3
        score += fsp3 * 2.0
        score += (num_hetero / max(1, mol.GetNumAtoms())) * 1.5
        score -= num_aromatic * 0.2
        
        return min(5, max(-5, score - 1))
    
    def calculate_lipinski(self, mol) -> Dict[str, Any]:
        """Calculate Lipinski's Rule of Five compliance"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        rules = {}
        
        rules["molecular_weight"] = {
            "value": round(mw, 2),
            "limit": "≤ 500 Da",
            "passes": mw <= 500
        }
        if mw > 500: violations += 1
        
        rules["logP"] = {
            "value": round(logp, 2),
            "limit": "≤ 5",
            "passes": logp <= 5
        }
        if logp > 5: violations += 1
        
        rules["h_bond_donors"] = {
            "value": hbd,
            "limit": "≤ 5",
            "passes": hbd <= 5
        }
        if hbd > 5: violations += 1
        
        rules["h_bond_acceptors"] = {
            "value": hba,
            "limit": "≤ 10",
            "passes": hba <= 10
        }
        if hba > 10: violations += 1
        
        passes = violations <= 1
        
        if violations == 0:
            interpretation = "Fully compliant - excellent oral bioavailability expected"
        elif violations == 1:
            interpretation = "1 violation - good oral bioavailability likely"
        elif violations == 2:
            interpretation = "2 violations - oral bioavailability may be limited"
        else:
            interpretation = f"{violations} violations - poor oral bioavailability expected"
        
        return {
            "violations": violations,
            "passes_ro5": passes,
            "rules": rules,
            "interpretation": interpretation,
            "normalized_score": max(0, 100 - (violations * 25))
        }
    
    def calculate_fto_score(self, smiles: str) -> Dict[str, Any]:
        """Calculate Freedom-to-Operate score based on similarity to patented compounds"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        fp_query = FingerprintMols.FingerprintMol(mol)
        
        similarities = []
        closest_match = None
        max_similarity = 0.0
        
        for ref in self.patent_reference_compounds:
            ref_mol = Chem.MolFromSmiles(ref["smiles"])
            if ref_mol:
                ref_fp = FingerprintMols.FingerprintMol(ref_mol)
                similarity = DataStructs.TanimotoSimilarity(fp_query, ref_fp)
                similarities.append({
                    "compound": ref["name"],
                    "similarity": round(similarity, 3),
                    "patent_status": ref["patent_status"],
                    "expiry": ref["expiry"]
                })
                if similarity > max_similarity:
                    max_similarity = similarity
                    closest_match = ref
        
        similarities.sort(key=lambda x: x["similarity"], reverse=True)
        
        if max_similarity < 0.7:
            fto_score = 95
            risk_level = "Very Low"
            recommendation = "Novel structure - excellent patent opportunity"
        elif max_similarity < 0.8:
            fto_score = 75
            risk_level = "Low"
            recommendation = "Distinct from known compounds - good patent potential"
        elif max_similarity < 0.9:
            fto_score = 50
            risk_level = "Medium"
            recommendation = "Some similarity to known compounds - requires careful claims drafting"
        elif max_similarity < 0.95:
            fto_score = 25
            risk_level = "High"
            recommendation = "High similarity - may require licensing or alternative approach"
        else:
            fto_score = 10
            risk_level = "Very High"
            recommendation = "Very similar to existing compound - limited patent opportunity"
        
        if closest_match and closest_match["patent_status"] in ["Expired", "Natural compound"]:
            fto_score = min(100, fto_score + 20)
            recommendation += f" (Note: closest match '{closest_match['name']}' has {closest_match['patent_status'].lower()} status)"
        
        return {
            "fto_score": round(fto_score, 1),
            "risk_level": risk_level,
            "recommendation": recommendation,
            "closest_match": closest_match["name"] if closest_match else None,
            "max_similarity": round(max_similarity, 3),
            "similar_compounds": similarities[:5],
            "normalized_score": fto_score
        }
    
    def calculate_overall_viability(self, smiles: str, include_admet: bool = True) -> Dict[str, Any]:
        """Calculate comprehensive viability score (0-100)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES string"}
        
        sa_data = self.calculate_sa_score(mol)
        qed_data = self.calculate_qed_score(mol)
        np_data = self.calculate_np_score(mol)
        lipinski_data = self.calculate_lipinski(mol)
        fto_data = self.calculate_fto_score(smiles)
        
        scores = {}
        weights = {}
        
        scores["synthetic_accessibility"] = sa_data["normalized_score"]
        weights["synthetic_accessibility"] = 0.20
        
        scores["drug_likeness"] = qed_data["normalized_score"]
        weights["drug_likeness"] = 0.15
        
        scores["natural_product_likeness"] = np_data["normalized_score"]
        weights["natural_product_likeness"] = 0.10
        
        scores["lipinski_compliance"] = lipinski_data["normalized_score"]
        weights["lipinski_compliance"] = 0.15
        
        scores["patent_freedom"] = fto_data["normalized_score"]
        weights["patent_freedom"] = 0.20
        
        if include_admet:
            try:
                from admet_ml_predictor import ADMETMLPredictor
                admet_predictor = ADMETMLPredictor()
                admet_results = admet_predictor.predict_all(smiles)
                
                safety_score = 100
                if admet_results.get('hERG_inhibition', {}).get('prediction') == 'Inhibitor':
                    safety_score -= 30
                if admet_results.get('hepatotoxicity', {}).get('prediction') == 'Toxic':
                    safety_score -= 25
                if admet_results.get('ames_mutagenicity', {}).get('prediction') == 'Mutagenic':
                    safety_score -= 20
                
                scores["safety_profile"] = max(0, safety_score)
                weights["safety_profile"] = 0.20
            except:
                scores["safety_profile"] = 75
                weights["safety_profile"] = 0.20
        else:
            scores["safety_profile"] = 75
            weights["safety_profile"] = 0.20
        
        total_weight = sum(weights.values())
        viability_score = sum(scores[k] * weights[k] for k in scores) / total_weight
        
        if viability_score >= 80:
            priority = "Very High"
            recommendation = "Excellent candidate - proceed to synthesis and patent filing"
        elif viability_score >= 65:
            priority = "High"
            recommendation = "Good candidate - further analysis recommended before filing"
        elif viability_score >= 50:
            priority = "Medium"
            recommendation = "Moderate candidate - consider optimization before proceeding"
        else:
            priority = "Low"
            recommendation = "Poor candidate - significant improvements needed"
        
        return {
            "smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "overall_score": round(viability_score, 1),
            "priority": priority,
            "recommendation": recommendation,
            "component_analysis": {
                "synthetic_accessibility": sa_data,
                "drug_likeness": qed_data,
                "natural_product_likeness": np_data,
                "lipinski_compliance": lipinski_data,
                "patent_freedom": fto_data
            },
            "component_scores": {k: round(v, 1) for k, v in scores.items()},
            "weights": weights,
            "timestamp": datetime.now().isoformat()
        }
    
    def batch_viability_analysis(self, compounds: List[Dict]) -> Dict[str, Any]:
        """Analyze multiple compounds and rank by viability"""
        results = []
        
        for compound in compounds:
            smiles = compound.get('smiles', '')
            name = compound.get('name', 'Unknown')
            
            analysis = self.calculate_overall_viability(smiles)
            if 'error' not in analysis:
                analysis['compound_name'] = name
                results.append(analysis)
        
        results.sort(key=lambda x: x['overall_score'], reverse=True)
        
        for i, result in enumerate(results, 1):
            result['rank'] = i
        
        return {
            "total_analyzed": len(results),
            "high_priority_count": sum(1 for r in results if r['priority'] in ['Very High', 'High']),
            "results": results,
            "top_candidates": results[:10],
            "timestamp": datetime.now().isoformat()
        }


class PatentDocumentGenerator:
    """Generate provisional patent document drafts"""
    
    def __init__(self):
        self.viability_analyzer = ComprehensiveViabilityAnalyzer()
    
    def generate_provisional_patent(self, 
                                     smiles: str, 
                                     compound_name: str,
                                     parent_compound: str = None,
                                     therapeutic_area: str = None,
                                     inventor_name: str = "PharmaSight AI Discovery Engine",
                                     modification_type: str = None) -> Dict[str, Any]:
        """Generate a provisional patent document draft"""
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES string"}
        
        canonical_smiles = Chem.MolToSmiles(mol)
        viability = self.viability_analyzer.calculate_overall_viability(smiles)
        
        mw = round(Descriptors.MolWt(mol), 2)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        logp = round(Descriptors.MolLogP(mol), 2)
        tpsa = round(Descriptors.TPSA(mol), 2)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        if not therapeutic_area:
            therapeutic_area = self._predict_therapeutic_area(mol, parent_compound)
        
        title = f"{modification_type or 'Novel'} Derivative of {parent_compound or 'Tryptamine'} for {therapeutic_area}"
        
        abstract = f"""The invention relates to {compound_name}, a novel {"prodrug" if "ester" in (modification_type or "").lower() else "analog"} """
        abstract += f"""designed to improve pharmacokinetics, stability, and therapeutic utility in {therapeutic_area.lower()} applications. """
        abstract += f"""This compound has molecular formula {formula} and molecular weight {mw} Da. """
        abstract += f"""The invention provides synthesis methods, formulation guidance, pharmacodynamic rationale, and therapeutic use claims."""
        
        background = f"""Tryptamine derivatives such as {parent_compound or "the parent compound"} have shown therapeutic promise """
        background += f"""in treating various {therapeutic_area.lower().replace(" therapy", "")} disorders. """
        background += f"""However, limitations such as rapid metabolism, poor oral bioavailability, or regulatory constraints """
        background += f"""can be addressed using {"prodrug" if "ester" in (modification_type or "").lower() else "structural modification"} strategies. """
        background += f"""The {modification_type or "novel modification"} may offer improved onset/duration profiles and stability."""
        
        summary = f"""This patent discloses a novel {modification_type or "derivative"} of {parent_compound or "tryptamine"} """
        summary += f"""with potential clinical use in {therapeutic_area.lower()}. """
        summary += f"""It includes synthesis methodology, predicted pharmacokinetic improvements, and therapeutic uses. """
        summary += f"""Viability Score: {viability['overall_score']}/100 ({viability['priority']} priority)."""
        
        claims = [
            f"1. A compound comprising {compound_name} with SMILES notation: {canonical_smiles}.",
            f"2. A pharmaceutical composition comprising the compound of claim 1 and a pharmaceutically acceptable carrier.",
            f"3. A method of treating {therapeutic_area.lower().replace(' therapy', '')} disorders using the composition of claim 2.",
            f"4. The compound of claim 1, wherein the compound is formulated for oral administration.",
            f"5. The compound of claim 1, wherein the compound is formulated for parenteral administration.",
            f"6. A method of synthesis comprising the steps of preparing {compound_name} from {parent_compound or 'a suitable precursor'}."
        ]
        
        patent_doc = {
            "document_type": "Provisional Patent Application",
            "filing_date": datetime.now().strftime("%B %d, %Y"),
            "title": title,
            "inventor": inventor_name,
            "compound_info": {
                "name": compound_name,
                "smiles": canonical_smiles,
                "molecular_formula": formula,
                "molecular_weight": mw,
                "logP": logp,
                "tpsa": tpsa,
                "h_bond_donors": hbd,
                "h_bond_acceptors": hba
            },
            "sections": {
                "abstract": abstract,
                "background": background,
                "summary": summary,
                "claims": claims
            },
            "viability_analysis": viability,
            "therapeutic_area": therapeutic_area,
            "parent_compound": parent_compound,
            "modification_type": modification_type,
            "priority_recommendation": self._get_filing_recommendation(viability),
            "generated_by": "PharmaSight™ Patent Document Generator",
            "disclaimer": "This is an AI-generated draft for review. Consult a patent attorney before filing."
        }
        
        return patent_doc
    
    def _predict_therapeutic_area(self, mol, parent_compound: str = None) -> str:
        """Predict likely therapeutic area based on structure"""
        smiles = Chem.MolToSmiles(mol)
        
        if "[nH]c" in smiles.lower() or "c1c[nH]c2" in smiles.lower():
            return "Psychedelic Therapy"
        elif "OCO" in smiles:
            return "PTSD Therapy"
        elif "N(C)C" in smiles:
            return "Neuropsychiatric Therapy"
        elif parent_compound:
            therapeutic_map = {
                "psilocybin": "Psychedelic Therapy",
                "dmt": "Psychedelic Therapy",
                "5-meo-dmt": "Psychedelic Therapy",
                "mdma": "PTSD Therapy",
                "ketamine": "Depression Therapy",
                "lsd": "Psychedelic Therapy"
            }
            return therapeutic_map.get(parent_compound.lower(), "Neuropsychiatric Therapy")
        else:
            return "Neuropsychiatric Therapy"
    
    def _get_filing_recommendation(self, viability: Dict) -> Dict[str, Any]:
        """Get patent filing recommendation based on viability"""
        score = viability.get('overall_score', 0)
        fto = viability.get('component_analysis', {}).get('patent_freedom', {})
        
        if score >= 80 and fto.get('fto_score', 0) >= 70:
            return {
                "urgency": "Immediate",
                "recommendation": "File provisional patent application immediately",
                "reason": "High viability and excellent patent freedom",
                "estimated_value": "$50,000 - $200,000"
            }
        elif score >= 65 and fto.get('fto_score', 0) >= 50:
            return {
                "urgency": "High",
                "recommendation": "Prepare patent application within 2 weeks",
                "reason": "Good viability with moderate patent freedom",
                "estimated_value": "$25,000 - $100,000"
            }
        elif score >= 50:
            return {
                "urgency": "Medium",
                "recommendation": "Consider filing after optimization",
                "reason": "Moderate viability - may benefit from improvements",
                "estimated_value": "$10,000 - $50,000"
            }
        else:
            return {
                "urgency": "Low",
                "recommendation": "Focus on optimization before filing",
                "reason": "Low viability - significant improvements needed",
                "estimated_value": "< $10,000"
            }
    
    def generate_patent_text(self, patent_doc: Dict) -> str:
        """Generate formatted patent document text"""
        sections = patent_doc['sections']
        compound = patent_doc['compound_info']
        
        text = f"""
================================================================================
                    PROVISIONAL PATENT APPLICATION
================================================================================

Title: {patent_doc['title']}

Filing Date: {patent_doc['filing_date']}
Inventor(s): {patent_doc['inventor']}

================================================================================
                              ABSTRACT
================================================================================

{sections['abstract']}

================================================================================
                          CHEMICAL STRUCTURE
================================================================================

Compound Name: {compound['name']}
SMILES: {compound['smiles']}
Molecular Formula: {compound['molecular_formula']}
Molecular Weight: {compound['molecular_weight']} Da

Figure 1. Structure of {compound['name']}

================================================================================
                           BACKGROUND
================================================================================

{sections['background']}

================================================================================
                      SUMMARY OF THE INVENTION
================================================================================

{sections['summary']}

================================================================================
                             CLAIMS
================================================================================

"""
        for claim in sections['claims']:
            text += f"{claim}\n\n"
        
        text += f"""
================================================================================
                      VIABILITY ANALYSIS
================================================================================

Overall Viability Score: {patent_doc['viability_analysis']['overall_score']}/100
Priority Level: {patent_doc['viability_analysis']['priority']}
Recommendation: {patent_doc['viability_analysis']['recommendation']}

Filing Recommendation: {patent_doc['priority_recommendation']['recommendation']}
Urgency: {patent_doc['priority_recommendation']['urgency']}
Estimated IP Value: {patent_doc['priority_recommendation']['estimated_value']}

================================================================================
                           DISCLAIMER
================================================================================

{patent_doc['disclaimer']}

Generated by: {patent_doc['generated_by']}
"""
        return text


def analyze_compound_viability(smiles: str) -> Dict[str, Any]:
    """Convenience function for quick viability analysis"""
    analyzer = ComprehensiveViabilityAnalyzer()
    return analyzer.calculate_overall_viability(smiles)


def generate_patent_draft(smiles: str, compound_name: str, 
                          parent_compound: str = None,
                          therapeutic_area: str = None,
                          modification_type: str = None) -> Dict[str, Any]:
    """Convenience function for generating patent drafts"""
    generator = PatentDocumentGenerator()
    return generator.generate_provisional_patent(
        smiles=smiles,
        compound_name=compound_name,
        parent_compound=parent_compound,
        therapeutic_area=therapeutic_area,
        modification_type=modification_type
    )
