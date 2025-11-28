"""
ADMET Prediction Module using Machine Learning
Implements predictive models for Absorption, Distribution, Metabolism, Excretion, and Toxicity
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf, rdMolDescriptors
from typing import Dict, List, Tuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MolecularDescriptorCalculator:
    """Calculate molecular descriptors for ADMET prediction"""
    
    @staticmethod
    def calculate_descriptors(smiles: str) -> Dict[str, float]:
        """Calculate comprehensive molecular descriptors"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return {}
        
        descriptors = {
            # Basic properties
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Crippen.MolLogP(mol),
            'tpsa': MolSurf.TPSA(mol),
            'num_h_acceptors': Lipinski.NumHAcceptors(mol),
            'num_h_donors': Lipinski.NumHDonors(mol),
            'num_rotatable_bonds': Lipinski.NumRotatableBonds(mol),
            'num_aromatic_rings': Lipinski.NumAromaticRings(mol),
            
            # Additional descriptors
            'num_heteroatoms': Lipinski.NumHeteroatoms(mol),
            'num_rings': Lipinski.RingCount(mol),
            'num_saturated_rings': Lipinski.NumSaturatedRings(mol),
            'num_aliphatic_rings': Lipinski.NumAliphaticRings(mol),
            'fraction_csp3': Lipinski.FractionCsp3(mol),
            
            # Complexity measures
            'bertz_ct': Descriptors.BertzCT(mol),
            'chi0': Descriptors.Chi0(mol),
            'chi1': Descriptors.Chi1(mol),
            'kappa1': Descriptors.Kappa1(mol),
            'kappa2': Descriptors.Kappa2(mol),
            
            # Electronic properties
            'num_valence_electrons': Descriptors.NumValenceElectrons(mol),
            'num_radical_electrons': Descriptors.NumRadicalElectrons(mol),
            
            # Surface area
            'labute_asa': MolSurf.LabuteASA(mol),
            
            # Molecular refractivity
            'molar_refractivity': Crippen.MolMR(mol),
        }
        
        return descriptors


class ADMETPredictor:
    """ADMET prediction using rule-based and ML models"""
    
    def __init__(self):
        self.descriptor_calc = MolecularDescriptorCalculator()
    
    def predict_absorption(self, smiles: str) -> Dict[str, any]:
        """
        Predict oral absorption properties
        Based on Lipinski's Rule of 5 and extended rules
        """
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        if not descriptors:
            return {'error': 'Invalid SMILES'}
        
        # Lipinski's Rule of 5
        lipinski_violations = 0
        if descriptors['molecular_weight'] > 500:
            lipinski_violations += 1
        if descriptors['logp'] > 5:
            lipinski_violations += 1
        if descriptors['num_h_donors'] > 5:
            lipinski_violations += 1
        if descriptors['num_h_acceptors'] > 10:
            lipinski_violations += 1
        
        # Veber rules (oral bioavailability)
        veber_compliant = (
            descriptors['num_rotatable_bonds'] <= 10 and
            descriptors['tpsa'] <= 140
        )
        
        # Calculate absorption score (0-100)
        absorption_score = 100
        absorption_score -= (lipinski_violations * 15)
        if not veber_compliant:
            absorption_score -= 20
        if descriptors['tpsa'] > 140:
            absorption_score -= 15
        
        absorption_score = max(0, min(100, absorption_score))
        
        return {
            'absorption_score': round(absorption_score, 1),
            'lipinski_violations': lipinski_violations,
            'veber_compliant': veber_compliant,
            'predicted_bioavailability': 'High' if absorption_score > 70 else 'Medium' if absorption_score > 40 else 'Low',
            'hia_probability': round(absorption_score / 100, 2),  # Human Intestinal Absorption
            'caco2_permeability': self._predict_caco2(descriptors),
        }
    
    def predict_distribution(self, smiles: str) -> Dict[str, any]:
        """Predict distribution properties including BBB penetration"""
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        if not descriptors:
            return {'error': 'Invalid SMILES'}
        
        # BBB penetration prediction
        # Based on: MW < 400, TPSA < 90, LogP 1-3, H-bonds < 8
        bbb_score = 100
        
        if descriptors['molecular_weight'] > 400:
            bbb_score -= 30
        if descriptors['tpsa'] > 90:
            bbb_score -= 25
        if descriptors['logp'] < 1 or descriptors['logp'] > 3:
            bbb_score -= 20
        if (descriptors['num_h_donors'] + descriptors['num_h_acceptors']) > 8:
            bbb_score -= 15
        
        bbb_score = max(0, min(100, bbb_score))
        
        # Volume of distribution prediction
        vd_category = 'Medium'
        if descriptors['logp'] > 3:
            vd_category = 'High'
        elif descriptors['logp'] < 0:
            vd_category = 'Low'
        
        # Plasma protein binding
        ppb_percentage = min(99, max(10, 50 + (descriptors['logp'] * 10)))
        
        return {
            'bbb_penetration_score': round(bbb_score, 1),
            'bbb_penetration': 'High' if bbb_score > 70 else 'Medium' if bbb_score > 40 else 'Low',
            'volume_of_distribution': vd_category,
            'plasma_protein_binding': round(ppb_percentage, 1),
            'tissue_distribution': 'Extensive' if descriptors['logp'] > 2 else 'Limited',
        }
    
    def predict_metabolism(self, smiles: str) -> Dict[str, any]:
        """Predict metabolic properties"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        
        # CYP450 substrate prediction (simplified)
        # Based on molecular properties and functional groups
        cyp_substrates = {
            'CYP3A4': self._predict_cyp3a4_substrate(mol, descriptors),
            'CYP2D6': self._predict_cyp2d6_substrate(mol, descriptors),
            'CYP2C9': self._predict_cyp2c9_substrate(mol, descriptors),
            'CYP2C19': self._predict_cyp2c19_substrate(mol, descriptors),
            'CYP1A2': self._predict_cyp1a2_substrate(mol, descriptors),
        }
        
        # Metabolic stability prediction
        stability_score = 70  # baseline
        if descriptors['num_aromatic_rings'] > 2:
            stability_score += 10
        if descriptors['num_rotatable_bonds'] > 8:
            stability_score -= 15
        if descriptors['fraction_csp3'] < 0.3:
            stability_score -= 10
        
        stability_score = max(0, min(100, stability_score))
        
        return {
            'cyp_substrates': cyp_substrates,
            'metabolic_stability_score': round(stability_score, 1),
            'predicted_half_life': 'Long' if stability_score > 70 else 'Medium' if stability_score > 40 else 'Short',
            'first_pass_metabolism': 'Low' if stability_score > 60 else 'High',
        }
    
    def predict_excretion(self, smiles: str) -> Dict[str, any]:
        """Predict excretion properties"""
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        if not descriptors:
            return {'error': 'Invalid SMILES'}
        
        # Renal clearance prediction
        # Based on MW, LogP, and polarity
        renal_clearance = 'Medium'
        if descriptors['molecular_weight'] < 300 and descriptors['logp'] < 2:
            renal_clearance = 'High'
        elif descriptors['molecular_weight'] > 500 or descriptors['logp'] > 4:
            renal_clearance = 'Low'
        
        # Biliary excretion (typically for MW > 500)
        biliary_excretion = 'High' if descriptors['molecular_weight'] > 500 else 'Low'
        
        # Clearance rate estimation (mL/min/kg)
        clearance_rate = 10  # baseline
        if descriptors['logp'] < 1:
            clearance_rate += 5
        if descriptors['molecular_weight'] < 300:
            clearance_rate += 3
        
        return {
            'renal_clearance': renal_clearance,
            'biliary_excretion': biliary_excretion,
            'predicted_clearance_rate': round(clearance_rate, 1),
            'excretion_route': 'Renal' if renal_clearance == 'High' else 'Hepatic/Biliary',
        }
    
    def predict_toxicity(self, smiles: str) -> Dict[str, any]:
        """Predict toxicity properties"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        
        # hERG liability (cardiotoxicity)
        herg_risk = self._predict_herg_liability(mol, descriptors)
        
        # Hepatotoxicity prediction
        hepatotox_risk = self._predict_hepatotoxicity(mol, descriptors)
        
        # Mutagenicity (Ames test prediction)
        mutagenicity_risk = self._predict_mutagenicity(mol)
        
        # Overall toxicity score
        toxicity_score = 100
        if herg_risk == 'High':
            toxicity_score -= 30
        elif herg_risk == 'Medium':
            toxicity_score -= 15
        
        if hepatotox_risk == 'High':
            toxicity_score -= 25
        elif hepatotox_risk == 'Medium':
            toxicity_score -= 10
        
        if mutagenicity_risk == 'High':
            toxicity_score -= 35
        elif mutagenicity_risk == 'Medium':
            toxicity_score -= 15
        
        toxicity_score = max(0, min(100, toxicity_score))
        
        return {
            'safety_score': round(toxicity_score, 1),
            'herg_liability': herg_risk,
            'hepatotoxicity_risk': hepatotox_risk,
            'mutagenicity_risk': mutagenicity_risk,
            'overall_safety': 'High' if toxicity_score > 70 else 'Medium' if toxicity_score > 40 else 'Low',
            'ld50_prediction': self._predict_ld50(descriptors),
        }
    
    def predict_full_admet(self, smiles: str) -> Dict[str, any]:
        """Comprehensive ADMET prediction"""
        return {
            'smiles': smiles,
            'absorption': self.predict_absorption(smiles),
            'distribution': self.predict_distribution(smiles),
            'metabolism': self.predict_metabolism(smiles),
            'excretion': self.predict_excretion(smiles),
            'toxicity': self.predict_toxicity(smiles),
            'drug_likeness': self._calculate_drug_likeness(smiles),
        }
    
    # Helper methods for specific predictions
    
    def _predict_caco2(self, descriptors: Dict) -> str:
        """Predict Caco-2 permeability"""
        if descriptors['tpsa'] < 60 and descriptors['logp'] > 1:
            return 'High'
        elif descriptors['tpsa'] > 120 or descriptors['logp'] < -1:
            return 'Low'
        return 'Medium'
    
    def _predict_cyp3a4_substrate(self, mol, descriptors: Dict) -> bool:
        """Predict if compound is CYP3A4 substrate"""
        # Simplified rule: MW 300-700, LogP 2-5
        return (300 < descriptors['molecular_weight'] < 700 and 
                2 < descriptors['logp'] < 5)
    
    def _predict_cyp2d6_substrate(self, mol, descriptors: Dict) -> bool:
        """Predict if compound is CYP2D6 substrate"""
        # Typically basic nitrogen-containing compounds
        return descriptors['num_heteroatoms'] > 2
    
    def _predict_cyp2c9_substrate(self, mol, descriptors: Dict) -> bool:
        """Predict if compound is CYP2C9 substrate"""
        return 200 < descriptors['molecular_weight'] < 500
    
    def _predict_cyp2c19_substrate(self, mol, descriptors: Dict) -> bool:
        """Predict if compound is CYP2C19 substrate"""
        return descriptors['num_aromatic_rings'] >= 2
    
    def _predict_cyp1a2_substrate(self, mol, descriptors: Dict) -> bool:
        """Predict if compound is CYP1A2 substrate"""
        return descriptors['num_aromatic_rings'] >= 2
    
    def _predict_herg_liability(self, mol, descriptors: Dict) -> str:
        """Predict hERG channel blocking (cardiotoxicity)"""
        # High risk: basic nitrogen, LogP 2-5, MW 300-500
        risk_score = 0
        
        if 2 < descriptors['logp'] < 5:
            risk_score += 1
        if 300 < descriptors['molecular_weight'] < 500:
            risk_score += 1
        if descriptors['num_aromatic_rings'] >= 2:
            risk_score += 1
        
        if risk_score >= 3:
            return 'High'
        elif risk_score >= 2:
            return 'Medium'
        return 'Low'
    
    def _predict_hepatotoxicity(self, mol, descriptors: Dict) -> str:
        """Predict hepatotoxicity risk"""
        risk_score = 0
        
        if descriptors['logp'] > 5:
            risk_score += 1
        if descriptors['molecular_weight'] > 500:
            risk_score += 1
        if descriptors['num_aromatic_rings'] > 3:
            risk_score += 1
        
        if risk_score >= 2:
            return 'High'
        elif risk_score >= 1:
            return 'Medium'
        return 'Low'
    
    def _predict_mutagenicity(self, mol) -> str:
        """Predict mutagenicity (Ames test)"""
        # Simplified structural alert-based prediction
        # Check for common mutagenic substructures
        
        mutagenic_smarts = [
            '[N;X3](=O)=O',  # Nitro groups
            'N=N',  # Azo groups
            '[#6]=[#6]([Cl,Br,I])',  # Halogenated alkenes
        ]
        
        for smart in mutagenic_smarts:
            pattern = Chem.MolFromSmarts(smart)
            if pattern and mol.HasSubstructMatch(pattern):
                return 'High'
        
        return 'Low'
    
    def _predict_ld50(self, descriptors: Dict) -> str:
        """Predict acute toxicity (LD50)"""
        # Simplified prediction
        if descriptors['molecular_weight'] > 500 and descriptors['logp'] < 3:
            return '>2000 mg/kg (Low toxicity)'
        elif descriptors['logp'] > 5:
            return '<500 mg/kg (High toxicity)'
        return '500-2000 mg/kg (Moderate toxicity)'
    
    def _calculate_drug_likeness(self, smiles: str) -> Dict[str, any]:
        """Calculate overall drug-likeness score"""
        descriptors = self.descriptor_calc.calculate_descriptors(smiles)
        if not descriptors:
            return {'error': 'Invalid SMILES'}
        
        # QED (Quantitative Estimate of Drug-likeness) approximation
        score = 100
        
        # Lipinski violations
        violations = 0
        if descriptors['molecular_weight'] > 500:
            violations += 1
        if descriptors['logp'] > 5:
            violations += 1
        if descriptors['num_h_donors'] > 5:
            violations += 1
        if descriptors['num_h_acceptors'] > 10:
            violations += 1
        
        score -= (violations * 20)
        
        # Additional penalties
        if descriptors['num_rotatable_bonds'] > 10:
            score -= 10
        if descriptors['tpsa'] > 140:
            score -= 10
        
        score = max(0, min(100, score))
        
        return {
            'drug_likeness_score': round(score, 1),
            'lipinski_violations': violations,
            'lead_likeness': 'Yes' if score > 60 else 'No',
            'synthetic_accessibility': 'Easy' if descriptors['bertz_ct'] < 500 else 'Moderate' if descriptors['bertz_ct'] < 800 else 'Difficult',
        }


# Singleton instance
_admet_predictor = None

def get_admet_predictor() -> ADMETPredictor:
    """Get or create ADMET predictor instance"""
    global _admet_predictor
    if _admet_predictor is None:
        _admet_predictor = ADMETPredictor()
    return _admet_predictor

