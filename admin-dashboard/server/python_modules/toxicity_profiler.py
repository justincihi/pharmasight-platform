"""
Comprehensive Toxicity Profiling Module
Predicts multiple toxicity endpoints for drug candidates
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
import numpy as np
from typing import Dict, List, Optional

class ToxicityProfiler:
    """
    Comprehensive toxicity prediction for drug candidates
    """
    
    def __init__(self):
        self.toxicity_rules = self._load_toxicity_rules()
    
    def _load_toxicity_rules(self) -> Dict:
        """Load structural alerts and toxicity rules"""
        return {
            'hERG_blockers': [
                # Tertiary amines
                '[NX3;H0;!$(NC=O)]',
                # Basic nitrogen with aromatic rings
                '[nX3;H1,H0;!$(nc=O)]',
            ],
            'hepatotoxic_patterns': [
                # Nitro aromatics
                '[c][N+](=O)[O-]',
                # Halogenated aromatics
                '[c][F,Cl,Br,I]',
                # Aromatic amines
                '[c][NH2]',
            ],
            'mutagenic_patterns': [
                # Aromatic nitro compounds
                '[c][N+](=O)[O-]',
                # Aromatic amines
                '[c][NH2]',
                # Epoxides
                'C1OC1',
                # Aziridines
                'C1NC1',
            ],
            'carcinogenic_patterns': [
                # Aromatic amines
                '[c][NH2]',
                # N-nitroso compounds
                '[N;X3][N+](=O)[O-]',
                # Polycyclic aromatic hydrocarbons (simplified)
                'c1ccc2c(c1)ccc3c2cccc3',
            ]
        }
    
    def predict_hERG_toxicity(self, smiles: str) -> Dict:
        """
        Predict hERG channel blocking potential
        hERG blocking can cause cardiac arrhythmias
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Calculate relevant descriptors
        logP = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        num_basic_groups = rdMolDescriptors.CalcNumLipinskiHBA(mol)
        
        # Check for structural alerts
        alerts = []
        for pattern in self.toxicity_rules['hERG_blockers']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                alerts.append(pattern)
        
        # Risk scoring (0-100)
        risk_score = 0
        
        # LogP contribution (hydrophobic drugs more likely to block hERG)
        if logP > 3:
            risk_score += min(30, (logP - 3) * 10)
        
        # MW contribution
        if 300 < mw < 500:
            risk_score += 20
        
        # Basic groups (cationic at physiological pH)
        if num_basic_groups > 0:
            risk_score += min(25, num_basic_groups * 12)
        
        # Aromatic rings
        if num_aromatic_rings >= 2:
            risk_score += min(15, num_aromatic_rings * 5)
        
        # Structural alerts
        risk_score += len(alerts) * 10
        
        risk_score = min(100, risk_score)
        
        # Classify risk level
        if risk_score < 30:
            risk_level = 'Low'
            recommendation = 'Acceptable hERG profile'
        elif risk_score < 60:
            risk_level = 'Medium'
            recommendation = 'Consider hERG testing'
        else:
            risk_level = 'High'
            recommendation = 'High risk - structural modification recommended'
        
        return {
            'risk_score': round(risk_score, 1),
            'risk_level': risk_level,
            'recommendation': recommendation,
            'logP': round(logP, 2),
            'molecular_weight': round(mw, 2),
            'tpsa': round(tpsa, 2),
            'aromatic_rings': num_aromatic_rings,
            'basic_groups': num_basic_groups,
            'structural_alerts': len(alerts)
        }
    
    def predict_hepatotoxicity(self, smiles: str) -> Dict:
        """
        Predict liver toxicity potential
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        logP = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        # Check for hepatotoxic structural alerts
        alerts = []
        for pattern in self.toxicity_rules['hepatotoxic_patterns']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                alerts.append(pattern)
        
        risk_score = 0
        
        # High lipophilicity increases hepatotoxicity risk
        if logP > 5:
            risk_score += min(35, (logP - 5) * 15)
        
        # Structural alerts
        risk_score += len(alerts) * 20
        
        # Reactive metabolites (simplified)
        if mw > 500:
            risk_score += 15
        
        risk_score = min(100, risk_score)
        
        if risk_score < 25:
            risk_level = 'Low'
            recommendation = 'Low hepatotoxicity risk'
        elif risk_score < 55:
            risk_level = 'Medium'
            recommendation = 'Monitor liver function in preclinical studies'
        else:
            risk_level = 'High'
            recommendation = 'High risk - consider structural modifications'
        
        return {
            'risk_score': round(risk_score, 1),
            'risk_level': risk_level,
            'recommendation': recommendation,
            'logP': round(logP, 2),
            'molecular_weight': round(mw, 2),
            'structural_alerts': len(alerts)
        }
    
    def predict_mutagenicity(self, smiles: str) -> Dict:
        """
        Predict mutagenic potential (Ames test prediction)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Check for mutagenic structural alerts
        alerts = []
        alert_types = []
        for pattern in self.toxicity_rules['mutagenic_patterns']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                alerts.append(pattern)
                if 'N+' in pattern:
                    alert_types.append('Nitro compound')
                elif 'NH2' in pattern:
                    alert_types.append('Aromatic amine')
                elif 'C1OC1' in pattern:
                    alert_types.append('Epoxide')
                elif 'C1NC1' in pattern:
                    alert_types.append('Aziridine')
        
        # Risk based primarily on structural alerts
        risk_score = len(alerts) * 30
        risk_score = min(100, risk_score)
        
        if risk_score == 0:
            risk_level = 'Low'
            recommendation = 'No mutagenic alerts detected'
            prediction = 'Likely non-mutagenic'
        elif risk_score < 50:
            risk_level = 'Medium'
            recommendation = 'Conduct Ames test'
            prediction = 'Possible mutagenic activity'
        else:
            risk_level = 'High'
            recommendation = 'High risk - structural modification strongly recommended'
            prediction = 'Likely mutagenic'
        
        return {
            'risk_score': round(risk_score, 1),
            'risk_level': risk_level,
            'prediction': prediction,
            'recommendation': recommendation,
            'structural_alerts': len(alerts),
            'alert_types': list(set(alert_types))
        }
    
    def predict_carcinogenicity(self, smiles: str) -> Dict:
        """
        Predict carcinogenic potential
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Check for carcinogenic structural alerts
        alerts = []
        alert_types = []
        for pattern in self.toxicity_rules['carcinogenic_patterns']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                alerts.append(pattern)
                if 'NH2' in pattern:
                    alert_types.append('Aromatic amine')
                elif 'N+' in pattern:
                    alert_types.append('N-nitroso compound')
                elif 'ccc' in pattern:
                    alert_types.append('Polycyclic aromatic')
        
        mw = Descriptors.MolWt(mol)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        
        risk_score = len(alerts) * 25
        
        # Polycyclic aromatics increase risk
        if num_aromatic_rings >= 3:
            risk_score += 20
        
        risk_score = min(100, risk_score)
        
        if risk_score == 0:
            risk_level = 'Low'
            recommendation = 'No carcinogenic alerts detected'
            prediction = 'Likely non-carcinogenic'
        elif risk_score < 50:
            risk_level = 'Medium'
            recommendation = 'Conduct carcinogenicity studies'
            prediction = 'Possible carcinogenic activity'
        else:
            risk_level = 'High'
            recommendation = 'High risk - avoid or modify structure'
            prediction = 'Likely carcinogenic'
        
        return {
            'risk_score': round(risk_score, 1),
            'risk_level': risk_level,
            'prediction': prediction,
            'recommendation': recommendation,
            'structural_alerts': len(alerts),
            'alert_types': list(set(alert_types)),
            'aromatic_rings': num_aromatic_rings
        }
    
    def get_comprehensive_profile(self, smiles: str) -> Dict:
        """
        Get complete toxicity profile for a molecule
        """
        return {
            'smiles': smiles,
            'hERG': self.predict_hERG_toxicity(smiles),
            'hepatotoxicity': self.predict_hepatotoxicity(smiles),
            'mutagenicity': self.predict_mutagenicity(smiles),
            'carcinogenicity': self.predict_carcinogenicity(smiles)
        }

def predict_toxicity_profile(smiles: str) -> Dict:
    """
    Main function to get toxicity profile
    """
    profiler = ToxicityProfiler()
    return profiler.get_comprehensive_profile(smiles)

if __name__ == '__main__':
    # Test with aspirin
    test_smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    profile = predict_toxicity_profile(test_smiles)
    print('Toxicity Profile for Aspirin:')
    print(f"hERG Risk: {profile['hERG']['risk_level']} ({profile['hERG']['risk_score']})")
    print(f"Hepatotoxicity Risk: {profile['hepatotoxicity']['risk_level']} ({profile['hepatotoxicity']['risk_score']})")
    print(f"Mutagenicity Risk: {profile['mutagenicity']['risk_level']} ({profile['mutagenicity']['risk_score']})")
    print(f"Carcinogenicity Risk: {profile['carcinogenicity']['risk_level']} ({profile['carcinogenicity']['risk_score']})")
