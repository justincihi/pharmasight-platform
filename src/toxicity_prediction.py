#!/usr/bin/env python3
"""
PharmaSight™ Toxicity and Safety Profiling Module
Phase 4 Implementation: hERG, Hepatotoxicity, Ames Test, CYP450 Inhibition

Implements rule-based predictors using structural alerts and molecular descriptors for:
- hERG channel inhibition (cardiac toxicity risk)
- Hepatotoxicity (liver toxicity)
- Mutagenicity (Ames test)
- CYP450 enzyme inhibition (drug-drug interaction potential)

Note: These predictions use established SAR rules and toxicophore patterns 
(similar to Derek Nexus, OECD QSAR Toolbox) rather than ML models.
"""

import copy
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from datetime import datetime


@dataclass
class ToxicityResult:
    """Container for toxicity prediction result"""
    prediction: str
    probability: float
    confidence: str
    risk_level: str
    details: Dict
    recommendations: List[str]


class MolecularDescriptorCalculator:
    """Calculate molecular descriptors for toxicity prediction"""
    
    @staticmethod
    def calculate_descriptors(mol) -> Dict:
        """Calculate comprehensive molecular descriptors"""
        try:
            descriptors = {
                'mw': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'aliphatic_rings': rdMolDescriptors.CalcNumAliphaticRings(mol),
                'total_rings': rdMolDescriptors.CalcNumRings(mol),
                'fraction_sp3': rdMolDescriptors.CalcFractionCSP3(mol),
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'formal_charge': Chem.GetFormalCharge(mol),
                'num_heteroatoms': rdMolDescriptors.CalcNumHeteroatoms(mol),
                'num_halogens': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']),
                'num_nitrogen': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'),
                'num_oxygen': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'),
                'num_sulfur': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'S'),
                'balaban_j': Descriptors.BalabanJ(mol) if mol.GetNumBonds() > 0 else 0,
                'bertz_ct': Descriptors.BertzCT(mol),
                'chi0': Descriptors.Chi0(mol),
                'chi1': Descriptors.Chi1(mol),
                'kappa1': Descriptors.Kappa1(mol),
                'kappa2': Descriptors.Kappa2(mol),
                'hall_kier_alpha': Descriptors.HallKierAlpha(mol),
                'labute_asa': Descriptors.LabuteASA(mol),
                'peoe_vsa': sum([
                    Descriptors.PEOE_VSA1(mol),
                    Descriptors.PEOE_VSA2(mol),
                    Descriptors.PEOE_VSA3(mol)
                ]),
                'smr_vsa': sum([
                    Descriptors.SMR_VSA1(mol),
                    Descriptors.SMR_VSA2(mol),
                    Descriptors.SMR_VSA3(mol)
                ]),
                'slogp_vsa': sum([
                    Descriptors.SlogP_VSA1(mol),
                    Descriptors.SlogP_VSA2(mol),
                    Descriptors.SlogP_VSA3(mol)
                ]),
                'qed': Descriptors.qed(mol),
                'num_stereocenters': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
            }
            return descriptors
        except Exception as e:
            return {'error': str(e)}
    
    @staticmethod
    def check_substructures(mol) -> Dict:
        """Check for toxicity-related substructures"""
        alerts = {}
        
        toxicophores = {
            'nitro_aromatic': '[N+](=O)[O-]c',
            'nitroso': 'N=O',
            'azide': '[N-]=[N+]=N',
            'diazo': 'N=N',
            'epoxide': 'C1OC1',
            'aziridine': 'C1NC1',
            'alpha_beta_unsaturated_carbonyl': 'C=CC(=O)',
            'polycyclic_aromatic': 'c1ccc2c(c1)cccc2',
            'quinone': 'O=C1C=CC(=O)C=C1',
            'hydrazine': 'NN',
            'hydroxylamine': 'NO',
            'acyl_halide': 'C(=O)[F,Cl,Br,I]',
            'sulfonate_ester': 'OS(=O)(=O)',
            'aromatic_amine': 'c[NH2]',
            'michael_acceptor': 'C=C-C(=O)',
            'aldehyde': '[CH1](=O)',
        }
        
        for name, smarts in toxicophores.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                alerts[name] = mol.HasSubstructMatch(pattern)
        
        return alerts


class HERGPredictor:
    """
    hERG Channel Inhibition Predictor
    Predicts cardiac toxicity risk via hERG potassium channel inhibition
    """
    
    def __init__(self):
        self.threshold_ic50 = 10.0  # µM
        self.high_risk_features = [
            'basic_nitrogen',
            'aromatic_rings_3plus',
            'logp_high',
            'molecular_flexibility'
        ]
    
    def predict(self, smiles: str) -> ToxicityResult:
        """Predict hERG inhibition risk"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ToxicityResult(
                prediction='Error',
                probability=0.0,
                confidence='N/A',
                risk_level='Unknown',
                details={'error': 'Invalid SMILES'},
                recommendations=[]
            )
        
        descriptors = MolecularDescriptorCalculator.calculate_descriptors(mol)
        
        risk_score = self._calculate_herg_risk(mol, descriptors)
        predicted_ic50 = self._estimate_ic50(risk_score)
        
        if risk_score > 0.7:
            prediction = 'hERG+'
            risk_level = 'Critical'
            confidence = 'High' if risk_score > 0.85 else 'Moderate'
        elif risk_score > 0.5:
            prediction = 'Borderline'
            risk_level = 'High'
            confidence = 'Moderate'
        elif risk_score > 0.3:
            prediction = 'Low Risk'
            risk_level = 'Moderate'
            confidence = 'Moderate'
        else:
            prediction = 'hERG-'
            risk_level = 'Low'
            confidence = 'High' if risk_score < 0.15 else 'Moderate'
        
        recommendations = self._generate_recommendations(risk_score, descriptors)
        
        return ToxicityResult(
            prediction=prediction,
            probability=round(risk_score, 3),
            confidence=confidence,
            risk_level=risk_level,
            details={
                'estimated_ic50_um': round(predicted_ic50, 2),
                'threshold_um': self.threshold_ic50,
                'key_risk_factors': self._identify_risk_factors(mol, descriptors),
                'molecular_weight': round(descriptors.get('mw', 0), 1),
                'logp': round(descriptors.get('logp', 0), 2),
                'aromatic_rings': descriptors.get('aromatic_rings', 0),
                'basic_nitrogen_count': self._count_basic_nitrogens(mol)
            },
            recommendations=recommendations
        )
    
    def _calculate_herg_risk(self, mol, descriptors: Dict) -> float:
        """Calculate hERG risk score based on molecular features"""
        risk = 0.0
        
        logp = descriptors.get('logp', 0)
        if logp > 3.5:
            risk += 0.25 * min((logp - 3.5) / 2, 1.0)
        
        mw = descriptors.get('mw', 0)
        if 350 < mw < 600:
            risk += 0.15
        
        aromatic = descriptors.get('aromatic_rings', 0)
        if aromatic >= 3:
            risk += 0.20
        elif aromatic >= 2:
            risk += 0.10
        
        basic_n = self._count_basic_nitrogens(mol)
        if basic_n >= 1:
            risk += 0.20 * min(basic_n, 2)
        
        tpsa = descriptors.get('tpsa', 0)
        if tpsa < 60:
            risk += 0.10
        
        rotatable = descriptors.get('rotatable_bonds', 0)
        if rotatable > 8:
            risk += 0.05
        
        has_piperidine = mol.HasSubstructMatch(Chem.MolFromSmarts('C1CCNCC1'))
        has_piperazine = mol.HasSubstructMatch(Chem.MolFromSmarts('C1CNCCN1'))
        if has_piperidine or has_piperazine:
            risk += 0.15
        
        return min(risk, 1.0)
    
    def _estimate_ic50(self, risk_score: float) -> float:
        """Estimate hERG IC50 from risk score"""
        if risk_score > 0.8:
            return 1.0 + (1 - risk_score) * 5
        elif risk_score > 0.5:
            return 5 + (0.8 - risk_score) * 15
        elif risk_score > 0.3:
            return 20 + (0.5 - risk_score) * 50
        else:
            return 100 + (0.3 - risk_score) * 200
    
    def _count_basic_nitrogens(self, mol) -> int:
        """Count basic nitrogen atoms"""
        basic_n_pattern = Chem.MolFromSmarts('[N;H2,H1,H0;!$(N=*);!$(N#*)]')
        if basic_n_pattern:
            return len(mol.GetSubstructMatches(basic_n_pattern))
        return 0
    
    def _identify_risk_factors(self, mol, descriptors: Dict) -> List[str]:
        """Identify specific hERG risk factors"""
        factors = []
        
        if descriptors.get('logp', 0) > 4:
            factors.append('High lipophilicity (LogP > 4)')
        if descriptors.get('aromatic_rings', 0) >= 3:
            factors.append('Multiple aromatic rings')
        if self._count_basic_nitrogens(mol) >= 2:
            factors.append('Multiple basic nitrogen atoms')
        if descriptors.get('mw', 0) > 500:
            factors.append('High molecular weight')
        if mol.HasSubstructMatch(Chem.MolFromSmarts('C1CCNCC1')):
            factors.append('Contains piperidine ring')
        
        return factors if factors else ['No major risk factors identified']
    
    def _generate_recommendations(self, risk_score: float, descriptors: Dict) -> List[str]:
        """Generate risk mitigation recommendations"""
        recommendations = []
        
        if risk_score > 0.5:
            recommendations.append('Conduct hERG patch clamp assay before advancing')
            
        if descriptors.get('logp', 0) > 4:
            recommendations.append('Consider reducing lipophilicity through polar group addition')
        
        if descriptors.get('aromatic_rings', 0) >= 3:
            recommendations.append('Consider scaffold modifications to reduce planarity')
        
        if risk_score > 0.7:
            recommendations.append('This compound may not be suitable for further development due to cardiac safety concerns')
        elif risk_score > 0.3:
            recommendations.append('Monitor QT interval in preclinical and clinical studies')
        
        return recommendations if recommendations else ['Standard safety monitoring recommended']


class HepatotoxicityPredictor:
    """
    Hepatotoxicity (Drug-Induced Liver Injury) Predictor
    Predicts liver toxicity risk based on molecular features
    """
    
    def __init__(self):
        self.reactive_metabolite_alerts = [
            'quinone',
            'epoxide',
            'alpha_beta_unsaturated_carbonyl',
            'michael_acceptor'
        ]
    
    def predict(self, smiles: str) -> ToxicityResult:
        """Predict hepatotoxicity risk"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ToxicityResult(
                prediction='Error',
                probability=0.0,
                confidence='N/A',
                risk_level='Unknown',
                details={'error': 'Invalid SMILES'},
                recommendations=[]
            )
        
        descriptors = MolecularDescriptorCalculator.calculate_descriptors(mol)
        alerts = MolecularDescriptorCalculator.check_substructures(mol)
        
        risk_score = self._calculate_hepatotox_risk(mol, descriptors, alerts)
        
        if risk_score > 0.7:
            prediction = 'Hepatotoxic'
            risk_level = 'High'
            confidence = 'High' if risk_score > 0.85 else 'Moderate'
        elif risk_score > 0.4:
            prediction = 'Moderate Risk'
            risk_level = 'Moderate'
            confidence = 'Moderate'
        else:
            prediction = 'Low Risk'
            risk_level = 'Low'
            confidence = 'High' if risk_score < 0.2 else 'Moderate'
        
        return ToxicityResult(
            prediction=prediction,
            probability=round(risk_score, 3),
            confidence=confidence,
            risk_level=risk_level,
            details={
                'dili_severity': self._estimate_dili_severity(risk_score),
                'reactive_metabolite_alerts': [a for a, v in alerts.items() if v],
                'lipinski_violations': self._count_lipinski_violations(descriptors),
                'daily_dose_consideration': self._dose_risk_assessment(descriptors),
                'molecular_weight': round(descriptors.get('mw', 0), 1),
                'logp': round(descriptors.get('logp', 0), 2),
            },
            recommendations=self._generate_recommendations(risk_score, alerts, descriptors)
        )
    
    def _calculate_hepatotox_risk(self, mol, descriptors: Dict, alerts: Dict) -> float:
        """Calculate hepatotoxicity risk score"""
        risk = 0.0
        
        reactive_count = sum(1 for a in self.reactive_metabolite_alerts if alerts.get(a, False))
        risk += 0.15 * reactive_count
        
        logp = descriptors.get('logp', 0)
        if logp > 3:
            risk += 0.10 * min((logp - 3) / 2, 1.0)
        
        mw = descriptors.get('mw', 0)
        if mw > 400:
            risk += 0.10
        
        lipinski_violations = self._count_lipinski_violations(descriptors)
        if lipinski_violations >= 2:
            risk += 0.10
        
        if alerts.get('aromatic_amine', False):
            risk += 0.15
        
        if alerts.get('nitro_aromatic', False):
            risk += 0.20
        
        if alerts.get('hydrazine', False) or alerts.get('hydroxylamine', False):
            risk += 0.20
        
        tpsa = descriptors.get('tpsa', 0)
        if tpsa > 140:
            risk += 0.10
        
        return min(risk, 1.0)
    
    def _count_lipinski_violations(self, descriptors: Dict) -> int:
        """Count Lipinski's Rule of Five violations"""
        violations = 0
        if descriptors.get('mw', 0) > 500:
            violations += 1
        if descriptors.get('logp', 0) > 5:
            violations += 1
        if descriptors.get('hbd', 0) > 5:
            violations += 1
        if descriptors.get('hba', 0) > 10:
            violations += 1
        return violations
    
    def _estimate_dili_severity(self, risk_score: float) -> str:
        """Estimate DILI severity classification"""
        if risk_score > 0.7:
            return 'Severe (Most-DILI-Concern)'
        elif risk_score > 0.4:
            return 'Moderate (Less-DILI-Concern)'
        else:
            return 'Low (No-DILI-Concern)'
    
    def _dose_risk_assessment(self, descriptors: Dict) -> str:
        """Assess dose-dependent hepatotoxicity risk"""
        mw = descriptors.get('mw', 0)
        if mw > 600:
            return 'High daily dose compounds (>100 mg/day) show increased DILI risk'
        return 'Monitor if daily dose exceeds 100 mg'
    
    def _generate_recommendations(self, risk_score: float, alerts: Dict, descriptors: Dict) -> List[str]:
        """Generate hepatotoxicity risk mitigation recommendations"""
        recommendations = []
        
        if risk_score > 0.4:
            recommendations.append('Perform reactive metabolite screening')
            recommendations.append('Consider glutathione conjugation studies')
        
        if any(alerts.get(a, False) for a in self.reactive_metabolite_alerts):
            recommendations.append('Structural modification to remove reactive groups recommended')
        
        if risk_score > 0.7:
            recommendations.append('Conduct hepatocyte viability assays')
            recommendations.append('Monitor liver enzymes (ALT, AST) in preclinical studies')
        
        if descriptors.get('logp', 0) > 4:
            recommendations.append('Consider reducing lipophilicity to improve hepatic clearance')
        
        return recommendations if recommendations else ['Standard hepatotoxicity monitoring recommended']


class AmesTestPredictor:
    """
    Ames Test (Mutagenicity) Predictor
    Predicts genotoxicity risk based on structural alerts and molecular features
    """
    
    def __init__(self):
        self.mutagenic_alerts = {
            'aromatic_nitro': ('[N+](=O)[O-]c', 0.4),
            'aromatic_amine': ('c[NH2]', 0.25),
            'polycyclic_aromatic': ('c1ccc2c(c1)cccc2', 0.3),
            'aziridine': ('C1NC1', 0.5),
            'epoxide': ('C1OC1', 0.35),
            'aldehyde': ('[CH1](=O)', 0.2),
            'nitroso': ('N=O', 0.45),
            'azide': ('[N-]=[N+]=N', 0.4),
            'hydrazine': ('NN', 0.35),
            'alkyl_halide': ('[CH2][F,Cl,Br,I]', 0.3),
            'vinyl_halide': ('C=C[F,Cl,Br,I]', 0.35),
        }
    
    def predict(self, smiles: str) -> ToxicityResult:
        """Predict Ames test result (mutagenicity)"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ToxicityResult(
                prediction='Error',
                probability=0.0,
                confidence='N/A',
                risk_level='Unknown',
                details={'error': 'Invalid SMILES'},
                recommendations=[]
            )
        
        descriptors = MolecularDescriptorCalculator.calculate_descriptors(mol)
        
        risk_score, triggered_alerts = self._calculate_mutagenicity_risk(mol)
        
        if risk_score > 0.6:
            prediction = 'Ames+'
            risk_level = 'High'
            confidence = 'High' if risk_score > 0.8 else 'Moderate'
        elif risk_score > 0.3:
            prediction = 'Borderline'
            risk_level = 'Moderate'
            confidence = 'Moderate'
        else:
            prediction = 'Ames-'
            risk_level = 'Low'
            confidence = 'High' if risk_score < 0.15 else 'Moderate'
        
        return ToxicityResult(
            prediction=prediction,
            probability=round(risk_score, 3),
            confidence=confidence,
            risk_level=risk_level,
            details={
                'mutagenic_alerts': triggered_alerts,
                'alert_count': len(triggered_alerts),
                'bacterial_strains_concern': self._identify_strain_concerns(triggered_alerts),
                'metabolic_activation_required': self._needs_s9_activation(mol),
                'molecular_weight': round(descriptors.get('mw', 0), 1),
            },
            recommendations=self._generate_recommendations(risk_score, triggered_alerts)
        )
    
    def _calculate_mutagenicity_risk(self, mol) -> Tuple[float, List[str]]:
        """Calculate mutagenicity risk from structural alerts"""
        total_risk = 0.0
        triggered = []
        
        for alert_name, (smarts, weight) in self.mutagenic_alerts.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                total_risk += weight
                triggered.append(alert_name)
        
        return min(total_risk, 1.0), triggered
    
    def _identify_strain_concerns(self, alerts: List[str]) -> List[str]:
        """Identify Ames strains likely to show positive results"""
        strains = []
        
        if 'aromatic_amine' in alerts or 'polycyclic_aromatic' in alerts:
            strains.extend(['TA98', 'TA100'])
        if 'aromatic_nitro' in alerts:
            strains.extend(['TA98', 'TA1537'])
        if 'epoxide' in alerts or 'aziridine' in alerts:
            strains.extend(['TA100', 'TA1535'])
        if 'aldehyde' in alerts:
            strains.append('TA104')
        
        return list(set(strains)) if strains else ['No specific strain concerns']
    
    def _needs_s9_activation(self, mol) -> str:
        """Determine if metabolic activation is needed"""
        aromatic_amine = mol.HasSubstructMatch(Chem.MolFromSmarts('c[NH2]'))
        polycyclic = mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccc2c(c1)cccc2'))
        
        if aromatic_amine or polycyclic:
            return 'Yes - S9 metabolic activation likely required'
        return 'No - Direct-acting mutagen possible'
    
    def _generate_recommendations(self, risk_score: float, alerts: List[str]) -> List[str]:
        """Generate mutagenicity risk mitigation recommendations"""
        recommendations = []
        
        if risk_score > 0.3:
            recommendations.append('Conduct Ames test with multiple bacterial strains')
            recommendations.append('Include S9 metabolic activation in testing')
        
        if 'aromatic_nitro' in alerts:
            recommendations.append('Consider replacing nitro group with alternative electron-withdrawing group')
        
        if 'aromatic_amine' in alerts:
            recommendations.append('Consider N-alkylation or amide formation to reduce mutagenic potential')
        
        if risk_score > 0.6:
            recommendations.append('Structural modification strongly recommended before advancement')
            recommendations.append('Consider additional genotoxicity assays (micronucleus, chromosomal aberration)')
        
        return recommendations if recommendations else ['Standard genotoxicity screening recommended']


class CYP450InhibitionPredictor:
    """
    CYP450 Enzyme Inhibition Predictor
    Predicts drug-drug interaction potential via CYP450 enzyme inhibition
    """
    
    CYP_ISOFORMS = ['1A2', '2C9', '2C19', '2D6', '3A4']
    
    CYP_INHIBITOR_PATTERNS = {
        '1A2': {
            'patterns': ['c1ccc2[nH]ccc2c1', 'Nc1ccccc1', 'c1ccc2ccccc2c1'],
            'weight_logp': 0.15,
            'weight_aromatic': 0.20,
            'weight_nitrogen': 0.15
        },
        '2C9': {
            'patterns': ['c1ccc(cc1)S(=O)(=O)', 'C(=O)Nc', 'c1ccc(cc1)C(=O)'],
            'weight_logp': 0.20,
            'weight_mw': 0.15,
            'weight_acidic': 0.15
        },
        '2C19': {
            'patterns': ['c1ccc2[nH]c(nc2c1)', 'c1cnc2ccccc2n1'],
            'weight_logp': 0.15,
            'weight_aromatic': 0.15,
            'weight_nitrogen': 0.20
        },
        '2D6': {
            'patterns': ['C1CCNCC1', 'CCNCC', '[NH1]CC'],
            'weight_basic_n': 0.30,
            'weight_aromatic': 0.15,
            'weight_logp': 0.10
        },
        '3A4': {
            'patterns': ['c1ccc2c(c1)cccc2', 'c1ccccc1Cc2ccccc2'],
            'weight_mw': 0.20,
            'weight_logp': 0.20,
            'weight_aromatic': 0.15
        }
    }
    
    def __init__(self):
        self.threshold_ic50 = 10.0  # µM threshold for significant inhibition
    
    def predict(self, smiles: str) -> Dict:
        """Predict CYP450 inhibition profile"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        descriptors = MolecularDescriptorCalculator.calculate_descriptors(mol)
        
        inhibition_profile = {}
        overall_ddi_risk = 0.0
        
        for isoform in self.CYP_ISOFORMS:
            result = self._predict_isoform_inhibition(mol, isoform, descriptors)
            inhibition_profile[f'CYP{isoform}'] = result
            overall_ddi_risk = max(overall_ddi_risk, result['probability'])
        
        phenotyping = self._identify_phenotyping_needs(inhibition_profile)
        clinical_relevance = self._assess_clinical_relevance(inhibition_profile)
        
        return {
            'smiles': smiles,
            'inhibition_profile': inhibition_profile,
            'overall_ddi_risk': round(overall_ddi_risk, 3),
            'ddi_risk_level': self._categorize_ddi_risk(overall_ddi_risk),
            'primary_cyp_concerns': self._identify_primary_concerns(inhibition_profile),
            'phenotyping_recommendations': phenotyping,
            'clinical_relevance': clinical_relevance,
            'victim_drug_examples': self._get_victim_drug_examples(inhibition_profile),
            'recommendations': self._generate_recommendations(inhibition_profile, overall_ddi_risk)
        }
    
    def _predict_isoform_inhibition(self, mol, isoform: str, descriptors: Dict) -> Dict:
        """Predict inhibition for a specific CYP isoform"""
        config = self.CYP_INHIBITOR_PATTERNS.get(isoform, {})
        
        risk = 0.0
        matched_patterns = []
        
        for pattern_smarts in config.get('patterns', []):
            pattern = Chem.MolFromSmarts(pattern_smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                risk += 0.15
                matched_patterns.append(pattern_smarts)
        
        logp = descriptors.get('logp', 0)
        if logp > 2 and config.get('weight_logp', 0) > 0:
            risk += config['weight_logp'] * min((logp - 2) / 3, 1.0)
        
        aromatic = descriptors.get('aromatic_rings', 0)
        if aromatic >= 2 and config.get('weight_aromatic', 0) > 0:
            risk += config['weight_aromatic']
        
        mw = descriptors.get('mw', 0)
        if mw > 400 and config.get('weight_mw', 0) > 0:
            risk += config['weight_mw']
        
        risk = min(risk, 1.0)
        
        if risk > 0.6:
            prediction = 'Inhibitor'
            ic50_estimate = 1 + (1 - risk) * 9
        elif risk > 0.3:
            prediction = 'Weak Inhibitor'
            ic50_estimate = 10 + (0.6 - risk) * 40
        else:
            prediction = 'Non-Inhibitor'
            ic50_estimate = 50 + (0.3 - risk) * 150
        
        return {
            'prediction': prediction,
            'probability': round(risk, 3),
            'estimated_ic50_um': round(ic50_estimate, 1),
            'matched_patterns': len(matched_patterns),
            'confidence': 'High' if risk > 0.7 or risk < 0.2 else 'Moderate'
        }
    
    def _categorize_ddi_risk(self, risk: float) -> str:
        """Categorize overall DDI risk level"""
        if risk > 0.7:
            return 'High'
        elif risk > 0.4:
            return 'Moderate'
        elif risk > 0.2:
            return 'Low'
        else:
            return 'Minimal'
    
    def _identify_primary_concerns(self, profile: Dict) -> List[str]:
        """Identify CYPs with highest inhibition risk"""
        concerns = []
        for isoform, data in profile.items():
            if data.get('probability', 0) > 0.5:
                concerns.append(f"{isoform}: {data['prediction']}")
        return concerns if concerns else ['No major CYP inhibition concerns']
    
    def _identify_phenotyping_needs(self, profile: Dict) -> List[str]:
        """Identify CYP phenotyping recommendations"""
        needs = []
        for isoform, data in profile.items():
            if data.get('probability', 0) > 0.4:
                needs.append(f"Include {isoform} substrates in DDI studies")
        return needs if needs else ['Standard CYP panel screening recommended']
    
    def _assess_clinical_relevance(self, profile: Dict) -> Dict:
        """Assess clinical relevance of inhibition profile"""
        high_risk_cyps = [iso for iso, data in profile.items() if data.get('probability', 0) > 0.6]
        
        relevance = {
            'polypharmacy_concern': len(high_risk_cyps) >= 2,
            'narrow_therapeutic_index_caution': any(
                iso in ['CYP2C9', 'CYP2D6', 'CYP3A4'] for iso in high_risk_cyps
            ),
            'genetic_polymorphism_relevance': any(
                iso in ['CYP2C9', 'CYP2C19', 'CYP2D6'] for iso in high_risk_cyps
            )
        }
        
        return relevance
    
    def _get_victim_drug_examples(self, profile: Dict) -> Dict:
        """Get examples of victim drugs for inhibited CYPs"""
        victim_drugs = {
            'CYP1A2': ['Caffeine', 'Theophylline', 'Clozapine', 'Melatonin'],
            'CYP2C9': ['Warfarin', 'Phenytoin', 'NSAIDs', 'Losartan'],
            'CYP2C19': ['Omeprazole', 'Clopidogrel', 'Diazepam', 'Citalopram'],
            'CYP2D6': ['Codeine', 'Tamoxifen', 'Metoprolol', 'Fluoxetine'],
            'CYP3A4': ['Midazolam', 'Cyclosporine', 'Statins', 'Many HIV drugs']
        }
        
        result = {}
        for isoform, data in profile.items():
            if data.get('probability', 0) > 0.4:
                result[isoform] = victim_drugs.get(isoform, [])
        
        return result if result else {'note': 'No significant victim drug concerns'}
    
    def _generate_recommendations(self, profile: Dict, overall_risk: float) -> List[str]:
        """Generate DDI risk mitigation recommendations"""
        recommendations = []
        
        if overall_risk > 0.5:
            recommendations.append('Conduct in vitro CYP inhibition studies with recombinant enzymes')
            recommendations.append('Consider clinical DDI study with probe substrates')
        
        high_risk_cyps = [iso for iso, data in profile.items() if data.get('probability', 0) > 0.6]
        
        if 'CYP3A4' in high_risk_cyps:
            recommendations.append('Avoid co-administration with CYP3A4-metabolized narrow TI drugs')
        
        if 'CYP2D6' in high_risk_cyps:
            recommendations.append('Consider CYP2D6 genotype-guided dosing in clinical trials')
        
        if overall_risk > 0.7:
            recommendations.append('Consider structural modifications to reduce CYP inhibition liability')
        
        return recommendations if recommendations else ['Standard DDI assessment in development']


class ToxicityProfiler:
    """
    Comprehensive Toxicity Profiler
    Combines all toxicity predictors for a complete safety assessment
    """
    
    def __init__(self):
        self.herg_predictor = HERGPredictor()
        self.hepatotox_predictor = HepatotoxicityPredictor()
        self.ames_predictor = AmesTestPredictor()
        self.cyp_predictor = CYP450InhibitionPredictor()
    
    def get_full_profile(self, smiles: str) -> Dict:
        """Get comprehensive toxicity profile"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        herg_result = self.herg_predictor.predict(smiles)
        hepatotox_result = self.hepatotox_predictor.predict(smiles)
        ames_result = self.ames_predictor.predict(smiles)
        cyp_result = self.cyp_predictor.predict(smiles)
        
        overall_risk = self._calculate_overall_toxicity_risk(
            herg_result, hepatotox_result, ames_result, cyp_result
        )
        
        return {
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(mol, canonical=True),
            'timestamp': datetime.now().isoformat(),
            'herg_assessment': {
                'prediction': herg_result.prediction,
                'probability': herg_result.probability,
                'risk_level': herg_result.risk_level,
                'details': herg_result.details,
                'recommendations': herg_result.recommendations
            },
            'hepatotoxicity_assessment': {
                'prediction': hepatotox_result.prediction,
                'probability': hepatotox_result.probability,
                'risk_level': hepatotox_result.risk_level,
                'details': hepatotox_result.details,
                'recommendations': hepatotox_result.recommendations
            },
            'mutagenicity_assessment': {
                'prediction': ames_result.prediction,
                'probability': ames_result.probability,
                'risk_level': ames_result.risk_level,
                'details': ames_result.details,
                'recommendations': ames_result.recommendations
            },
            'cyp450_assessment': cyp_result,
            'overall_toxicity_score': round(overall_risk, 3),
            'overall_risk_category': self._categorize_overall_risk(overall_risk),
            'development_recommendation': self._get_development_recommendation(overall_risk),
            'priority_concerns': self._identify_priority_concerns(
                herg_result, hepatotox_result, ames_result, cyp_result
            )
        }
    
    def _calculate_overall_toxicity_risk(
        self, herg, hepatotox, ames, cyp
    ) -> float:
        """Calculate weighted overall toxicity risk"""
        weights = {
            'herg': 0.30,
            'hepatotox': 0.25,
            'ames': 0.25,
            'cyp': 0.20
        }
        
        total = (
            weights['herg'] * herg.probability +
            weights['hepatotox'] * hepatotox.probability +
            weights['ames'] * ames.probability +
            weights['cyp'] * cyp.get('overall_ddi_risk', 0)
        )
        
        return min(total, 1.0)
    
    def _categorize_overall_risk(self, risk: float) -> str:
        """Categorize overall toxicity risk"""
        if risk > 0.7:
            return 'Critical'
        elif risk > 0.5:
            return 'High'
        elif risk > 0.3:
            return 'Moderate'
        elif risk > 0.15:
            return 'Low'
        else:
            return 'Minimal'
    
    def _get_development_recommendation(self, risk: float) -> str:
        """Get development stage recommendation"""
        if risk > 0.7:
            return 'Not recommended for advancement - significant structural modifications required'
        elif risk > 0.5:
            return 'Proceed with caution - address major safety liabilities before advancement'
        elif risk > 0.3:
            return 'Acceptable for early development with comprehensive safety monitoring'
        else:
            return 'Favorable safety profile - suitable for advancement with standard monitoring'
    
    def _identify_priority_concerns(
        self, herg, hepatotox, ames, cyp
    ) -> List[Dict]:
        """Identify and prioritize safety concerns"""
        concerns = []
        
        if herg.probability > 0.5:
            concerns.append({
                'category': 'Cardiac Safety',
                'concern': 'hERG inhibition',
                'severity': herg.risk_level,
                'probability': herg.probability
            })
        
        if hepatotox.probability > 0.5:
            concerns.append({
                'category': 'Hepatic Safety',
                'concern': 'Hepatotoxicity risk',
                'severity': hepatotox.risk_level,
                'probability': hepatotox.probability
            })
        
        if ames.probability > 0.5:
            concerns.append({
                'category': 'Genetic Toxicity',
                'concern': 'Mutagenicity',
                'severity': ames.risk_level,
                'probability': ames.probability
            })
        
        cyp_risk = cyp.get('overall_ddi_risk', 0)
        if cyp_risk > 0.5:
            concerns.append({
                'category': 'Drug Interactions',
                'concern': 'CYP450 inhibition',
                'severity': cyp.get('ddi_risk_level', 'Unknown'),
                'probability': cyp_risk
            })
        
        concerns.sort(key=lambda x: x['probability'], reverse=True)
        return concerns if concerns else [{'category': 'None', 'concern': 'No major safety concerns identified'}]


def get_toxicity_predictor():
    """Get singleton toxicity profiler instance"""
    return ToxicityProfiler()
