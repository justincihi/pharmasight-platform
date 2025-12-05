"""
Enhanced Multi-Receptor Docking Scorer
Provides advanced docking scoring with receptor family affinity profiling,
therapeutic relevance weighting, and integration with viability analysis.
Includes performance caching for repeated calculations.
"""

import copy
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Fragments
from rdkit import DataStructs
import json
import hashlib
from functools import lru_cache
from datetime import datetime, timedelta
import threading


@dataclass
class ReceptorProfile:
    """Receptor family affinity profile"""
    family: str
    affinity_score: float
    confidence: float
    predicted_ki: str
    mechanism: str
    therapeutic_relevance: List[str]
    safety_concerns: List[str]


THERAPEUTIC_WEIGHTS = {
    'Psychiatric': 1.2,
    'Neurological': 1.15,
    'Pain': 1.1,
    'Cardiovascular': 1.0,
    'Metabolic': 1.0,
    'Immunological': 0.95,
    'Oncology': 1.25,
    'Infectious': 1.1,
}

RECEPTOR_THERAPEUTIC_MAP = {
    'Serotonin': ['Psychiatric', 'Neurological', 'Pain'],
    'Dopamine': ['Psychiatric', 'Neurological', 'Pain'],
    'Opioid': ['Pain', 'Psychiatric'],
    'GABA': ['Psychiatric', 'Neurological', 'Pain'],
    'Glutamate': ['Neurological', 'Psychiatric', 'Pain'],
    'Adrenergic': ['Cardiovascular', 'Psychiatric', 'Metabolic'],
    'Cannabinoid': ['Pain', 'Psychiatric', 'Neurological'],
    'Muscarinic': ['Neurological', 'Cardiovascular'],
    'Nicotinic': ['Neurological', 'Psychiatric'],
    'Histamine': ['Immunological', 'Neurological'],
    'Adenosine': ['Cardiovascular', 'Neurological'],
    'TRP': ['Pain', 'Neurological'],
    'Purinergic': ['Pain', 'Cardiovascular'],
    'Orexin': ['Psychiatric', 'Metabolic'],
    'Melatonin': ['Psychiatric', 'Neurological'],
    'Chemokine': ['Immunological', 'Oncology'],
}

SAFETY_RECEPTORS = {
    '5-HT2B': {'concern': 'Cardiac valvulopathy', 'severity': 'Critical'},
    'hERG': {'concern': 'QT prolongation', 'severity': 'Critical'},
    'D2': {'concern': 'Extrapyramidal effects', 'severity': 'High'},
    'H1': {'concern': 'Sedation/drowsiness', 'severity': 'Moderate'},
    'M1': {'concern': 'Cognitive impairment', 'severity': 'Moderate'},
    'Alpha1A': {'concern': 'Orthostatic hypotension', 'severity': 'Moderate'},
    'Sigma1': {'concern': 'Psychotomimetic effects', 'severity': 'High'},
}


class DockingCache:
    """Thread-safe cache for docking results with TTL"""
    
    def __init__(self, max_size: int = 1000, ttl_minutes: int = 60):
        self._cache = {}
        self._timestamps = {}
        self._lock = threading.Lock()
        self.max_size = max_size
        self.ttl = timedelta(minutes=ttl_minutes)
        self._hits = 0
        self._misses = 0
    
    def _make_key(self, smiles: str, target_families: List[str] = None, include_safety: bool = True) -> str:
        """Generate cache key from inputs"""
        families_str = ','.join(sorted(target_families)) if target_families else 'all'
        key_str = f"{smiles}|{families_str}|{include_safety}"
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, smiles: str, target_families: List[str] = None, include_safety: bool = True) -> Optional[Dict]:
        """Get cached result if valid (deep copy to prevent mutation leaks)"""
        key = self._make_key(smiles, target_families, include_safety)
        with self._lock:
            if key in self._cache:
                if datetime.now() - self._timestamps[key] < self.ttl:
                    self._hits += 1
                    return copy.deepcopy(self._cache[key])
                else:
                    del self._cache[key]
                    del self._timestamps[key]
            self._misses += 1
            return None
    
    def set(self, smiles: str, result: Dict, target_families: List[str] = None, include_safety: bool = True):
        """Cache a result (deep copy to isolate from caller mutations)"""
        key = self._make_key(smiles, target_families, include_safety)
        with self._lock:
            if len(self._cache) >= self.max_size:
                oldest_key = min(self._timestamps, key=self._timestamps.get)
                del self._cache[oldest_key]
                del self._timestamps[oldest_key]
            
            self._cache[key] = copy.deepcopy(result)
            self._timestamps[key] = datetime.now()
    
    def clear(self):
        """Clear all cached results"""
        with self._lock:
            self._cache.clear()
            self._timestamps.clear()
    
    def stats(self) -> Dict:
        """Get cache statistics"""
        with self._lock:
            total = self._hits + self._misses
            hit_rate = (self._hits / total * 100) if total > 0 else 0
            return {
                'size': len(self._cache),
                'max_size': self.max_size,
                'hits': self._hits,
                'misses': self._misses,
                'hit_rate': f"{hit_rate:.1f}%"
            }


_docking_cache = DockingCache(max_size=1000, ttl_minutes=60)


class EnhancedDockingScorer:
    """Advanced multi-receptor docking scorer with therapeutic weighting and caching"""
    
    def __init__(self, use_cache: bool = True):
        self.receptor_profiles = {}
        self.pharmacophore_weights = self._init_pharmacophore_weights()
        self.use_cache = use_cache
        self.cache = _docking_cache
        
    def _init_pharmacophore_weights(self) -> Dict:
        """Initialize pharmacophore feature weights for each receptor family"""
        return {
            'Serotonin': {
                'indole': 0.30, 'basic_nitrogen': 0.25, 'aromatic': 0.15,
                'logp_optimal': (1.5, 4.0), 'mw_optimal': (180, 400),
                'tpsa_optimal': (20, 80), 'rotatable_max': 8
            },
            'Dopamine': {
                'catechol': 0.25, 'basic_nitrogen': 0.30, 'phenol': 0.15,
                'logp_optimal': (1.0, 4.5), 'mw_optimal': (150, 450),
                'tpsa_optimal': (30, 90), 'rotatable_max': 10
            },
            'Opioid': {
                'piperidine': 0.25, 'basic_nitrogen': 0.25, 'rings_min': 3,
                'logp_optimal': (1.5, 5.0), 'mw_optimal': (250, 550),
                'tpsa_optimal': (40, 100), 'rotatable_max': 8
            },
            'GABA': {
                'amide': 0.20, 'halogen': 0.15, 'hbd_min': 2,
                'logp_optimal': (-1, 3), 'mw_optimal': (150, 400),
                'tpsa_optimal': (50, 120), 'rotatable_max': 6
            },
            'Cannabinoid': {
                'phenol': 0.20, 'alkyl_chain': 0.25, 'aromatic': 0.15,
                'logp_optimal': (3.0, 7.0), 'mw_optimal': (280, 500),
                'tpsa_optimal': (40, 80), 'rotatable_max': 12
            },
            'Adrenergic': {
                'catechol': 0.25, 'basic_nitrogen': 0.25, 'phenol': 0.20,
                'logp_optimal': (0.5, 3.5), 'mw_optimal': (150, 400),
                'tpsa_optimal': (40, 100), 'rotatable_max': 8
            },
        }
    
    def calculate_enhanced_docking_score(
        self, 
        smiles: str,
        target_families: List[str] = None,
        include_safety: bool = True
    ) -> Dict:
        """
        Calculate enhanced docking scores for a compound with caching
        
        Args:
            smiles: SMILES string
            target_families: Optional list of receptor families to focus on
            include_safety: Include safety receptor screening
            
        Returns:
            Dict with docking scores, profiles, and recommendations
        """
        if self.use_cache:
            cached_result = self.cache.get(smiles, target_families, include_safety)
            if cached_result:
                cached_result['from_cache'] = True
                return cached_result
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        features = self._extract_molecular_features(mol)
        
        receptor_scores = {}
        for family in RECEPTOR_THERAPEUTIC_MAP.keys():
            if target_families and family not in target_families:
                continue
            score_data = self._score_receptor_family(mol, family, features)
            receptor_scores[family] = score_data
        
        affinity_profile = self._generate_affinity_profile(receptor_scores)
        
        safety_profile = {}
        if include_safety:
            safety_profile = self._assess_safety_targets(mol, features)
        
        therapeutic_scores = self._calculate_therapeutic_scores(receptor_scores)
        
        lead_potential = self._assess_lead_potential(
            receptor_scores, safety_profile, features
        )
        
        result = {
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(mol, canonical=True),
            'molecular_features': features,
            'receptor_family_scores': receptor_scores,
            'affinity_profile': affinity_profile,
            'safety_assessment': safety_profile,
            'therapeutic_potential': therapeutic_scores,
            'lead_potential': lead_potential,
            'top_targets': self._get_top_targets(receptor_scores, 5),
            'recommendations': self._generate_recommendations(
                receptor_scores, safety_profile, therapeutic_scores
            ),
            'from_cache': False
        }
        
        if self.use_cache:
            self.cache.set(smiles, result, target_families, include_safety)
        
        return result
    
    def get_cache_stats(self) -> Dict:
        """Get cache performance statistics"""
        return self.cache.stats()
    
    def clear_cache(self):
        """Clear the docking results cache"""
        self.cache.clear()
    
    def _extract_molecular_features(self, mol) -> Dict:
        """Extract comprehensive molecular features for scoring"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        aliphatic_rings = rdMolDescriptors.CalcNumAliphaticRings(mol)
        total_rings = rdMolDescriptors.CalcNumRings(mol)
        fraction_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        
        has_basic_nitrogen = self._check_basic_nitrogen(mol)
        has_indole = mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccc2[nH]ccc2c1')) if Chem.MolFromSmarts('c1ccc2[nH]ccc2c1') else False
        has_phenol = mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccccc1O')) if Chem.MolFromSmarts('c1ccccc1O') else False
        has_catechol = mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccc(O)c(O)c1')) if Chem.MolFromSmarts('c1ccc(O)c(O)c1') else False
        has_piperidine = mol.HasSubstructMatch(Chem.MolFromSmarts('C1CCNCC1')) if Chem.MolFromSmarts('C1CCNCC1') else False
        has_morpholine = mol.HasSubstructMatch(Chem.MolFromSmarts('C1COCCN1')) if Chem.MolFromSmarts('C1COCCN1') else False
        has_amide = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)N')) if Chem.MolFromSmarts('C(=O)N') else False
        has_sulfonamide = mol.HasSubstructMatch(Chem.MolFromSmarts('S(=O)(=O)N')) if Chem.MolFromSmarts('S(=O)(=O)N') else False
        has_halogen = any(mol.HasSubstructMatch(Chem.MolFromSmarts(f'[{x}]')) for x in ['F', 'Cl', 'Br', 'I'])
        
        num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        num_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        
        return {
            'molecular_weight': round(mw, 2),
            'logP': round(logp, 2),
            'hbd': hbd,
            'hba': hba,
            'tpsa': round(tpsa, 2),
            'rotatable_bonds': rotatable,
            'aromatic_rings': aromatic_rings,
            'aliphatic_rings': aliphatic_rings,
            'total_rings': total_rings,
            'fraction_sp3': round(fraction_sp3, 3),
            'num_heteroatoms': num_heteroatoms,
            'num_nitrogen': num_nitrogen,
            'num_oxygen': num_oxygen,
            'has_basic_nitrogen': has_basic_nitrogen,
            'has_indole': has_indole,
            'has_phenol': has_phenol,
            'has_catechol': has_catechol,
            'has_piperidine': has_piperidine,
            'has_morpholine': has_morpholine,
            'has_amide': has_amide,
            'has_sulfonamide': has_sulfonamide,
            'has_halogen': has_halogen,
        }
    
    def _check_basic_nitrogen(self, mol) -> bool:
        """Check for basic (protonatable) nitrogen"""
        patterns = [
            '[N;X3;!$(N-C=O);!$(N-S=O)]',
            '[NH1;!$(N-C=O);!$(N-S=O)]',
            '[NH2;!$(N-C=O)]'
        ]
        for pattern in patterns:
            patt = Chem.MolFromSmarts(pattern)
            if patt and mol.HasSubstructMatch(patt):
                return True
        return False
    
    def _score_receptor_family(self, mol, family: str, features: Dict) -> Dict:
        """Score compound against a receptor family"""
        base_score = 0.20
        confidence = 0.5
        
        mw = features['molecular_weight']
        logp = features['logP']
        tpsa = features['tpsa']
        
        if family == 'Serotonin':
            if features['has_indole']:
                base_score += 0.35
                confidence += 0.25
            elif features['aromatic_rings'] >= 2 and features['has_basic_nitrogen']:
                base_score += 0.20
                confidence += 0.15
            if features['has_basic_nitrogen'] and 1.5 <= logp <= 4.0:
                base_score += 0.15
                confidence += 0.1
            if 180 < mw < 400:
                base_score += 0.10
            if 20 < tpsa < 80:
                base_score += 0.05
                
        elif family == 'Dopamine':
            if features['has_catechol'] and features['has_basic_nitrogen']:
                base_score += 0.35
                confidence += 0.25
            elif features['has_phenol'] and features['has_basic_nitrogen']:
                base_score += 0.25
                confidence += 0.15
            elif features['aromatic_rings'] >= 1 and features['has_basic_nitrogen']:
                base_score += 0.15
            if 1.0 <= logp <= 4.5:
                base_score += 0.10
            if 150 < mw < 450:
                base_score += 0.08
            if features['has_piperidine'] or features['has_morpholine']:
                base_score += 0.10
                confidence += 0.1
                
        elif family == 'Opioid':
            if features['total_rings'] >= 3 and features['has_basic_nitrogen']:
                base_score += 0.30
                confidence += 0.2
            elif features['has_piperidine'] and features['aromatic_rings'] >= 1:
                base_score += 0.25
                confidence += 0.15
            if 250 < mw < 550:
                base_score += 0.12
            if 1.5 <= logp <= 5.0:
                base_score += 0.10
            if features['hbd'] >= 1 and features['hba'] >= 2:
                base_score += 0.08
                
        elif family == 'GABA':
            if features['has_amide'] or (features['hbd'] >= 2 and features['hba'] >= 3):
                base_score += 0.28
                confidence += 0.2
            if 150 < mw < 400:
                base_score += 0.12
            if -1 <= logp <= 3:
                base_score += 0.12
            if 50 < tpsa < 120:
                base_score += 0.10
            if features['has_halogen'] and features['aromatic_rings'] >= 1:
                base_score += 0.10
                confidence += 0.1
                
        elif family == 'Glutamate':
            if features['hbd'] >= 2 and features['hba'] >= 4:
                base_score += 0.25
                confidence += 0.2
            if 100 < mw < 350:
                base_score += 0.15
            if -2 <= logp <= 2:
                base_score += 0.15
            if tpsa > 80:
                base_score += 0.10
                
        elif family == 'Cannabinoid':
            if features['aromatic_rings'] >= 1 and features['aliphatic_rings'] >= 1:
                base_score += 0.22
                confidence += 0.15
            if 3 <= logp <= 7:
                base_score += 0.18
                confidence += 0.1
            if 280 < mw < 500:
                base_score += 0.12
            if features['has_phenol']:
                base_score += 0.10
            if features['rotatable_bonds'] >= 4:
                base_score += 0.08
                
        elif family == 'Adrenergic':
            if features['has_catechol'] and features['has_basic_nitrogen']:
                base_score += 0.35
                confidence += 0.25
            elif features['has_phenol'] and features['has_basic_nitrogen'] and features['hbd'] >= 2:
                base_score += 0.28
                confidence += 0.2
            if 0.5 <= logp <= 3.5:
                base_score += 0.12
            if 150 < mw < 400:
                base_score += 0.10
                
        elif family == 'Muscarinic':
            if features['has_basic_nitrogen'] and features['total_rings'] >= 2:
                base_score += 0.28
                confidence += 0.2
            if 200 < mw < 450:
                base_score += 0.12
            if 1 <= logp <= 4:
                base_score += 0.12
            if features['has_piperidine'] or features['aliphatic_rings'] >= 1:
                base_score += 0.08
                
        elif family == 'Nicotinic':
            pyridine_patt = Chem.MolFromSmarts('c1ccncc1')
            has_pyridine = pyridine_patt and mol.HasSubstructMatch(pyridine_patt)
            if has_pyridine and features['has_basic_nitrogen']:
                base_score += 0.32
                confidence += 0.25
            elif features['aromatic_rings'] >= 1 and features['has_basic_nitrogen']:
                base_score += 0.18
            if 120 < mw < 350:
                base_score += 0.12
            if 0 <= logp <= 3:
                base_score += 0.12
                
        elif family == 'Histamine':
            if features['has_basic_nitrogen'] and features['aromatic_rings'] >= 1:
                base_score += 0.22
                confidence += 0.15
            if 200 < mw < 400:
                base_score += 0.15
            if 1 <= logp <= 4:
                base_score += 0.12
            if features['has_piperidine'] or features['has_morpholine']:
                base_score += 0.08
                
        elif family == 'Adenosine':
            if features['total_rings'] >= 2 and features['num_heteroatoms'] >= 5:
                base_score += 0.30
                confidence += 0.25
            if 220 < mw < 420:
                base_score += 0.12
            if -1 <= logp <= 2.5:
                base_score += 0.10
            if features['hba'] >= 5:
                base_score += 0.08
                
        elif family == 'TRP':
            if features['aromatic_rings'] >= 1 and 2.5 <= logp <= 5.5:
                base_score += 0.22
                confidence += 0.15
            if 280 < mw < 480:
                base_score += 0.12
            if features['hbd'] >= 1 and features['hba'] >= 2 and features['has_amide']:
                base_score += 0.15
                confidence += 0.1
                
        elif family == 'Purinergic':
            if features['num_heteroatoms'] >= 5 and features['total_rings'] >= 2:
                base_score += 0.30
                confidence += 0.2
            if 220 < mw < 480:
                base_score += 0.10
            if -1 <= logp <= 2.5:
                base_score += 0.12
            if tpsa > 90:
                base_score += 0.10
                
        elif family == 'Orexin':
            if features['has_amide'] and features['aromatic_rings'] >= 2 and features['total_rings'] >= 3:
                base_score += 0.30
                confidence += 0.25
            if 380 < mw < 580:
                base_score += 0.15
            if 2.5 <= logp <= 5:
                base_score += 0.10
                
        elif family == 'Melatonin':
            if features['has_indole'] and features['has_amide']:
                base_score += 0.40
                confidence += 0.3
            elif features['has_indole']:
                base_score += 0.28
                confidence += 0.2
            if 200 < mw < 320:
                base_score += 0.10
            if 1 <= logp <= 3:
                base_score += 0.08
                
        elif family == 'Chemokine':
            if features['has_basic_nitrogen'] and features['aromatic_rings'] >= 2 and features['total_rings'] >= 3:
                base_score += 0.28
                confidence += 0.2
            if 380 < mw < 600:
                base_score += 0.12
            if 2.5 <= logp <= 5.5:
                base_score += 0.10
        
        else:
            if features['aromatic_rings'] >= 1 and features['has_basic_nitrogen']:
                base_score += 0.15
            if 200 < mw < 500:
                base_score += 0.10
            if 0 <= logp <= 5:
                base_score += 0.10
        
        final_score = max(0.0, min(1.0, base_score))
        confidence = max(0.0, min(1.0, confidence))
        
        return {
            'score': round(final_score, 3),
            'confidence': round(confidence, 3),
            'predicted_ki': self._score_to_ki(final_score),
            'affinity_class': self._classify_affinity(final_score),
            'therapeutic_areas': RECEPTOR_THERAPEUTIC_MAP.get(family, [])
        }
    
    def _score_to_ki(self, score: float) -> str:
        """Convert score to predicted Ki range"""
        if score > 0.9:
            return "<1 nM"
        elif score > 0.85:
            return "1-10 nM"
        elif score > 0.75:
            return "10-100 nM"
        elif score > 0.65:
            return "100-500 nM"
        elif score > 0.55:
            return "0.5-5 µM"
        elif score > 0.45:
            return "5-50 µM"
        else:
            return ">50 µM"
    
    def _classify_affinity(self, score: float) -> str:
        """Classify affinity level"""
        if score > 0.85:
            return "Very High"
        elif score > 0.70:
            return "High"
        elif score > 0.55:
            return "Moderate"
        elif score > 0.40:
            return "Low"
        else:
            return "Very Low"
    
    def _generate_affinity_profile(self, receptor_scores: Dict) -> Dict:
        """Generate overall affinity profile"""
        if not receptor_scores:
            return {'status': 'No scores available'}
        
        scores = [data['score'] for data in receptor_scores.values()]
        high_affinity = [f for f, data in receptor_scores.items() if data['score'] > 0.7]
        moderate_affinity = [f for f, data in receptor_scores.items() if 0.5 <= data['score'] <= 0.7]
        
        selectivity_index = 1.0
        if len(high_affinity) == 1:
            selectivity_index = 10.0
        elif len(high_affinity) == 2:
            selectivity_index = 5.0
        elif len(high_affinity) >= 3:
            selectivity_index = 2.0
        
        return {
            'max_score': round(max(scores), 3),
            'mean_score': round(np.mean(scores), 3),
            'std_score': round(np.std(scores), 3),
            'high_affinity_targets': high_affinity,
            'moderate_affinity_targets': moderate_affinity,
            'selectivity_index': selectivity_index,
            'selectivity_class': 'Selective' if selectivity_index >= 5 else 'Moderate' if selectivity_index >= 2 else 'Promiscuous',
            'polypharmacology_score': len(high_affinity) * 15 + len(moderate_affinity) * 5
        }
    
    def _assess_safety_targets(self, mol, features: Dict) -> Dict:
        """Assess binding to safety-critical receptors"""
        safety_concerns = []
        overall_risk = 'Low'
        
        for receptor, info in SAFETY_RECEPTORS.items():
            risk_score = self._estimate_safety_receptor_binding(mol, receptor, features)
            if risk_score > 0.6:
                severity = info['severity']
                safety_concerns.append({
                    'receptor': receptor,
                    'concern': info['concern'],
                    'severity': severity,
                    'estimated_binding': round(risk_score, 3)
                })
                if severity == 'Critical':
                    overall_risk = 'Critical'
                elif severity == 'High' and overall_risk not in ['Critical']:
                    overall_risk = 'High'
                elif severity == 'Moderate' and overall_risk == 'Low':
                    overall_risk = 'Moderate'
        
        return {
            'concerns': safety_concerns,
            'overall_risk': overall_risk,
            'recommendation': self._get_safety_recommendation(overall_risk, safety_concerns)
        }
    
    def _estimate_safety_receptor_binding(self, mol, receptor: str, features: Dict) -> float:
        """Estimate binding to specific safety receptors"""
        score = 0.2
        
        if receptor == '5-HT2B':
            if features['has_indole'] and features['has_basic_nitrogen']:
                score += 0.4
            elif features['aromatic_rings'] >= 2 and features['has_basic_nitrogen']:
                score += 0.25
            if 200 < features['molecular_weight'] < 400:
                score += 0.1
                
        elif receptor == 'hERG':
            if features['has_basic_nitrogen'] and features['aromatic_rings'] >= 2:
                score += 0.3
            if features['logP'] > 3.5:
                score += 0.15
            if features['molecular_weight'] > 400:
                score += 0.1
                
        elif receptor == 'D2':
            if features['has_basic_nitrogen'] and features['aromatic_rings'] >= 1:
                score += 0.25
            if features['has_piperidine']:
                score += 0.15
                
        elif receptor == 'H1':
            if features['has_basic_nitrogen'] and features['aromatic_rings'] >= 1:
                score += 0.25
            if 2 <= features['logP'] <= 4:
                score += 0.1
                
        elif receptor == 'M1':
            if features['has_basic_nitrogen'] and features['total_rings'] >= 2:
                score += 0.25
                
        elif receptor == 'Alpha1A':
            if features['has_basic_nitrogen'] and features['has_phenol']:
                score += 0.25
                
        elif receptor == 'Sigma1':
            if features['has_basic_nitrogen'] and features['logP'] > 3:
                score += 0.3
        
        return min(1.0, score)
    
    def _get_safety_recommendation(self, risk: str, concerns: List) -> str:
        """Generate safety recommendation"""
        if risk == 'Critical':
            return "Critical safety concerns identified. Structural modification strongly recommended before proceeding."
        elif risk == 'High':
            return "High safety risk. Consider assay confirmation and structural optimization."
        elif risk == 'Moderate':
            return "Moderate safety concerns. Monitor in preclinical studies."
        else:
            return "No major safety concerns identified based on receptor profiling."
    
    def _calculate_therapeutic_scores(self, receptor_scores: Dict) -> Dict:
        """Calculate therapeutic potential scores by area"""
        area_scores = {}
        
        for family, data in receptor_scores.items():
            areas = data.get('therapeutic_areas', [])
            for area in areas:
                if area not in area_scores:
                    area_scores[area] = []
                weighted_score = data['score'] * THERAPEUTIC_WEIGHTS.get(area, 1.0)
                area_scores[area].append(weighted_score)
        
        therapeutic_potential = {}
        for area, scores in area_scores.items():
            therapeutic_potential[area] = {
                'score': round(max(scores), 3),
                'mean_score': round(np.mean(scores), 3),
                'target_count': len(scores),
                'priority': 'High' if max(scores) > 0.7 else 'Medium' if max(scores) > 0.5 else 'Low'
            }
        
        sorted_areas = sorted(
            therapeutic_potential.items(),
            key=lambda x: x[1]['score'],
            reverse=True
        )
        
        return {
            'by_area': therapeutic_potential,
            'top_areas': [area for area, _ in sorted_areas[:3]],
            'overall_potential': 'High' if any(d['score'] > 0.7 for d in therapeutic_potential.values()) else 'Moderate'
        }
    
    def _assess_lead_potential(
        self, 
        receptor_scores: Dict, 
        safety_profile: Dict, 
        features: Dict
    ) -> Dict:
        """Assess overall lead compound potential"""
        max_score = max(d['score'] for d in receptor_scores.values()) if receptor_scores else 0
        
        safety_penalty = 0
        if safety_profile.get('overall_risk') == 'Critical':
            safety_penalty = 0.3
        elif safety_profile.get('overall_risk') == 'High':
            safety_penalty = 0.15
        elif safety_profile.get('overall_risk') == 'Moderate':
            safety_penalty = 0.05
        
        lipinski_score = 1.0
        if features['molecular_weight'] > 500:
            lipinski_score -= 0.25
        if features['logP'] > 5:
            lipinski_score -= 0.25
        if features['hbd'] > 5:
            lipinski_score -= 0.25
        if features['hba'] > 10:
            lipinski_score -= 0.25
        
        lead_score = (max_score * 0.5 + lipinski_score * 0.3) * (1 - safety_penalty)
        
        return {
            'lead_score': round(lead_score, 3),
            'classification': 'Excellent' if lead_score > 0.75 else 'Good' if lead_score > 0.55 else 'Moderate' if lead_score > 0.35 else 'Poor',
            'activity_score': round(max_score, 3),
            'safety_penalty': round(safety_penalty, 3),
            'drug_likeness_score': round(lipinski_score, 3),
            'recommendation': self._get_lead_recommendation(lead_score, safety_penalty)
        }
    
    def _get_lead_recommendation(self, lead_score: float, safety_penalty: float) -> str:
        """Generate lead optimization recommendation"""
        if lead_score > 0.75 and safety_penalty < 0.1:
            return "Excellent lead candidate. Proceed to detailed profiling."
        elif lead_score > 0.55:
            if safety_penalty > 0.1:
                return "Good activity but safety concerns. Consider structural modifications."
            return "Promising lead. Further optimization recommended."
        elif lead_score > 0.35:
            return "Moderate potential. Significant optimization needed."
        else:
            return "Limited lead potential. Consider alternative scaffolds."
    
    def _get_top_targets(self, receptor_scores: Dict, n: int = 5) -> List[Dict]:
        """Get top N receptor targets by score"""
        sorted_targets = sorted(
            receptor_scores.items(),
            key=lambda x: x[1]['score'],
            reverse=True
        )[:n]
        
        return [
            {
                'family': family,
                'score': data['score'],
                'predicted_ki': data['predicted_ki'],
                'affinity_class': data['affinity_class'],
                'confidence': data['confidence']
            }
            for family, data in sorted_targets
        ]
    
    def _generate_recommendations(
        self, 
        receptor_scores: Dict, 
        safety_profile: Dict, 
        therapeutic_scores: Dict
    ) -> List[str]:
        """Generate actionable recommendations"""
        recommendations = []
        
        top_score = max(d['score'] for d in receptor_scores.values()) if receptor_scores else 0
        
        if top_score > 0.7:
            recommendations.append("Strong receptor affinity detected. Consider in vitro validation.")
        elif top_score > 0.5:
            recommendations.append("Moderate affinity. Structure optimization may improve binding.")
        
        if safety_profile.get('overall_risk') in ['Critical', 'High']:
            recommendations.append("Address safety liabilities before advancement.")
        
        if therapeutic_scores.get('overall_potential') == 'High':
            top_areas = therapeutic_scores.get('top_areas', [])
            if top_areas:
                recommendations.append(f"Strongest therapeutic potential in: {', '.join(top_areas[:2])}")
        
        return recommendations


def get_enhanced_docking_scores(smiles: str, target_families: List[str] = None) -> Dict:
    """Convenience function to get enhanced docking scores"""
    scorer = EnhancedDockingScorer()
    return scorer.calculate_enhanced_docking_score(smiles, target_families)
