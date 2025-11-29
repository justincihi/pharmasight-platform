#!/usr/bin/env python3
"""
Virtual High-Throughput Screening (vHTS) Pipeline
Screens compounds against all 82 receptor targets automatically
"""

import json
import numpy as np
from typing import Dict, List, Tuple
from datetime import datetime
import concurrent.futures
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs
import sqlite3
import hashlib

class VirtualScreeningPipeline:
    """Automated virtual screening against multiple receptor targets"""
    
    def __init__(self):
        self.load_receptors()
        self.screening_results = []
        self.polypharmacology_profiles = {}
        
    def load_receptors(self):
        """Load all receptor targets from database"""
        from src.receptor_database import get_all_receptor_targets
        from src.receptor_subtypes import RECEPTOR_SUBTYPES
        
        self.all_receptors = get_all_receptor_targets()
        self.receptor_families = self._classify_receptors()
        
    def _classify_receptors(self) -> Dict:
        """Classify receptors by family and therapeutic area"""
        families = {}
        for name, data in self.all_receptors.items():
            family = data.get('family', 'Unknown')
            if family not in families:
                families[family] = []
            families[family].append(name)
        return families
    
    def screen_compound(self, smiles: str, compound_name: str = None) -> Dict:
        """Screen a single compound against all receptors"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        compound_id = compound_name or f"CPD_{hashlib.md5(smiles.encode()).hexdigest()[:8]}"
        
        # Calculate molecular fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        
        screening_result = {
            "compound_id": compound_id,
            "smiles": smiles,
            "timestamp": datetime.now().isoformat(),
            "receptor_hits": [],
            "selectivity_profile": {},
            "polypharmacology_score": 0,
            "therapeutic_areas": set(),
            "off_targets": [],
            "safety_flags": []
        }
        
        # Screen against all receptors
        hit_count = 0
        high_affinity_hits = []
        
        for receptor_name, receptor_data in self.all_receptors.items():
            # Simulate binding affinity (in production, use real docking)
            binding_score = self._calculate_binding_score(mol, receptor_data)
            
            if binding_score > 0.5:  # Threshold for hit
                hit_count += 1
                hit_data = {
                    "receptor": receptor_name,
                    "family": receptor_data.get('family', 'Unknown'),
                    "binding_score": round(binding_score, 3),
                    "predicted_ki": self._score_to_ki(binding_score),
                    "activity_type": self._predict_activity_type(binding_score)
                }
                
                screening_result["receptor_hits"].append(hit_data)
                
                # Track therapeutic areas
                if 'therapeutic_relevance' in receptor_data:
                    screening_result["therapeutic_areas"].add(
                        receptor_data['therapeutic_relevance']
                    )
                
                if binding_score > 0.8:
                    high_affinity_hits.append(receptor_name)
        
        # Calculate selectivity and polypharmacology
        screening_result["selectivity_profile"] = self._calculate_selectivity(
            screening_result["receptor_hits"]
        )
        screening_result["polypharmacology_score"] = self._calculate_polypharmacology(
            screening_result["receptor_hits"]
        )
        
        # Identify off-targets and safety issues
        screening_result["off_targets"] = self._identify_off_targets(
            screening_result["receptor_hits"]
        )
        screening_result["safety_flags"] = self._assess_safety(
            screening_result["receptor_hits"]
        )
        
        # Convert set to list for JSON serialization
        screening_result["therapeutic_areas"] = list(screening_result["therapeutic_areas"])
        
        return screening_result
    
    def _calculate_binding_score(self, mol, receptor_data: Dict) -> float:
        """Calculate binding score using advanced pharmacophore matching for all receptor families"""
        import random
        from rdkit.Chem import rdMolDescriptors, Fragments
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        aliphatic_rings = rdMolDescriptors.CalcNumAliphaticRings(mol)
        total_rings = rdMolDescriptors.CalcNumRings(mol)
        num_heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
        fraction_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        
        has_basic_nitrogen = self._has_basic_nitrogen(mol)
        has_phenol = Fragments.fr_phenol(mol) > 0
        has_ether = Fragments.fr_ether(mol) > 0
        has_amine = Fragments.fr_NH2(mol) > 0 or Fragments.fr_NH1(mol) > 0 or Fragments.fr_NH0(mol) > 0
        has_amide = Fragments.fr_amide(mol) > 0
        has_halogen = Fragments.fr_halogen(mol) > 0
        has_sulfur = Fragments.fr_sulfide(mol) > 0 or Fragments.fr_sulfonamd(mol) > 0
        has_benzene = Fragments.fr_benzene(mol) > 0
        has_piperidine = Fragments.fr_piperdine(mol) > 0
        has_morpholine = Fragments.fr_morpholine(mol) > 0
        has_furan = Fragments.fr_furan(mol) > 0
        has_thiophene = Fragments.fr_thiophene(mol) > 0
        has_pyridine = Fragments.fr_pyridine(mol) > 0
        
        indole_pattern = Chem.MolFromSmarts('c1ccc2[nH]ccc2c1')
        has_indole = mol.HasSubstructMatch(indole_pattern) if indole_pattern else False
        
        base_score = 0.25
        family = receptor_data.get('family', '')
        receptor_name = receptor_data.get('name', '')
        receptor_type = receptor_data.get('type', 'GPCR')
        
        pharmacophore_score = self._get_pharmacophore_score(
            family, mol, mw, logp, hbd, hba, tpsa, aromatic_rings, aliphatic_rings,
            total_rings, rotatable_bonds, has_basic_nitrogen, has_phenol, has_ether,
            has_amine, has_amide, has_halogen, has_indole, has_benzene, has_piperidine,
            has_morpholine, has_furan, has_thiophene, has_pyridine, fraction_sp3,
            num_heteroatoms, has_sulfur, receptor_type
        )
        
        base_score += pharmacophore_score
        
        noise = random.gauss(0, 0.06)
        base_score += noise
        
        return max(0.0, min(1.0, base_score))
    
    def _has_basic_nitrogen(self, mol) -> bool:
        """Check if molecule contains a basic nitrogen (protonatable at physiological pH)"""
        from rdkit.Chem import AllChem
        pattern_tertiary = Chem.MolFromSmarts('[N;X3;!$(N-C=O);!$(N-S=O)]')
        pattern_secondary = Chem.MolFromSmarts('[NH1;!$(N-C=O);!$(N-S=O)]')
        pattern_primary = Chem.MolFromSmarts('[NH2;!$(N-C=O)]')
        
        has_basic = False
        if pattern_tertiary and mol.HasSubstructMatch(pattern_tertiary):
            has_basic = True
        if pattern_secondary and mol.HasSubstructMatch(pattern_secondary):
            has_basic = True
        if pattern_primary and mol.HasSubstructMatch(pattern_primary):
            has_basic = True
        return has_basic
    
    def _get_pharmacophore_score(self, family: str, mol, mw, logp, hbd, hba, tpsa,
                                  aromatic_rings, aliphatic_rings, total_rings,
                                  rotatable_bonds, has_basic_nitrogen, has_phenol,
                                  has_ether, has_amine, has_amide, has_halogen,
                                  has_indole, has_benzene, has_piperidine,
                                  has_morpholine, has_furan, has_thiophene,
                                  has_pyridine, fraction_sp3, num_heteroatoms,
                                  has_sulfur, receptor_type) -> float:
        """Calculate pharmacophore-based score for each receptor family"""
        
        score = 0.0
        
        if family == 'Serotonin':
            if has_indole:
                score += 0.20
            elif aromatic_rings >= 2 and has_basic_nitrogen:
                score += 0.12
            if has_basic_nitrogen and 2 <= logp <= 4:
                score += 0.10
            if 180 < mw < 380:
                score += 0.08
            if 20 < tpsa < 80:
                score += 0.05
            if has_phenol and aromatic_rings >= 1:
                score += 0.05
                
        elif family == 'Dopamine':
            if has_phenol and has_basic_nitrogen:
                score += 0.22
            elif aromatic_rings >= 1 and has_basic_nitrogen:
                score += 0.12
            if 1.5 <= logp <= 4.5:
                score += 0.10
            if 150 < mw < 450:
                score += 0.08
            if has_piperidine or has_morpholine:
                score += 0.08
            if 30 < tpsa < 90:
                score += 0.05
                
        elif family == 'Adrenergic':
            if has_phenol and has_basic_nitrogen and hbd >= 2:
                score += 0.22
            elif aromatic_rings >= 1 and has_basic_nitrogen:
                score += 0.10
            if 0.5 <= logp <= 3.5:
                score += 0.10
            if 150 < mw < 400:
                score += 0.08
            if 40 < tpsa < 100:
                score += 0.06
            if has_ether:
                score += 0.04
                
        elif family == 'Opioid':
            if total_rings >= 3 and has_basic_nitrogen:
                score += 0.22
            elif has_piperidine and aromatic_rings >= 1:
                score += 0.18
            if 250 < mw < 550:
                score += 0.10
            if 1.5 <= logp <= 5.0:
                score += 0.08
            if hbd >= 1 and hba >= 2:
                score += 0.06
            if aliphatic_rings >= 2:
                score += 0.06
                
        elif family == 'GABA':
            if receptor_type == 'Ion channel':
                if has_amide or (hbd >= 2 and hba >= 3):
                    score += 0.22
                if 150 < mw < 400:
                    score += 0.10
                if -1 <= logp <= 3:
                    score += 0.10
                if 50 < tpsa < 120:
                    score += 0.08
                if has_halogen and aromatic_rings >= 1:
                    score += 0.08
            else:
                if has_basic_nitrogen and aromatic_rings >= 1:
                    score += 0.18
                if 200 < mw < 500:
                    score += 0.08
                    
        elif family == 'Glutamate':
            if hbd >= 2 and hba >= 4:
                score += 0.20
            if 100 < mw < 350:
                score += 0.12
            if -2 <= logp <= 2:
                score += 0.12
            if tpsa > 80:
                score += 0.08
            if has_amine and not has_benzene:
                score += 0.08
                
        elif family == 'Cannabinoid':
            if aromatic_rings >= 1 and aliphatic_rings >= 1:
                score += 0.18
            if 3 <= logp <= 7:
                score += 0.15
            if 280 < mw < 500:
                score += 0.10
            if has_phenol:
                score += 0.08
            if rotatable_bonds >= 4:
                score += 0.06
            if 40 < tpsa < 80:
                score += 0.05
                
        elif family == 'Muscarinic':
            if has_basic_nitrogen and (has_ether or total_rings >= 2):
                score += 0.22
            if 200 < mw < 450:
                score += 0.10
            if 1 <= logp <= 4:
                score += 0.10
            if 30 < tpsa < 80:
                score += 0.06
            if has_piperidine or aliphatic_rings >= 1:
                score += 0.06
                
        elif family == 'Nicotinic':
            if has_pyridine and has_basic_nitrogen:
                score += 0.25
            elif aromatic_rings >= 1 and has_basic_nitrogen:
                score += 0.12
            if 120 < mw < 350:
                score += 0.10
            if 0 <= logp <= 3:
                score += 0.10
            if 20 < tpsa < 70:
                score += 0.05
                
        elif family == 'Histamine':
            if has_basic_nitrogen and aromatic_rings >= 1:
                score += 0.18
            if 200 < mw < 400:
                score += 0.12
            if 1 <= logp <= 4:
                score += 0.10
            if 30 < tpsa < 80:
                score += 0.08
            if has_piperidine or has_morpholine:
                score += 0.06
            if has_halogen:
                score += 0.04
                
        elif family == 'Adenosine':
            if mw < 180 or mw > 500:
                return 0.0
            if total_rings >= 2 and num_heteroatoms >= 5:
                score += 0.25
            elif total_rings >= 2 and num_heteroatoms >= 4:
                score += 0.15
            if has_furan and num_heteroatoms >= 4:
                score += 0.12
            elif has_ether and num_heteroatoms >= 3:
                score += 0.06
            if 220 < mw < 420:
                score += 0.08
            if -1 <= logp <= 2.5:
                score += 0.08
            if hba >= 5:
                score += 0.06
                
        elif family == 'TRP':
            if mw < 220 or mw > 550:
                return 0.0
            if aromatic_rings >= 1 and 2.5 <= logp <= 5.5:
                score += 0.18
            if 280 < mw < 480:
                score += 0.10
            if hbd >= 1 and hba >= 2 and has_amide:
                score += 0.12
            elif hbd >= 1 and hba >= 2:
                score += 0.06
            if 50 < tpsa < 100:
                score += 0.05
                
        elif family == 'Purinergic':
            if mw < 180 or mw > 550:
                return 0.0
            if num_heteroatoms >= 5 and total_rings >= 2:
                score += 0.25
            elif num_heteroatoms >= 4 and total_rings >= 2:
                score += 0.15
            if 220 < mw < 480:
                score += 0.08
            if -1 <= logp <= 2.5:
                score += 0.10
            if tpsa > 90:
                score += 0.08
            if hbd >= 2 and hba >= 4:
                score += 0.06
                
        elif family == 'Orexin':
            if mw < 320 or mw > 650:
                return 0.0
            if has_amide and aromatic_rings >= 2 and total_rings >= 3:
                score += 0.25
            if 380 < mw < 580:
                score += 0.12
            if 2.5 <= logp <= 5:
                score += 0.08
            if 70 < tpsa < 130:
                score += 0.06
            if has_sulfur and has_amide:
                score += 0.08
                
        elif family == 'Melatonin':
            if not has_indole:
                if mw < 180 or mw > 380:
                    return 0.0
            if has_indole and has_amide:
                score += 0.35
            elif has_indole:
                score += 0.22
            elif aromatic_rings >= 1 and has_amide and 200 < mw < 350:
                score += 0.12
            if 200 < mw < 320:
                score += 0.08
            if 1 <= logp <= 3:
                score += 0.06
            if 45 < tpsa < 75:
                score += 0.04
                
        elif family == 'Chemokine':
            if mw < 320 or mw > 700:
                return 0.0
            if tpsa < 50 or tpsa > 160:
                return 0.0
            if has_basic_nitrogen and aromatic_rings >= 2 and total_rings >= 3:
                score += 0.22
            if 380 < mw < 600:
                score += 0.10
            if 2.5 <= logp <= 5.5:
                score += 0.08
            if 70 < tpsa < 140:
                score += 0.06
            if has_amide and has_piperidine:
                score += 0.08
                
        else:
            if aromatic_rings >= 1 and has_basic_nitrogen:
                score += 0.12
            if 200 < mw < 500:
                score += 0.08
            if 0 <= logp <= 5:
                score += 0.08
            if 30 < tpsa < 100:
                score += 0.05
        
        return score
    
    def _score_to_ki(self, score: float) -> str:
        """Convert binding score to predicted Ki"""
        if score > 0.9:
            return "<10 nM"
        elif score > 0.8:
            return "10-100 nM"
        elif score > 0.7:
            return "100-500 nM"
        elif score > 0.6:
            return "0.5-5 µM"
        elif score > 0.5:
            return "5-50 µM"
        else:
            return ">50 µM"
    
    def _predict_activity_type(self, score: float) -> str:
        """Predict if compound is agonist, antagonist, or modulator"""
        if score > 0.85:
            return "Strong binder (activity unknown)"
        elif score > 0.7:
            return "Moderate binder"
        elif score > 0.5:
            return "Weak binder"
        else:
            return "Non-binder"
    
    def _calculate_selectivity(self, hits: List[Dict]) -> Dict:
        """Calculate selectivity profile"""
        if not hits:
            return {"status": "No hits", "selectivity_index": 0}
        
        # Group hits by family
        family_hits = {}
        for hit in hits:
            family = hit['family']
            if family not in family_hits:
                family_hits[family] = []
            family_hits[family].append(hit['binding_score'])
        
        # Calculate selectivity index
        if len(family_hits) == 1:
            selectivity = "Highly selective"
            index = 10.0
        elif len(family_hits) == 2:
            selectivity = "Moderately selective"
            index = 5.0
        else:
            selectivity = "Promiscuous"
            index = 1.0
        
        return {
            "status": selectivity,
            "selectivity_index": index,
            "family_distribution": {k: len(v) for k, v in family_hits.items()}
        }
    
    def _calculate_polypharmacology(self, hits: List[Dict]) -> float:
        """Calculate polypharmacology score"""
        if not hits:
            return 0.0
        
        # Score based on number of targets and diversity
        num_hits = len(hits)
        families = set(hit['family'] for hit in hits)
        
        # Polypharmacology score (0-100)
        score = min(100, (num_hits * 5) + (len(families) * 10))
        
        return score
    
    def _identify_off_targets(self, hits: List[Dict]) -> List[Dict]:
        """Identify potential off-target effects"""
        off_targets = []
        
        critical_off_targets = {
            "5-HT2B": "Cardiac valvulopathy risk",
            "hERG": "QT prolongation risk",
            "M1": "Cognitive side effects",
            "H1": "Sedation",
            "Alpha1A": "Orthostatic hypotension",
            "D2": "Extrapyramidal symptoms"
        }
        
        for hit in hits:
            receptor = hit['receptor']
            if receptor in critical_off_targets and hit['binding_score'] > 0.6:
                off_targets.append({
                    "receptor": receptor,
                    "concern": critical_off_targets[receptor],
                    "binding_score": hit['binding_score'],
                    "risk_level": "High" if hit['binding_score'] > 0.8 else "Moderate"
                })
        
        return off_targets
    
    def _assess_safety(self, hits: List[Dict]) -> List[str]:
        """Assess safety based on receptor profile"""
        safety_flags = []
        
        # Check for problematic combinations
        receptor_names = [hit['receptor'] for hit in hits if hit['binding_score'] > 0.7]
        
        # Serotonin syndrome risk
        serotonin_count = sum(1 for r in receptor_names if '5-HT' in r)
        if serotonin_count > 3:
            safety_flags.append("Serotonin syndrome risk (multiple 5-HT receptors)")
        
        # CNS depression risk
        if any(r in receptor_names for r in ['GABA-A', 'MOR', 'H1']):
            safety_flags.append("CNS depression risk")
        
        # Cardiovascular risk
        if '5-HT2B' in receptor_names:
            safety_flags.append("Cardiac valvulopathy risk (5-HT2B agonism)")
        
        return safety_flags
    
    def batch_screen(self, compound_list: List[Tuple[str, str]], 
                     parallel: bool = True) -> List[Dict]:
        """Screen multiple compounds in parallel"""
        results = []
        
        if parallel:
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                futures = []
                for smiles, name in compound_list:
                    future = executor.submit(self.screen_compound, smiles, name)
                    futures.append(future)
                
                for future in concurrent.futures.as_completed(futures):
                    results.append(future.result())
        else:
            for smiles, name in compound_list:
                results.append(self.screen_compound(smiles, name))
        
        return results
    
    def find_selective_ligands(self, target_receptor: str, 
                               avoid_receptors: List[str] = None) -> List[Dict]:
        """Find compounds selective for a specific receptor"""
        selective_compounds = []
        
        for result in self.screening_results:
            hits = result['receptor_hits']
            target_hit = None
            avoid_hits = []
            
            for hit in hits:
                if hit['receptor'] == target_receptor:
                    target_hit = hit
                elif avoid_receptors and hit['receptor'] in avoid_receptors:
                    avoid_hits.append(hit)
            
            # Check selectivity criteria
            if target_hit and target_hit['binding_score'] > 0.7:
                if not avoid_hits or all(h['binding_score'] < 0.4 for h in avoid_hits):
                    selective_compounds.append({
                        "compound": result['compound_id'],
                        "smiles": result['smiles'],
                        "target_affinity": target_hit['binding_score'],
                        "selectivity": "High" if not avoid_hits else "Moderate"
                    })
        
        return selective_compounds
    
    def identify_drug_repurposing(self) -> List[Dict]:
        """Identify repurposing opportunities from screening results"""
        repurposing_candidates = []
        
        for result in self.screening_results:
            if result['polypharmacology_score'] > 30:
                # Look for unexpected therapeutic areas
                areas = result['therapeutic_areas']
                if len(areas) > 1:
                    repurposing_candidates.append({
                        "compound": result['compound_id'],
                        "original_target": result['receptor_hits'][0]['receptor'] if result['receptor_hits'] else None,
                        "new_therapeutic_areas": list(areas),
                        "confidence": result['polypharmacology_score']
                    })
        
        return repurposing_candidates
    
    def generate_screening_report(self, results: List[Dict] = None) -> Dict:
        """Generate comprehensive screening report"""
        if results:
            self.screening_results = results
        
        report = {
            "screening_date": datetime.now().isoformat(),
            "total_compounds": len(self.screening_results),
            "summary": {
                "compounds_with_hits": 0,
                "selective_compounds": 0,
                "promiscuous_compounds": 0,
                "safety_concerns": 0
            },
            "top_hits": [],
            "selective_ligands": {},
            "polypharmacology_profiles": [],
            "repurposing_opportunities": []
        }
        
        # Analyze results
        for result in self.screening_results:
            if result.get('receptor_hits'):
                report['summary']['compounds_with_hits'] += 1
                
                # Check selectivity
                selectivity = result.get('selectivity_profile', {})
                if selectivity.get('status') == 'Highly selective':
                    report['summary']['selective_compounds'] += 1
                elif selectivity.get('status') == 'Promiscuous':
                    report['summary']['promiscuous_compounds'] += 1
                
                # Check safety
                if result.get('safety_flags'):
                    report['summary']['safety_concerns'] += 1
                
                # Top hits (highest binding scores)
                for hit in result['receptor_hits']:
                    if hit['binding_score'] > 0.85:
                        report['top_hits'].append({
                            "compound": result['compound_id'],
                            "receptor": hit['receptor'],
                            "score": hit['binding_score']
                        })
        
        # Sort top hits by score
        report['top_hits'] = sorted(
            report['top_hits'], 
            key=lambda x: x['score'], 
            reverse=True
        )[:20]
        
        # Add repurposing opportunities
        report['repurposing_opportunities'] = self.identify_drug_repurposing()
        
        return report