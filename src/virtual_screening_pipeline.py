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
        """Calculate binding score using pharmacophore matching and ML"""
        # Simplified scoring based on molecular properties
        # In production, this would use real docking or ML models
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Receptor-specific scoring
        base_score = 0.5
        
        # Check known ligands similarity
        if 'agonists' in receptor_data:
            for known_ligand in receptor_data.get('agonists', [])[:3]:
                # Simulate similarity check
                base_score += 0.1
        
        # Apply receptor family-specific rules
        family = receptor_data.get('family', '')
        
        if family == 'GABA':
            if hbd > 1 and hba > 2:
                base_score += 0.2
        elif family == 'Serotonin':
            if 150 < mw < 400 and logp > 1:
                base_score += 0.15
        elif family == 'Dopamine':
            if hba > 1 and logp > 2:
                base_score += 0.15
        elif family == 'Opioid':
            if mw > 250 and hbd > 0:
                base_score += 0.2
        
        # Add some randomness for simulation
        import random
        base_score += random.uniform(-0.2, 0.3)
        
        return max(0, min(1, base_score))
    
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