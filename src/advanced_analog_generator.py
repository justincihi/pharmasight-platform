#!/usr/bin/env python3
"""
Advanced Analog Generator for PharmaSight™
Structure-based optimization with scaffold hopping, R-group enumeration,
matched molecular pair analysis, and multi-objective optimization.
Integrates with Enhanced Docking Scorer for target-aware analog design.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, DataStructs
from rdkit.Chem import BRICS, rdFMCS, RWMol, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from typing import List, Dict, Tuple, Optional, Set
import numpy as np
import hashlib
from datetime import datetime
from collections import defaultdict

try:
    from enhanced_docking_scorer import EnhancedDockingScorer
except ImportError:
    try:
        from src.enhanced_docking_scorer import EnhancedDockingScorer
    except ImportError:
        EnhancedDockingScorer = None

try:
    from comprehensive_viability import ComprehensiveViabilityAnalyzer
except ImportError:
    try:
        from src.comprehensive_viability import ComprehensiveViabilityAnalyzer
    except ImportError:
        ComprehensiveViabilityAnalyzer = None


class AdvancedAnalogGenerator:
    """
    Advanced analog generation with structure-based optimization.
    Integrates scaffold hopping, R-group enumeration, matched molecular pairs,
    and multi-objective scoring using enhanced docking profiles.
    """
    
    SCAFFOLD_BIOISOSTERES = {
        'benzene': ['pyridine', 'pyrimidine', 'thiophene', 'furan', 'pyrrole', 'cyclohexane'],
        'pyridine': ['benzene', 'pyrimidine', 'pyrazine', 'pyridazine', 'thiazole'],
        'phenyl': ['naphthyl', 'biphenyl', 'indolyl', 'benzofuranyl'],
        'indole': ['benzimidazole', 'benzofuran', 'benzothiophene', 'azaindole'],
        'piperidine': ['piperazine', 'morpholine', 'pyrrolidine', 'tropane'],
        'cyclohexane': ['benzene', 'tetrahydropyran', 'cyclohexene'],
    }
    
    SCAFFOLD_SMARTS = {
        'benzene': 'c1ccccc1',
        'pyridine': 'c1ccncc1',
        'pyrimidine': 'c1ncncc1',
        'thiophene': 'c1ccsc1',
        'furan': 'c1ccoc1',
        'pyrrole': 'c1cc[nH]c1',
        'cyclohexane': 'C1CCCCC1',
        'piperidine': 'C1CCNCC1',
        'piperazine': 'C1CNCCN1',
        'morpholine': 'C1COCCN1',
        'pyrrolidine': 'C1CCNC1',
        'indole': 'c1ccc2[nH]ccc2c1',
        'benzimidazole': 'c1ccc2nc[nH]c2c1',
        'pyrazine': 'c1cnccn1',
    }
    
    R_GROUP_LIBRARY = {
        'small_hydrophobic': ['C', 'CC', 'C(C)C', 'C1CC1', 'F', 'Cl'],
        'polar': ['O', 'N', 'CO', 'CN', 'C(=O)O', 'C(=O)N'],
        'electron_withdrawing': ['F', 'Cl', 'C(F)(F)F', 'C#N', '[N+](=O)[O-]', 'S(=O)(=O)C'],
        'electron_donating': ['OC', 'NC', 'N(C)C', 'SC'],
        'fluorinated': ['F', 'C(F)(F)F', 'C(F)F', 'OC(F)(F)F'],
        'heterocyclic': ['c1cccnc1', 'c1ccoc1', 'c1cc[nH]c1', 'C1CCNCC1'],
        'linkers': ['O', 'N', 'S', 'C=O', 'C(=O)N'],
    }
    
    MATCHED_PAIR_TRANSFORMS = [
        {'from': '[CH3:1]', 'to': '[CH2F:1]', 'name': 'Methyl to fluoromethyl', 'effect': 'metabolic_stability'},
        {'from': '[CH3:1]', 'to': '[CF3:1]', 'name': 'Methyl to trifluoromethyl', 'effect': 'lipophilicity'},
        {'from': '[cH:1]', 'to': '[c:1]F', 'name': 'Aromatic H to F', 'effect': 'metabolic_stability'},
        {'from': '[cH:1]', 'to': '[c:1]Cl', 'name': 'Aromatic H to Cl', 'effect': 'potency'},
        {'from': '[OH:1]', 'to': '[F:1]', 'name': 'Hydroxyl to fluorine', 'effect': 'stability'},
        {'from': '[NH2:1]', 'to': '[NHCH3:1]', 'name': 'Primary to secondary amine', 'effect': 'basicity'},
        {'from': 'C(=O)O', 'to': 'C(=O)N', 'name': 'Acid to amide', 'effect': 'metabolic_stability'},
        {'from': 'c1ccccc1', 'to': 'c1ccncc1', 'name': 'Phenyl to pyridyl', 'effect': 'solubility'},
        {'from': '[CH2:1][CH2:1]', 'to': '[CH2:1]O[CH2:1]', 'name': 'Ethylene to ether', 'effect': 'flexibility'},
        {'from': 'C1CCCCC1', 'to': 'c1ccccc1', 'name': 'Cyclohexyl to phenyl', 'effect': 'planarity'},
    ]
    
    def __init__(self):
        self.docking_scorer = EnhancedDockingScorer() if EnhancedDockingScorer else None
        self.viability_analyzer = ComprehensiveViabilityAnalyzer() if ComprehensiveViabilityAnalyzer else None
        self.generated_cache: Dict[str, Dict] = {}
        self._init_pains_filter()
    
    def _init_pains_filter(self):
        """Initialize PAINS (Pan-Assay Interference) filter"""
        try:
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            self.pains_catalog = FilterCatalog(params)
        except:
            self.pains_catalog = None
    
    def validate_molecule(self, smiles: str) -> Tuple[bool, Optional[Chem.Mol], str]:
        """Validate SMILES and return molecule"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, None, "Invalid SMILES structure"
            Chem.SanitizeMol(mol)
            canonical = Chem.MolToSmiles(mol)
            return True, mol, canonical
        except Exception as e:
            return False, None, str(e)
    
    def extract_murcko_scaffold(self, mol: Chem.Mol) -> Tuple[str, str]:
        """Extract Murcko scaffold (framework and generic)"""
        try:
            framework = MurckoScaffold.GetScaffoldForMol(mol)
            generic = MurckoScaffold.MakeScaffoldGeneric(framework)
            return Chem.MolToSmiles(framework), Chem.MolToSmiles(generic)
        except:
            return "", ""
    
    def decompose_brics(self, mol: Chem.Mol) -> List[str]:
        """Decompose molecule using BRICS fragmentation"""
        try:
            fragments = list(BRICS.BRICSDecompose(mol))
            return fragments
        except:
            return []
    
    def calculate_fingerprint(self, mol: Chem.Mol, fp_type: str = 'morgan') -> Optional[object]:
        """Calculate molecular fingerprint"""
        try:
            if fp_type == 'morgan':
                return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            elif fp_type == 'rdkit':
                return Chem.RDKFingerprint(mol)
            elif fp_type == 'maccs':
                return AllChem.GetMACCSKeysFingerprint(mol)
        except:
            return None
    
    def calculate_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> float:
        """Calculate Tanimoto similarity"""
        fp1 = self.calculate_fingerprint(mol1)
        fp2 = self.calculate_fingerprint(mol2)
        if fp1 and fp2:
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        return 0.0
    
    def check_pains(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """Check for PAINS alerts"""
        if not self.pains_catalog:
            return True, []
        
        entries = self.pains_catalog.GetMatches(mol)
        if entries:
            alerts = [entry.GetDescription() for entry in entries]
            return False, alerts
        return True, []
    
    def calculate_properties(self, mol: Chem.Mol) -> Dict:
        """Calculate comprehensive molecular properties"""
        try:
            props = {
                'mw': round(Descriptors.MolWt(mol), 2),
                'logp': round(Descriptors.MolLogP(mol), 2),
                'hbd': rdMolDescriptors.CalcNumHBD(mol),
                'hba': rdMolDescriptors.CalcNumHBA(mol),
                'tpsa': round(rdMolDescriptors.CalcTPSA(mol), 2),
                'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'rings': rdMolDescriptors.CalcNumRings(mol),
                'fraction_sp3': round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
                'num_stereocenters': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            }
            
            props['lipinski_violations'] = sum([
                props['mw'] > 500,
                props['logp'] > 5,
                props['hbd'] > 5,
                props['hba'] > 10
            ])
            
            props['veber_compliant'] = props['rotatable_bonds'] <= 10 and props['tpsa'] <= 140
            
            return props
        except Exception as e:
            return {'error': str(e)}
    
    def scaffold_hop(self, mol: Chem.Mol, max_hops: int = 5) -> List[Dict]:
        """
        Perform scaffold hopping by replacing core scaffolds with bioisosteres.
        Returns list of hopped molecules with metadata.
        """
        results = []
        smiles = Chem.MolToSmiles(mol)
        
        for scaffold_name, scaffold_smarts in self.SCAFFOLD_SMARTS.items():
            pattern = Chem.MolFromSmarts(scaffold_smarts)
            if not pattern:
                continue
            
            if mol.HasSubstructMatch(pattern):
                bioisosteres = self.SCAFFOLD_BIOISOSTERES.get(scaffold_name, [])
                
                for replacement_name in bioisosteres[:max_hops]:
                    replacement_smarts = self.SCAFFOLD_SMARTS.get(replacement_name)
                    if not replacement_smarts:
                        continue
                    
                    try:
                        rxn_smarts = f"[{scaffold_smarts}:1]>>[{replacement_smarts}:1]"
                        
                        new_smiles = smiles.replace(
                            scaffold_smarts.replace('c', 'C').replace('n', 'N'),
                            replacement_smarts.replace('c', 'C').replace('n', 'N')
                        )
                        
                        if new_smiles == smiles:
                            new_smiles = smiles.replace(scaffold_smarts, replacement_smarts)
                        
                        new_mol = Chem.MolFromSmiles(new_smiles)
                        if new_mol:
                            Chem.SanitizeMol(new_mol)
                            canonical = Chem.MolToSmiles(new_mol)
                            
                            results.append({
                                'smiles': canonical,
                                'mol': new_mol,
                                'original_scaffold': scaffold_name,
                                'new_scaffold': replacement_name,
                                'hop_type': 'scaffold_replacement',
                                'similarity': self.calculate_similarity(mol, new_mol)
                            })
                    except:
                        continue
        
        seen = set()
        unique_results = []
        for r in results:
            if r['smiles'] not in seen:
                seen.add(r['smiles'])
                unique_results.append(r)
        
        return unique_results
    
    def enumerate_r_groups(self, mol: Chem.Mol, position_smarts: str = "[CH3]",
                           r_group_category: str = 'small_hydrophobic',
                           max_analogs: int = 10) -> List[Dict]:
        """
        Enumerate R-groups at specified positions.
        """
        results = []
        r_groups = self.R_GROUP_LIBRARY.get(r_group_category, self.R_GROUP_LIBRARY['small_hydrophobic'])
        
        pattern = Chem.MolFromSmarts(position_smarts)
        if not pattern:
            return results
        
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            return results
        
        for match in matches[:3]:
            atom_idx = match[0]
            
            for r_group in r_groups:
                try:
                    rw_mol = RWMol(mol)
                    r_mol = Chem.MolFromSmiles(r_group)
                    if not r_mol:
                        continue
                    
                    combined = Chem.CombineMols(rw_mol, r_mol)
                    rw_combined = RWMol(combined)
                    
                    num_orig_atoms = mol.GetNumAtoms()
                    rw_combined.AddBond(atom_idx, num_orig_atoms, Chem.BondType.SINGLE)
                    
                    new_mol = rw_combined.GetMol()
                    Chem.SanitizeMol(new_mol)
                    canonical = Chem.MolToSmiles(new_mol)
                    
                    results.append({
                        'smiles': canonical,
                        'mol': new_mol,
                        'r_group': r_group,
                        'position': atom_idx,
                        'category': r_group_category,
                        'similarity': self.calculate_similarity(mol, new_mol)
                    })
                    
                    if len(results) >= max_analogs:
                        break
                except:
                    continue
            
            if len(results) >= max_analogs:
                break
        
        return results
    
    def matched_molecular_pair_transform(self, mol: Chem.Mol, 
                                          target_effect: Optional[str] = None) -> List[Dict]:
        """
        Apply matched molecular pair transformations.
        Optionally filter by target effect (metabolic_stability, potency, solubility, etc.)
        """
        results = []
        smiles = Chem.MolToSmiles(mol)
        
        transforms = self.MATCHED_PAIR_TRANSFORMS
        if target_effect:
            transforms = [t for t in transforms if t['effect'] == target_effect]
        
        for transform in transforms:
            try:
                from_pattern = transform['from']
                to_pattern = transform['to']
                
                rxn_smarts = f"{from_pattern}>>{to_pattern}"
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                
                if rxn:
                    products = rxn.RunReactants((mol,))
                    
                    for product_set in products[:3]:
                        for product in product_set:
                            try:
                                Chem.SanitizeMol(product)
                                product_smiles = Chem.MolToSmiles(product)
                                
                                if product_smiles != smiles:
                                    results.append({
                                        'smiles': product_smiles,
                                        'mol': product,
                                        'transform_name': transform['name'],
                                        'expected_effect': transform['effect'],
                                        'from_pattern': from_pattern,
                                        'to_pattern': to_pattern,
                                        'similarity': self.calculate_similarity(mol, product)
                                    })
                            except:
                                continue
            except:
                continue
        
        seen = set()
        unique_results = []
        for r in results:
            if r['smiles'] not in seen:
                seen.add(r['smiles'])
                unique_results.append(r)
        
        return unique_results
    
    def multi_objective_score(self, mol: Chem.Mol, 
                               target_receptor: Optional[str] = None,
                               weights: Optional[Dict] = None) -> Dict:
        """
        Calculate multi-objective score balancing:
        - Activity (docking affinity)
        - Drug-likeness (QED-like)
        - Safety (PAINS, toxicity alerts)
        - Synthesizability
        """
        default_weights = {
            'activity': 0.35,
            'drug_likeness': 0.25,
            'safety': 0.25,
            'synthesizability': 0.15
        }
        weights = weights or default_weights
        
        scores = {}
        docking_result = None
        
        if self.docking_scorer:
            try:
                smiles = Chem.MolToSmiles(mol)
                docking_result = self.docking_scorer.calculate_enhanced_docking_score(smiles)
                lead_potential = docking_result.get('lead_potential', {})
                scores['activity'] = lead_potential.get('activity_score', lead_potential.get('lead_score', 0.5))
                top_targets = docking_result.get('top_targets', [])
                scores['top_receptor'] = top_targets[0]['family'] if top_targets else 'Unknown'
            except Exception as e:
                scores['activity'] = 0.5
                scores['top_receptor'] = 'Unknown'
        else:
            scores['activity'] = 0.5
            scores['top_receptor'] = 'Unknown'
        
        props = self.calculate_properties(mol)
        drug_likeness = 1.0
        
        if props.get('mw', 0) > 500:
            drug_likeness -= 0.15
        if props.get('logp', 0) > 5 or props.get('logp', 0) < -1:
            drug_likeness -= 0.15
        if props.get('hbd', 0) > 5:
            drug_likeness -= 0.15
        if props.get('hba', 0) > 10:
            drug_likeness -= 0.15
        if props.get('tpsa', 0) > 140:
            drug_likeness -= 0.1
        if props.get('rotatable_bonds', 0) > 10:
            drug_likeness -= 0.1
        if props.get('fraction_sp3', 0) < 0.25:
            drug_likeness -= 0.1
        
        scores['drug_likeness'] = max(0, drug_likeness)
        
        pains_clean, pains_alerts = self.check_pains(mol)
        safety_score = 1.0 if pains_clean else max(0, 1.0 - len(pains_alerts) * 0.2)
        
        if docking_result:
            try:
                safety_assessment = docking_result.get('safety_assessment', {})
                safety_risk = safety_assessment.get('overall_risk', 'Low')
                if safety_risk == 'Critical':
                    safety_score *= 0.5
                elif safety_risk == 'High':
                    safety_score *= 0.7
                elif safety_risk == 'Moderate':
                    safety_score *= 0.85
            except:
                pass
        
        scores['safety'] = safety_score
        scores['pains_alerts'] = pains_alerts
        
        sa_score = 5.0
        num_atoms = mol.GetNumHeavyAtoms()
        num_rings = props.get('rings', 0)
        num_stereo = props.get('num_stereocenters', 0)
        
        sa_score += num_atoms * 0.03
        sa_score += num_rings * 0.3
        sa_score += num_stereo * 0.8
        
        sa_normalized = max(0, 1.0 - (sa_score - 1) / 9)
        scores['synthesizability'] = sa_normalized
        
        combined = (
            weights['activity'] * scores['activity'] +
            weights['drug_likeness'] * scores['drug_likeness'] +
            weights['safety'] * scores['safety'] +
            weights['synthesizability'] * scores['synthesizability']
        )
        
        scores['combined'] = round(combined, 4)
        scores['weights_used'] = weights
        scores['properties'] = props
        
        return scores
    
    def generate_optimized_analogs(self, smiles: str,
                                    optimization_target: str = 'balanced',
                                    num_analogs: int = 20,
                                    min_similarity: float = 0.4,
                                    therapeutic_area: Optional[str] = None) -> Dict:
        """
        Generate optimized analogs using multiple strategies.
        
        Args:
            smiles: Input compound SMILES
            optimization_target: 'balanced', 'potency', 'safety', 'synthesizability'
            num_analogs: Maximum analogs to return
            min_similarity: Minimum Tanimoto similarity to parent
            therapeutic_area: Optional therapeutic focus (oncology, cns, cardiovascular, etc.)
        
        Returns:
            Dictionary with parent info, analogs, and optimization summary
        """
        valid, parent_mol, canonical = self.validate_molecule(smiles)
        if not valid:
            return {'error': f'Invalid SMILES: {canonical}'}
        
        parent_props = self.calculate_properties(parent_mol)
        parent_score = self.multi_objective_score(parent_mol)
        parent_framework, parent_generic = self.extract_murcko_scaffold(parent_mol)
        
        weight_profiles = {
            'balanced': {'activity': 0.35, 'drug_likeness': 0.25, 'safety': 0.25, 'synthesizability': 0.15},
            'potency': {'activity': 0.5, 'drug_likeness': 0.2, 'safety': 0.2, 'synthesizability': 0.1},
            'safety': {'activity': 0.2, 'drug_likeness': 0.25, 'safety': 0.45, 'synthesizability': 0.1},
            'synthesizability': {'activity': 0.25, 'drug_likeness': 0.2, 'safety': 0.2, 'synthesizability': 0.35},
        }
        weights = weight_profiles.get(optimization_target, weight_profiles['balanced'])
        
        all_analogs = []
        seen_smiles = {canonical}
        
        scaffold_hops = self.scaffold_hop(parent_mol, max_hops=5)
        for hop in scaffold_hops:
            if hop['smiles'] not in seen_smiles and hop['similarity'] >= min_similarity:
                seen_smiles.add(hop['smiles'])
                all_analogs.append({
                    'smiles': hop['smiles'],
                    'mol': hop['mol'],
                    'strategy': 'scaffold_hop',
                    'detail': f"{hop['original_scaffold']} → {hop['new_scaffold']}",
                    'similarity': hop['similarity']
                })
        
        for category in ['electron_withdrawing', 'polar', 'fluorinated']:
            r_groups = self.enumerate_r_groups(parent_mol, "[CH3]", category, max_analogs=5)
            for rg in r_groups:
                if rg['smiles'] not in seen_smiles and rg['similarity'] >= min_similarity:
                    seen_smiles.add(rg['smiles'])
                    all_analogs.append({
                        'smiles': rg['smiles'],
                        'mol': rg['mol'],
                        'strategy': 'r_group_enum',
                        'detail': f"{category}: {rg['r_group']}",
                        'similarity': rg['similarity']
                    })
        
        target_effects = ['metabolic_stability', 'potency', 'solubility'] if optimization_target == 'balanced' else [optimization_target]
        for effect in target_effects:
            mmp = self.matched_molecular_pair_transform(parent_mol, target_effect=effect if effect in ['metabolic_stability', 'potency', 'solubility', 'stability'] else None)
            for m in mmp:
                if m['smiles'] not in seen_smiles and m['similarity'] >= min_similarity:
                    seen_smiles.add(m['smiles'])
                    all_analogs.append({
                        'smiles': m['smiles'],
                        'mol': m['mol'],
                        'strategy': 'matched_pair',
                        'detail': m['transform_name'],
                        'expected_effect': m['expected_effect'],
                        'similarity': m['similarity']
                    })
        
        scored_analogs = []
        for analog in all_analogs:
            try:
                score = self.multi_objective_score(analog['mol'], weights=weights)
                
                scored_analogs.append({
                    'smiles': analog['smiles'],
                    'strategy': analog['strategy'],
                    'modification': analog['detail'],
                    'expected_effect': analog.get('expected_effect', 'general'),
                    'similarity': round(analog['similarity'], 3),
                    'combined_score': score['combined'],
                    'activity_score': round(score['activity'], 3),
                    'drug_likeness_score': round(score['drug_likeness'], 3),
                    'safety_score': round(score['safety'], 3),
                    'synthesizability_score': round(score['synthesizability'], 3),
                    'top_receptor': score.get('top_receptor', 'Unknown'),
                    'pains_alerts': score.get('pains_alerts', []),
                    'properties': score['properties'],
                    'improvement_vs_parent': round(score['combined'] - parent_score['combined'], 4)
                })
            except Exception as e:
                continue
        
        scored_analogs.sort(key=lambda x: (-x['combined_score'], -x['similarity']))
        final_analogs = scored_analogs[:num_analogs]
        
        improved_count = len([a for a in final_analogs if a['improvement_vs_parent'] > 0])
        strategy_counts = defaultdict(int)
        for a in final_analogs:
            strategy_counts[a['strategy']] += 1
        
        avg_improvement = np.mean([a['improvement_vs_parent'] for a in final_analogs]) if final_analogs else 0
        
        return {
            'parent': {
                'smiles': canonical,
                'properties': parent_props,
                'multi_objective_score': parent_score['combined'],
                'scaffold': parent_framework,
                'generic_scaffold': parent_generic
            },
            'analogs': final_analogs,
            'optimization_summary': {
                'target': optimization_target,
                'weights_used': weights,
                'total_generated': len(all_analogs),
                'returned': len(final_analogs),
                'improved_count': improved_count,
                'improvement_rate': round(improved_count / len(final_analogs) * 100, 1) if final_analogs else 0,
                'average_improvement': round(avg_improvement, 4),
                'strategies_used': dict(strategy_counts),
                'min_similarity_used': min_similarity
            },
            'therapeutic_context': therapeutic_area,
            'timestamp': datetime.now().isoformat()
        }
    
    def screen_analog_series(self, smiles_list: List[str],
                              reference_smiles: Optional[str] = None) -> Dict:
        """
        Screen a series of analogs against multiple objectives.
        Useful for evaluating custom analog libraries.
        """
        results = []
        reference_mol = None
        
        if reference_smiles:
            valid, reference_mol, _ = self.validate_molecule(reference_smiles)
        
        for smiles in smiles_list[:100]:
            valid, mol, canonical = self.validate_molecule(smiles)
            if not valid:
                results.append({
                    'smiles': smiles,
                    'valid': False,
                    'error': canonical
                })
                continue
            
            score = self.multi_objective_score(mol)
            similarity = self.calculate_similarity(reference_mol, mol) if reference_mol else None
            
            results.append({
                'smiles': canonical,
                'valid': True,
                'combined_score': score['combined'],
                'activity': round(score['activity'], 3),
                'drug_likeness': round(score['drug_likeness'], 3),
                'safety': round(score['safety'], 3),
                'synthesizability': round(score['synthesizability'], 3),
                'similarity_to_reference': round(similarity, 3) if similarity else None,
                'properties': score['properties'],
                'pains_alerts': score.get('pains_alerts', [])
            })
        
        results.sort(key=lambda x: -x.get('combined_score', 0) if x['valid'] else float('-inf'))
        
        valid_results = [r for r in results if r['valid']]
        
        return {
            'results': results,
            'summary': {
                'total_screened': len(smiles_list),
                'valid_compounds': len(valid_results),
                'invalid_compounds': len(smiles_list) - len(valid_results),
                'top_score': max([r['combined_score'] for r in valid_results]) if valid_results else 0,
                'average_score': round(np.mean([r['combined_score'] for r in valid_results]), 4) if valid_results else 0,
                'pains_flagged': len([r for r in valid_results if r.get('pains_alerts')])
            },
            'reference_smiles': reference_smiles,
            'timestamp': datetime.now().isoformat()
        }


advanced_generator = AdvancedAnalogGenerator()


def generate_optimized_analogs(smiles: str, **kwargs) -> Dict:
    """Convenience function for analog generation"""
    return advanced_generator.generate_optimized_analogs(smiles, **kwargs)


def screen_analog_series(smiles_list: List[str], reference: str = None) -> Dict:
    """Convenience function for series screening"""
    return advanced_generator.screen_analog_series(smiles_list, reference)
