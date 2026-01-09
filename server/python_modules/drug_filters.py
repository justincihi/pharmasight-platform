"""
Drug Filters and Scoring for PharmaSight™
Implements PAINS/Brenk filters, CNS MPO scoring, and BBB permeability prediction
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, FilterCatalog
from typing import Dict, List, Tuple
import math


def check_pains_brenk_filters(smiles: str) -> Dict:
    """
    Check molecule against PAINS and Brenk structural alert filters
    
    Args:
        smiles: SMILES string of molecule
    
    Returns:
        Dictionary with filter results and alerts
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'error': 'Invalid SMILES'}
    
    results = {
        'smiles': smiles,
        'passes_pains': True,
        'passes_brenk': True,
        'pains_alerts': [],
        'brenk_alerts': [],
        'total_alerts': 0,
        'recommendation': 'PASS',
    }
    
    try:
        # PAINS filter
        params_pains = FilterCatalog.FilterCatalogParams()
        params_pains.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog_pains = FilterCatalog.FilterCatalog(params_pains)
        
        pains_matches = catalog_pains.GetMatches(mol)
        if pains_matches:
            results['passes_pains'] = False
            results['pains_alerts'] = [match.GetDescription() for match in pains_matches]
        
        # Brenk filter
        params_brenk = FilterCatalog.FilterCatalogParams()
        params_brenk.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        catalog_brenk = FilterCatalog.FilterCatalog(params_brenk)
        
        brenk_matches = catalog_brenk.GetMatches(mol)
        if brenk_matches:
            results['passes_brenk'] = False
            results['brenk_alerts'] = [match.GetDescription() for match in brenk_matches]
        
        results['total_alerts'] = len(results['pains_alerts']) + len(results['brenk_alerts'])
        
        if results['total_alerts'] > 0:
            results['recommendation'] = 'FAIL - Contains structural alerts'
        
    except Exception as e:
        results['error'] = str(e)
    
    return results


def calculate_cns_mpo_score(smiles: str) -> Dict:
    """
    Calculate CNS Multiparameter Optimization (MPO) score
    
    Based on Wager et al. (2010) ACS Chem Neurosci 1(6):435-449
    Scores 6 physicochemical properties on a 0-1 scale, sum = CNS MPO (0-6)
    
    Args:
        smiles: SMILES string of molecule
    
    Returns:
        Dictionary with CNS MPO score and component scores
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'error': 'Invalid SMILES'}
    
    # Calculate descriptors
    clogp = Descriptors.MolLogP(mol)
    clogd = clogp  # Approximation (true cLogD requires pH)
    mw = Descriptors.MolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    pka = 8.0  # Default assumption for basic compounds
    
    # CNS MPO scoring functions (0-1 for each property)
    def score_clogp(val):
        if val <= 3:
            return 1.0
        elif val >= 5:
            return 0.0
        else:
            return 1.0 - (val - 3) / 2
    
    def score_clogd(val):
        if val <= 2:
            return 1.0
        elif val >= 4:
            return 0.0
        else:
            return 1.0 - (val - 2) / 2
    
    def score_mw(val):
        if val <= 360:
            return 1.0
        elif val >= 500:
            return 0.0
        else:
            return 1.0 - (val - 360) / 140
    
    def score_tpsa(val):
        if 40 <= val <= 90:
            return 1.0
        elif val < 20 or val > 120:
            return 0.0
        elif val < 40:
            return (val - 20) / 20
        else:
            return 1.0 - (val - 90) / 30
    
    def score_hbd(val):
        if val <= 0.5:
            return 1.0
        elif val >= 3.5:
            return 0.0
        else:
            return 1.0 - (val - 0.5) / 3
    
    def score_pka(val):
        if 7 <= val <= 10:
            return 1.0
        elif val < 4 or val > 11:
            return 0.0
        elif val < 7:
            return (val - 4) / 3
        else:
            return 1.0 - (val - 10) / 1
    
    # Calculate component scores
    scores = {
        'clogp': round(score_clogp(clogp), 2),
        'clogd': round(score_clogd(clogd), 2),
        'mw': round(score_mw(mw), 2),
        'tpsa': round(score_tpsa(tpsa), 2),
        'hbd': round(score_hbd(hbd), 2),
        'pka': round(score_pka(pka), 2),
    }
    
    cns_mpo = sum(scores.values())
    
    # Interpretation
    if cns_mpo >= 4.0:
        interpretation = 'Excellent CNS drug-likeness'
    elif cns_mpo >= 3.0:
        interpretation = 'Good CNS drug-likeness'
    elif cns_mpo >= 2.0:
        interpretation = 'Moderate CNS drug-likeness'
    else:
        interpretation = 'Poor CNS drug-likeness'
    
    return {
        'smiles': smiles,
        'cns_mpo_score': round(cns_mpo, 2),
        'max_score': 6.0,
        'component_scores': scores,
        'descriptors': {
            'clogp': round(clogp, 2),
            'clogd': round(clogd, 2),
            'mw': round(mw, 2),
            'tpsa': round(tpsa, 2),
            'hbd': hbd,
            'pka': pka,
        },
        'interpretation': interpretation,
    }


def predict_bbb_permeability(smiles: str) -> Dict:
    """
    Predict blood-brain barrier (BBB) permeability using heuristic rules
    
    Based on multiple literature sources:
    - Lipinski-like rules for CNS penetration
    - TPSA < 90 Å² (Ertl et al.)
    - MW < 450 Da
    - LogP 1-5
    - HBD ≤ 3
    
    Args:
        smiles: SMILES string of molecule
    
    Returns:
        Dictionary with BBB permeability prediction
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'error': 'Invalid SMILES'}
    
    # Calculate descriptors
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # BBB permeability rules
    rules = {
        'mw_ok': mw <= 450,
        'logp_ok': 1 <= logp <= 5,
        'tpsa_ok': tpsa <= 90,
        'hbd_ok': hbd <= 3,
        'hba_ok': hba <= 9,
        'rotatable_bonds_ok': rotatable_bonds <= 8,
    }
    
    passes = sum(rules.values())
    total_rules = len(rules)
    
    # Prediction
    if passes == total_rules:
        prediction = 'HIGH'
        probability = 0.85
        recommendation = 'Likely to cross BBB'
    elif passes >= total_rules - 1:
        prediction = 'MODERATE'
        probability = 0.60
        recommendation = 'May cross BBB with optimization'
    elif passes >= total_rules - 2:
        prediction = 'LOW'
        probability = 0.30
        recommendation = 'Unlikely to cross BBB'
    else:
        prediction = 'VERY LOW'
        probability = 0.10
        recommendation = 'Very unlikely to cross BBB'
    
    return {
        'smiles': smiles,
        'bbb_permeability': prediction,
        'probability': probability,
        'recommendation': recommendation,
        'rules_passed': f'{passes}/{total_rules}',
        'rule_details': rules,
        'descriptors': {
            'mw': round(mw, 2),
            'logp': round(logp, 2),
            'tpsa': round(tpsa, 2),
            'hbd': hbd,
            'hba': hba,
            'rotatable_bonds': rotatable_bonds,
        },
    }


def comprehensive_drug_assessment(smiles: str) -> Dict:
    """
    Run all drug filters and scores on a molecule
    
    Args:
        smiles: SMILES string of molecule
    
    Returns:
        Dictionary with all assessment results
    """
    return {
        'smiles': smiles,
        'pains_brenk': check_pains_brenk_filters(smiles),
        'cns_mpo': calculate_cns_mpo_score(smiles),
        'bbb_permeability': predict_bbb_permeability(smiles),
    }


# Test function
if __name__ == '__main__':
    # Test with Aspirin
    test_smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    print('Testing with Aspirin:', test_smiles)
    print('\n=== PAINS/Brenk Filters ===')
    print(check_pains_brenk_filters(test_smiles))
    print('\n=== CNS MPO Score ===')
    print(calculate_cns_mpo_score(test_smiles))
    print('\n=== BBB Permeability ===')
    print(predict_bbb_permeability(test_smiles))
