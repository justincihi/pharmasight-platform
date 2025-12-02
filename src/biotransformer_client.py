"""
BioTransformer 3.0 Integration Client
Provides metabolism prediction via BioTransformer REST API

BioTransformer is open source (GPL v2.1):
https://github.com/Wishartlab-openscience/Biotransformer

API Modes:
- CYP450: Cytochrome P450 Phase I metabolism
- PHASEII: Phase II conjugation reactions (glucuronidation, sulfation, etc.)
- HGUT: Human gut microbiota metabolism
- ALLHUMAN: Combined human tissue + gut (recommended)
- SUPERBIO: Comprehensive 4-iteration metabolism prediction
- ENVMICRO: Environmental microbial degradation
"""

import os
import json
import time
import hashlib
import requests
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
from rdkit import Chem
from rdkit.Chem import Descriptors

BIOTRANSFORMER_API_BASE = "https://biotransformer.ca"
RATE_LIMIT_DELAY = 31
MAX_POLL_ATTEMPTS = 60
POLL_INTERVAL = 5

BIOTRANSFORMER_MODES = {
    'cyp450': 'CYP450',
    'phase2': 'PHASEII', 
    'phaseii': 'PHASEII',
    'gut': 'HGUT',
    'hgut': 'HGUT',
    'allhuman': 'ALLHUMAN',
    'human': 'ALLHUMAN',
    'superbio': 'SUPERBIO',
    'comprehensive': 'SUPERBIO',
    'envmicro': 'ENVMICRO',
    'environmental': 'ENVMICRO',
    'ecbased': 'EC-BASED',
}

_prediction_cache: Dict[str, Dict] = {}
_last_request_time: float = 0


def get_cache_key(smiles: str, mode: str) -> str:
    """Generate cache key from SMILES and mode"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            canonical = Chem.MolToSmiles(mol, canonical=True)
        else:
            canonical = smiles
    except:
        canonical = smiles
    
    key_str = f"{canonical}:{mode.upper()}"
    return hashlib.md5(key_str.encode()).hexdigest()


def get_cached_prediction(smiles: str, mode: str) -> Optional[Dict]:
    """Get cached prediction if available and not expired"""
    cache_key = get_cache_key(smiles, mode)
    if cache_key in _prediction_cache:
        cached = _prediction_cache[cache_key]
        cache_time = cached.get('timestamp', 0)
        if time.time() - cache_time < 86400:
            return cached.get('result')
    return None


def cache_prediction(smiles: str, mode: str, result: Dict):
    """Cache a prediction result"""
    cache_key = get_cache_key(smiles, mode)
    _prediction_cache[cache_key] = {
        'result': result,
        'timestamp': time.time()
    }


def rate_limit_wait():
    """Enforce rate limiting (2 requests per minute)"""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < RATE_LIMIT_DELAY:
        wait_time = RATE_LIMIT_DELAY - elapsed
        time.sleep(wait_time)
    _last_request_time = time.time()


def validate_smiles(smiles: str) -> Dict[str, Any]:
    """Validate SMILES and check molecular weight constraint"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'valid': False, 'error': 'Invalid SMILES string'}
        
        mw = Descriptors.MolWt(mol)
        if mw > 1500:
            return {
                'valid': False, 
                'error': f'Molecular weight ({mw:.1f} Da) exceeds BioTransformer limit of 1500 Da'
            }
        
        return {
            'valid': True,
            'canonical_smiles': Chem.MolToSmiles(mol, canonical=True),
            'molecular_weight': round(mw, 2)
        }
    except Exception as e:
        return {'valid': False, 'error': str(e)}


def submit_prediction(smiles: str, mode: str = 'ALLHUMAN', label: str = None) -> Dict[str, Any]:
    """
    Submit metabolism prediction request to BioTransformer API
    
    Args:
        smiles: SMILES string of compound
        mode: Prediction mode (CYP450, PHASEII, HGUT, ALLHUMAN, SUPERBIO, ENVMICRO)
        label: Optional query label
    
    Returns:
        Dict with query_id for polling, or error
    """
    mode_upper = BIOTRANSFORMER_MODES.get(mode.lower(), mode.upper())
    
    validation = validate_smiles(smiles)
    if not validation['valid']:
        return {'error': validation['error'], 'status': 'failed'}
    
    cached = get_cached_prediction(smiles, mode_upper)
    if cached:
        return {
            'status': 'completed',
            'cached': True,
            'result': cached
        }
    
    rate_limit_wait()
    
    try:
        data = {
            'biotransformer_option': mode_upper,
            'query_input': validation['canonical_smiles']
        }
        if label:
            data['query_label'] = label
        
        response = requests.post(
            f"{BIOTRANSFORMER_API_BASE}/queries.json",
            data=data,
            headers={'Accept': 'application/json'},
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            query_id = result.get('id')
            if query_id:
                return {
                    'status': 'submitted',
                    'query_id': query_id,
                    'mode': mode_upper,
                    'smiles': validation['canonical_smiles'],
                    'message': 'Prediction submitted. Poll for results.'
                }
            else:
                return {'error': 'No query ID returned', 'status': 'failed'}
        
        elif response.status_code == 429:
            return {
                'error': 'Rate limit exceeded. Please wait and retry.',
                'status': 'rate_limited',
                'retry_after': 60
            }
        
        else:
            return {
                'error': f'API error: {response.status_code}',
                'status': 'failed',
                'details': response.text[:500]
            }
    
    except requests.exceptions.Timeout:
        return {'error': 'Request timeout', 'status': 'timeout'}
    except requests.exceptions.RequestException as e:
        return {'error': f'Connection error: {str(e)}', 'status': 'failed'}


def get_prediction_status(query_id: int) -> Dict[str, Any]:
    """
    Poll for prediction results
    
    Args:
        query_id: Query ID from submit_prediction
    
    Returns:
        Dict with prediction status and results if complete
    """
    try:
        response = requests.get(
            f"{BIOTRANSFORMER_API_BASE}/queries/{query_id}.json",
            headers={'Accept': 'application/json'},
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            
            if 'prediction_errors' in data and data['prediction_errors']:
                return {
                    'status': 'failed',
                    'error': data['prediction_errors']
                }
            
            if 'metabolites' in data or 'results' in data:
                return {
                    'status': 'completed',
                    'result': data
                }
            
            return {
                'status': 'processing',
                'message': 'Prediction still in progress'
            }
        
        elif response.status_code == 404:
            return {'status': 'not_found', 'error': 'Query not found'}
        
        else:
            return {
                'status': 'error',
                'error': f'API error: {response.status_code}'
            }
    
    except requests.exceptions.RequestException as e:
        return {'status': 'error', 'error': str(e)}


def predict_metabolism_sync(
    smiles: str, 
    mode: str = 'ALLHUMAN',
    max_wait: int = 300,
    label: str = None
) -> Dict[str, Any]:
    """
    Synchronous metabolism prediction with polling
    
    Args:
        smiles: SMILES string
        mode: Prediction mode
        max_wait: Maximum wait time in seconds (default 5 minutes)
        label: Optional query label
    
    Returns:
        Complete prediction result or error
    """
    mode_upper = BIOTRANSFORMER_MODES.get(mode.lower(), mode.upper())
    
    cached = get_cached_prediction(smiles, mode_upper)
    if cached:
        return {
            'success': True,
            'cached': True,
            'mode': mode_upper,
            'metabolites': cached.get('metabolites', []),
            'reactions': cached.get('reactions', []),
            'raw_result': cached
        }
    
    submit_result = submit_prediction(smiles, mode, label)
    
    if submit_result.get('status') == 'completed':
        return {
            'success': True,
            'cached': True,
            'mode': mode_upper,
            'metabolites': submit_result['result'].get('metabolites', []),
            'reactions': submit_result['result'].get('reactions', []),
            'raw_result': submit_result['result']
        }
    
    if submit_result.get('status') != 'submitted':
        return {
            'success': False,
            'error': submit_result.get('error', 'Unknown error'),
            'mode': mode_upper
        }
    
    query_id = submit_result['query_id']
    start_time = time.time()
    poll_count = 0
    
    while (time.time() - start_time) < max_wait:
        time.sleep(POLL_INTERVAL)
        poll_count += 1
        
        status = get_prediction_status(query_id)
        
        if status['status'] == 'completed':
            result = status['result']
            cache_prediction(smiles, mode_upper, result)
            
            metabolites = parse_metabolites(result)
            
            return {
                'success': True,
                'cached': False,
                'mode': mode_upper,
                'query_id': query_id,
                'poll_count': poll_count,
                'processing_time': round(time.time() - start_time, 1),
                'metabolites': metabolites,
                'metabolite_count': len(metabolites),
                'raw_result': result
            }
        
        elif status['status'] == 'failed':
            return {
                'success': False,
                'error': status.get('error', 'Prediction failed'),
                'mode': mode_upper
            }
    
    return {
        'success': False,
        'error': f'Prediction timeout after {max_wait} seconds',
        'mode': mode_upper,
        'query_id': query_id,
        'message': 'Use query_id to poll for results later'
    }


def parse_metabolites(result: Dict) -> List[Dict]:
    """Parse metabolites from BioTransformer response"""
    metabolites = []
    
    raw_metabolites = result.get('metabolites', [])
    if not raw_metabolites:
        raw_metabolites = result.get('products', [])
    
    for met in raw_metabolites:
        if isinstance(met, dict):
            metabolite = {
                'smiles': met.get('smiles', met.get('SMILES', '')),
                'name': met.get('name', met.get('Common Name', 'Unknown')),
                'inchikey': met.get('inchikey', met.get('InChIKey', '')),
                'molecular_formula': met.get('molecular_formula', met.get('Molecular Formula', '')),
                'molecular_weight': met.get('molecular_weight', met.get('Molecular Weight', 0)),
                'reaction_type': met.get('reaction_type', met.get('Reaction', '')),
                'enzyme': met.get('enzyme', met.get('Enzyme(s)', '')),
                'biosystem': met.get('biosystem', met.get('Biosystem', '')),
            }
            metabolites.append(metabolite)
    
    return metabolites


def predict_cyp450(smiles: str, max_wait: int = 180) -> Dict[str, Any]:
    """Predict CYP450 Phase I metabolism"""
    return predict_metabolism_sync(smiles, 'CYP450', max_wait)


def predict_phase2(smiles: str, max_wait: int = 180) -> Dict[str, Any]:
    """Predict Phase II conjugation metabolism"""
    return predict_metabolism_sync(smiles, 'PHASEII', max_wait)


def predict_gut_metabolism(smiles: str, max_wait: int = 180) -> Dict[str, Any]:
    """Predict gut microbiota metabolism"""
    return predict_metabolism_sync(smiles, 'HGUT', max_wait)


def predict_human_metabolism(smiles: str, max_wait: int = 300) -> Dict[str, Any]:
    """Predict all human metabolism (liver + gut)"""
    return predict_metabolism_sync(smiles, 'ALLHUMAN', max_wait)


def predict_comprehensive(smiles: str, max_wait: int = 600) -> Dict[str, Any]:
    """Comprehensive metabolism prediction (SUPERBIO - 4 iterations)"""
    return predict_metabolism_sync(smiles, 'SUPERBIO', max_wait)


def predict_environmental(smiles: str, max_wait: int = 180) -> Dict[str, Any]:
    """Predict environmental microbial degradation"""
    return predict_metabolism_sync(smiles, 'ENVMICRO', max_wait)


def get_metabolism_summary(smiles: str) -> Dict[str, Any]:
    """
    Get a quick metabolism summary using local predictions
    Falls back to heuristic predictions if API unavailable
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        cyp_substrates = []
        cyp_smarts = {
            'CYP1A2': ['c1ccc2c(c1)cc[nH]2', 'c1ccc2c(c1)[nH]c3ccccc23'],
            'CYP2C9': ['c1ccc(cc1)C(=O)O', 'c1ccc(cc1)NS(=O)(=O)'],
            'CYP2C19': ['c1ccc2c(c1)[nH]cn2', 'c1ccc(cc1)C(=O)N'],
            'CYP2D6': ['c1ccc(cc1)CCN', 'c1ccc(cc1)CN'],
            'CYP3A4': ['C=O', 'CC(C)C'],
        }
        
        for cyp, patterns in cyp_smarts.items():
            for pattern in patterns:
                try:
                    patt = Chem.MolFromSmarts(pattern)
                    if patt and mol.HasSubstructMatch(patt):
                        cyp_substrates.append(cyp)
                        break
                except:
                    continue
        
        phase2_reactions = []
        phase2_smarts = {
            'Glucuronidation': ['[OH]', '[NH2]', 'C(=O)O'],
            'Sulfation': ['[OH]', '[NH2]'],
            'Acetylation': ['[NH2]', '[NH]'],
            'Methylation': ['[OH]', '[NH2]', '[SH]'],
            'Glutathione conjugation': ['C=C', '[Cl,Br,I]', 'C(=O)C'],
        }
        
        for reaction, patterns in phase2_smarts.items():
            for pattern in patterns:
                try:
                    patt = Chem.MolFromSmarts(pattern)
                    if patt and mol.HasSubstructMatch(patt):
                        phase2_reactions.append(reaction)
                        break
                except:
                    continue
        
        metabolism_rate = 'Moderate'
        if logp > 3 and mw < 500:
            metabolism_rate = 'High'
        elif logp < 0 or mw > 700:
            metabolism_rate = 'Low'
        
        return {
            'smiles': smiles,
            'predicted_cyp_substrates': list(set(cyp_substrates)),
            'predicted_phase2_reactions': list(set(phase2_reactions)),
            'estimated_metabolism_rate': metabolism_rate,
            'properties': {
                'molecular_weight': round(mw, 2),
                'logP': round(logp, 2),
                'tpsa': round(tpsa, 2),
                'hba': hba,
                'hbd': hbd,
                'rotatable_bonds': rotatable_bonds
            },
            'note': 'Local heuristic prediction. Use full API for detailed metabolite structures.'
        }
    
    except Exception as e:
        return {'error': str(e)}


class BioTransformerClient:
    """High-level client for BioTransformer API integration"""
    
    def __init__(self):
        self.base_url = BIOTRANSFORMER_API_BASE
        self.cache = {}
        self.job_queue = {}
    
    def predict(self, smiles: str, mode: str = 'allhuman', sync: bool = True, max_wait: int = 300):
        """
        Main prediction method
        
        Args:
            smiles: SMILES string
            mode: cyp450, phase2, gut, allhuman, superbio, envmicro
            sync: If True, wait for results. If False, return job ID.
            max_wait: Max seconds to wait for sync mode
        """
        if sync:
            return predict_metabolism_sync(smiles, mode, max_wait)
        else:
            result = submit_prediction(smiles, mode)
            if result.get('status') == 'submitted':
                job_id = result['query_id']
                self.job_queue[job_id] = {
                    'smiles': smiles,
                    'mode': mode,
                    'submitted': time.time()
                }
            return result
    
    def get_status(self, job_id: int):
        """Get status of async prediction"""
        return get_prediction_status(job_id)
    
    def get_summary(self, smiles: str):
        """Get quick metabolism summary (local prediction)"""
        return get_metabolism_summary(smiles)
    
    def predict_all_pathways(self, smiles: str, max_wait: int = 600):
        """Predict metabolism through all pathways (comprehensive)"""
        return predict_comprehensive(smiles, max_wait)


biotransformer = BioTransformerClient()
