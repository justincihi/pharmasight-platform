"""
Comprehensive Molecular Analysis Wrapper
Combines all Phase I & II analysis modules
"""

import json
import sys
from typing import Dict, List, Optional

# Import all analysis modules
from toxicity_profiler import predict_toxicity_profile
from synthetic_accessibility import calculate_sa_score
from structure_optimizer import optimize_structure, iterative_optimization

def run_comprehensive_analysis(smiles: str) -> Dict:
    """
    Run all analyses on a molecule
    """
    try:
        results = {
            'smiles': smiles,
            'toxicity_profile': predict_toxicity_profile(smiles),
            'synthetic_accessibility': calculate_sa_score(smiles),
            'optimization_suggestions': optimize_structure(smiles)[:5],  # Top 5 suggestions
            'status': 'success'
        }
        return results
    except Exception as e:
        return {
            'smiles': smiles,
            'status': 'error',
            'error': str(e)
        }

def run_toxicity_analysis(smiles: str) -> Dict:
    """
    Run only toxicity profiling
    """
    try:
        return predict_toxicity_profile(smiles)
    except Exception as e:
        return {'error': str(e)}

def run_sa_analysis(smiles: str) -> Dict:
    """
    Run only synthetic accessibility analysis
    """
    try:
        return calculate_sa_score(smiles)
    except Exception as e:
        return {'error': str(e)}

def run_optimization_analysis(smiles: str, target_property: Optional[str] = None) -> List[Dict]:
    """
    Run structure optimization analysis
    """
    try:
        return optimize_structure(smiles, target_property)
    except Exception as e:
        return [{'error': str(e)}]

def run_iterative_optimization(smiles: str, target_property: str, iterations: int = 3) -> List[Dict]:
    """
    Run iterative optimization
    """
    try:
        return iterative_optimization(smiles, target_property, iterations)
    except Exception as e:
        return [{'error': str(e)}]

def batch_analysis(smiles_list: List[str]) -> List[Dict]:
    """
    Run comprehensive analysis on multiple molecules
    """
    results = []
    for smiles in smiles_list:
        result = run_comprehensive_analysis(smiles)
        results.append(result)
    return results

if __name__ == '__main__':
    # CLI interface for Node.js integration
    if len(sys.argv) < 3:
        print(json.dumps({'error': 'Usage: python comprehensive_analysis.py <command> <smiles> [options]'}))
        sys.exit(1)
    
    command = sys.argv[1]
    smiles = sys.argv[2]
    
    if command == 'comprehensive':
        result = run_comprehensive_analysis(smiles)
    elif command == 'toxicity':
        result = run_toxicity_analysis(smiles)
    elif command == 'sa_score':
        result = run_sa_analysis(smiles)
    elif command == 'optimize':
        target = sys.argv[3] if len(sys.argv) > 3 else None
        result = run_optimization_analysis(smiles, target)
    elif command == 'iterative_optimize':
        target = sys.argv[3] if len(sys.argv) > 3 else 'reduce_lipophilicity'
        iterations = int(sys.argv[4]) if len(sys.argv) > 4 else 3
        result = run_iterative_optimization(smiles, target, iterations)
    else:
        result = {'error': f'Unknown command: {command}'}
    
    print(json.dumps(result, indent=2))
