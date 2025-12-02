#!/usr/bin/env python3
"""
Test script for PharmaSight™ viability analysis system
Tests the enhanced compound-analysis service with SA Score, QED, and viability scoring
"""

import requests
import json
from typing import Dict, Any

# Test compounds with known properties
TEST_COMPOUNDS = {
    "psilocybin": {
        "smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "name": "Psilocybin",
        "expected_qed": "> 0.5",  # Should have good drug-likeness
        "expected_sa": "< 5",  # Should be moderately easy to synthesize
    },
    "aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "name": "Aspirin",
        "expected_qed": "> 0.6",  # Well-known drug
        "expected_sa": "< 3",  # Very easy to synthesize
    },
    "caffeine": {
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "name": "Caffeine",
        "expected_qed": "> 0.5",
        "expected_sa": "< 4",
    },
}

def test_health_check(base_url: str = "http://localhost:8001"):
    """Test the health endpoint."""
    print("\n" + "="*60)
    print("TESTING HEALTH CHECK")
    print("="*60)
    
    try:
        response = requests.get(f"{base_url}/health")
        response.raise_for_status()
        data = response.json()
        
        print(f"✓ Service Status: {data['status']}")
        print(f"✓ Service Name: {data['service']}")
        print(f"✓ RDKit Version: {data['rdkit_version']}")
        print(f"✓ Redis Status: {data['redis']}")
        print("\nAvailable Features:")
        for feature, available in data.get('features', {}).items():
            status = "✓" if available else "✗"
            print(f"  {status} {feature}: {available}")
        
        return True
    except Exception as e:
        print(f"✗ Health check failed: {e}")
        return False

def test_basic_analysis(base_url: str = "http://localhost:8001"):
    """Test basic compound analysis."""
    print("\n" + "="*60)
    print("TESTING BASIC COMPOUND ANALYSIS")
    print("="*60)
    
    for compound_id, compound_info in TEST_COMPOUNDS.items():
        print(f"\nTesting: {compound_info['name']}")
        print(f"SMILES: {compound_info['smiles']}")
        
        try:
            response = requests.post(
                f"{base_url}/analyze",
                json={
                    "compound": compound_info['smiles'],
                    "include_viability": False
                }
            )
            response.raise_for_status()
            data = response.json()
            
            rdkit_analysis = data.get('rdkit_analysis', {})
            print(f"  ✓ Molecular Weight: {rdkit_analysis.get('molecular_weight')} Da")
            print(f"  ✓ LogP: {rdkit_analysis.get('logp')}")
            print(f"  ✓ TPSA: {rdkit_analysis.get('tpsa')} Ų")
            print(f"  ✓ H-Bond Donors: {rdkit_analysis.get('num_h_donors')}")
            print(f"  ✓ H-Bond Acceptors: {rdkit_analysis.get('num_h_acceptors')}")
            
        except Exception as e:
            print(f"  ✗ Analysis failed: {e}")

def test_viability_analysis(base_url: str = "http://localhost:8001"):
    """Test viability analysis with SA Score, QED, and viability scoring."""
    print("\n" + "="*60)
    print("TESTING VIABILITY ANALYSIS")
    print("="*60)
    
    for compound_id, compound_info in TEST_COMPOUNDS.items():
        print(f"\nAnalyzing: {compound_info['name']}")
        print(f"SMILES: {compound_info['smiles']}")
        
        try:
            response = requests.post(
                f"{base_url}/viability",
                json={
                    "smiles": compound_info['smiles'],
                    "include_3d": False
                }
            )
            response.raise_for_status()
            data = response.json()
            
            viability = data.get('viability_analysis', {})
            
            # Synthetic Accessibility
            sa_data = viability.get('synthetic_accessibility', {})
            print(f"\n  Synthetic Accessibility:")
            print(f"    Score: {sa_data.get('sa_score')} (1=easy, 10=difficult)")
            print(f"    Feasibility: {sa_data.get('synthetic_feasibility')}")
            
            # Drug-likeness (QED)
            qed_data = viability.get('drug_likeness', {})
            print(f"\n  Drug-likeness (QED):")
            print(f"    Score: {qed_data.get('qed_score')} (0-1, higher=better)")
            print(f"    Interpretation: {qed_data.get('qed_interpretation')}")
            
            # Lipinski's Rule of Five
            ro5_data = viability.get('lipinski_rule_of_five', {})
            print(f"\n  Lipinski's Rule of Five:")
            print(f"    Violations: {ro5_data.get('violations')}")
            print(f"    Passes: {ro5_data.get('passes_ro5')}")
            print(f"    {ro5_data.get('interpretation')}")
            
            # Overall Viability Score
            viability_score = viability.get('viability_score', {})
            print(f"\n  Overall Viability:")
            print(f"    Score: {viability_score.get('overall_score')}/100")
            print(f"    Priority: {viability_score.get('priority')}")
            print(f"    Recommendation: {viability_score.get('recommendation')}")
            
            # Component scores
            print(f"\n  Component Scores:")
            for component, score in viability_score.get('component_scores', {}).items():
                print(f"    {component}: {score}/100")
            
        except Exception as e:
            print(f"  ✗ Viability analysis failed: {e}")
            import traceback
            traceback.print_exc()

def test_batch_viability(base_url: str = "http://localhost:8001"):
    """Test batch viability analysis."""
    print("\n" + "="*60)
    print("TESTING BATCH VIABILITY ANALYSIS")
    print("="*60)
    
    smiles_list = [info['smiles'] for info in TEST_COMPOUNDS.values()]
    
    try:
        response = requests.post(
            f"{base_url}/batch_viability",
            json=smiles_list
        )
        response.raise_for_status()
        results = response.json()
        
        print(f"\nAnalyzed {len(results)} compounds")
        print("\nRanked by Viability Score:")
        print("-" * 60)
        
        for i, result in enumerate(results, 1):
            if 'error' in result:
                print(f"{i}. ERROR: {result['error']}")
            else:
                print(f"{i}. Score: {result['viability_score']}/100 | Priority: {result['priority']}")
                print(f"   SMILES: {result['smiles']}")
        
    except Exception as e:
        print(f"✗ Batch analysis failed: {e}")

def generate_viability_report(base_url: str = "http://localhost:8001"):
    """Generate a comprehensive viability report for all test compounds."""
    print("\n" + "="*60)
    print("GENERATING COMPREHENSIVE VIABILITY REPORT")
    print("="*60)
    
    report = {
        "timestamp": "2025-12-01",
        "compounds_analyzed": len(TEST_COMPOUNDS),
        "results": []
    }
    
    for compound_id, compound_info in TEST_COMPOUNDS.items():
        try:
            response = requests.post(
                f"{base_url}/viability",
                json={
                    "smiles": compound_info['smiles'],
                    "include_3d": False
                }
            )
            response.raise_for_status()
            data = response.json()
            
            report["results"].append({
                "name": compound_info['name'],
                "smiles": compound_info['smiles'],
                "viability_analysis": data.get('viability_analysis')
            })
            
        except Exception as e:
            report["results"].append({
                "name": compound_info['name'],
                "smiles": compound_info['smiles'],
                "error": str(e)
            })
    
    # Save report to file
    report_file = "/home/ubuntu/pharmasight-platform/viability_test_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n✓ Report saved to: {report_file}")
    print(f"✓ Total compounds analyzed: {len(report['results'])}")

def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("PHARMASIGHT™ VIABILITY ANALYSIS TEST SUITE")
    print("="*60)
    
    base_url = "http://localhost:8001"
    
    # Test 1: Health Check
    if not test_health_check(base_url):
        print("\n✗ Service not available. Make sure Docker Compose is running:")
        print("  docker-compose up compound-service")
        return
    
    # Test 2: Basic Analysis
    test_basic_analysis(base_url)
    
    # Test 3: Viability Analysis
    test_viability_analysis(base_url)
    
    # Test 4: Batch Analysis
    test_batch_viability(base_url)
    
    # Test 5: Generate Report
    generate_viability_report(base_url)
    
    print("\n" + "="*60)
    print("TEST SUITE COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
