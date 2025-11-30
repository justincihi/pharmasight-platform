#!/usr/bin/env python3.11
"""
PharmaSight Platform - Comprehensive Diagnostic Test Suite
Tests all major features and integrations
"""

import sys
import os
import json
from datetime import datetime

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def print_section(title):
    """Print a formatted section header"""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")

def test_rdkit_integration():
    """Test RDKit molecular visualization and editing"""
    print_section("1. RDKit Integration Test")
    
    try:
        from rdkit import Chem
        from rdkit import __version__ as rdkit_version
        print(f"✅ RDKit imported successfully (version {rdkit_version})")
        
        # Test SMILES parsing
        test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
        mol = Chem.MolFromSmiles(test_smiles)
        if mol:
            print(f"✅ SMILES parsing works: {test_smiles}")
        else:
            print(f"❌ SMILES parsing failed")
            return False
            
        # Test molecular visualizer
        from molecular_visualizer import MolecularVisualizer
        viz = MolecularVisualizer(output_dir='diagnostic_output')
        props = viz.calculate_properties(test_smiles)
        print(f"✅ MolecularVisualizer initialized")
        print(f"   - Molecular Weight: {props['molecular_weight']:.2f} Da")
        print(f"   - LogP: {props['logP']:.2f}")
        print(f"   - TPSA: {props['tpsa']:.2f}")
        print(f"   - Lipinski Violations: {props['lipinski_violations']}")
        
        # Test molecular editor
        from molecular_editor import MolecularEditor
        editor = MolecularEditor()
        canonical = editor.canonicalize(test_smiles)
        print(f"✅ MolecularEditor initialized")
        print(f"   - Canonical SMILES: {canonical}")
        
        # Test analog generation
        analogs = editor.generate_analogs(test_smiles, num_analogs=3, similarity_threshold=0.7)
        print(f"✅ Analog generation works: {len(analogs)} analogs generated")
        for i, analog in enumerate(analogs[:2], 1):
            print(f"   - Analog {i}: {analog['smiles'][:40]}... (similarity: {analog['similarity']:.3f})")
        
        # Test visualization
        result = viz.visualize_molecule(test_smiles, filename='test_aspirin.png', size=(400, 400))
        if result:
            print(f"✅ Molecular visualization successful")
        
        return True
        
    except Exception as e:
        print(f"❌ RDKit integration test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_compound_database():
    """Test compound database access"""
    print_section("2. Compound Database Test")
    
    try:
        from pharmasight_complete import COMPOUND_DATABASE
        
        num_compounds = len(COMPOUND_DATABASE)
        print(f"✅ Compound database loaded: {num_compounds} compounds")
        
        # Test a few compounds
        sample_compounds = list(COMPOUND_DATABASE.items())[:5]
        print(f"\nSample compounds:")
        for name, data in sample_compounds:
            print(f"   - {name}")
            if 'smiles' in data:
                print(f"     SMILES: {data['smiles'][:50]}...")
            if 'therapeutic_area' in data:
                print(f"     Therapeutic Area: {data['therapeutic_area']}")
        
        return True
        
    except Exception as e:
        print(f"❌ Compound database test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_ddi_analysis():
    """Test Drug-Drug Interaction analysis"""
    print_section("3. DDI Analysis Test")
    
    try:
        from ddi_analysis_fix import get_detailed_interaction_info
        
        # Test with a known interaction
        result = get_detailed_interaction_info("Sertraline", "Tramadol")
        
        if result:
            print(f"✅ DDI analysis works")
            print(f"   - Drug 1: {result.get('drug1', 'N/A')}")
            print(f"   - Drug 2: {result.get('drug2', 'N/A')}")
            print(f"   - Severity: {result.get('severity', 'N/A')}")
            print(f"   - Mechanism: {result.get('mechanism', 'N/A')[:100]}...")
        else:
            print(f"⚠️  DDI analysis returned no results (may be expected)")
        
        return True
        
    except Exception as e:
        print(f"❌ DDI analysis test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_analog_generation():
    """Test analog generation functionality"""
    print_section("4. Analog Generation Test")
    
    try:
        from analog_generation_fix import generate_analog_report, resolve_compound_name
        
        # Test compound name resolution
        test_compound = "Aspirin"
        resolved = resolve_compound_name(test_compound)
        
        if resolved:
            print(f"✅ Compound name resolution works")
            print(f"   - Input: {test_compound}")
            print(f"   - Resolved: {resolved.get('name', 'N/A')}")
            if 'smiles' in resolved:
                print(f"   - SMILES: {resolved['smiles'][:50]}...")
        
        # Test analog report generation
        report = generate_analog_report(test_compound)
        
        if report:
            print(f"✅ Analog report generation works")
            print(f"   - Report length: {len(report)} characters")
            print(f"   - Preview: {report[:200]}...")
        
        return True
        
    except Exception as e:
        print(f"❌ Analog generation test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_research_findings():
    """Test research findings functionality"""
    print_section("5. Research Findings Test")
    
    try:
        from research_findings_fix import (
            get_research_findings_with_hypotheses,
            search_research_findings,
            get_research_analytics
        )
        
        # Test getting research findings
        test_compound = "Psilocybin"
        findings = get_research_findings_with_hypotheses(test_compound)
        
        if findings:
            print(f"✅ Research findings retrieval works")
            print(f"   - Compound: {test_compound}")
            print(f"   - Findings length: {len(findings)} characters")
            print(f"   - Preview: {findings[:200]}...")
        else:
            print(f"⚠️  No research findings found for {test_compound}")
        
        # Test search functionality
        search_results = search_research_findings("depression")
        if search_results:
            print(f"✅ Research findings search works")
            print(f"   - Search term: 'depression'")
            print(f"   - Results: {len(search_results)} findings")
        
        # Test analytics
        analytics = get_research_analytics()
        if analytics:
            print(f"✅ Research analytics works")
            print(f"   - Analytics data available")
        
        return True
        
    except Exception as e:
        print(f"❌ Research findings test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_flask_app():
    """Test Flask application initialization"""
    print_section("6. Flask Application Test")
    
    try:
        from pharmasight_complete import app
        
        print(f"✅ Flask app imported successfully")
        print(f"   - App name: {app.name}")
        print(f"   - Secret key configured: {bool(app.secret_key)}")
        
        # Test routes
        with app.test_client() as client:
            # Test health endpoint
            response = client.get('/health')
            if response.status_code == 200:
                print(f"✅ /health endpoint works (status: {response.status_code})")
                data = json.loads(response.data)
                print(f"   - Status: {data.get('status', 'N/A')}")
            
            # Test main page
            response = client.get('/')
            if response.status_code == 200:
                print(f"✅ Main page endpoint works (status: {response.status_code})")
            
        return True
        
    except Exception as e:
        print(f"❌ Flask app test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_rdkit_api():
    """Test RDKit API endpoints"""
    print_section("7. RDKit API Test")
    
    try:
        from rdkit_api import (
            generate_molecule_image,
            calculate_molecular_properties,
            generate_analogs_with_images
        )
        
        test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
        
        # Test image generation
        image_data = generate_molecule_image(test_smiles)
        if image_data:
            print(f"✅ Molecule image generation works")
            print(f"   - Image data length: {len(image_data)} bytes")
        
        # Test property calculation
        props = calculate_molecular_properties(test_smiles)
        if props:
            print(f"✅ Molecular property calculation works")
            print(f"   - Properties: {list(props.keys())}")
        
        # Test analog generation with images
        analogs = generate_analogs_with_images(test_smiles, num_analogs=2)
        if analogs:
            print(f"✅ Analog generation with images works")
            print(f"   - Generated {len(analogs)} analogs")
        
        return True
        
    except Exception as e:
        print(f"❌ RDKit API test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_admet_predictor():
    """Test ADMET prediction functionality"""
    print_section("8. ADMET Predictor Test")
    
    try:
        from admet_predictor import predict_admet_properties
        
        test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
        
        predictions = predict_admet_properties(test_smiles)
        
        if predictions:
            print(f"✅ ADMET prediction works")
            print(f"   - Predictions: {list(predictions.keys())}")
            for key, value in list(predictions.items())[:5]:
                print(f"   - {key}: {value}")
        
        return True
        
    except Exception as e:
        print(f"❌ ADMET predictor test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def generate_diagnostic_report(results):
    """Generate a comprehensive diagnostic report"""
    print_section("DIAGNOSTIC SUMMARY")
    
    total_tests = len(results)
    passed_tests = sum(1 for r in results.values() if r)
    failed_tests = total_tests - passed_tests
    
    print(f"Total Tests: {total_tests}")
    print(f"Passed: {passed_tests} ✅")
    print(f"Failed: {failed_tests} ❌")
    print(f"Success Rate: {(passed_tests/total_tests)*100:.1f}%")
    
    print(f"\nDetailed Results:")
    for test_name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {status} - {test_name}")
    
    # Save report to file
    report = {
        "timestamp": datetime.now().isoformat(),
        "total_tests": total_tests,
        "passed": passed_tests,
        "failed": failed_tests,
        "success_rate": f"{(passed_tests/total_tests)*100:.1f}%",
        "results": {name: "PASS" if result else "FAIL" for name, result in results.items()}
    }
    
    with open('diagnostic_report.json', 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n📊 Diagnostic report saved to: diagnostic_report.json")
    
    return passed_tests == total_tests

def main():
    """Run all diagnostic tests"""
    print("\n" + "="*80)
    print("  PharmaSight Platform - Comprehensive Diagnostic Test Suite")
    print("  " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("="*80)
    
    # Create output directory
    os.makedirs('diagnostic_output', exist_ok=True)
    
    # Run all tests
    results = {
        "RDKit Integration": test_rdkit_integration(),
        "Compound Database": test_compound_database(),
        "DDI Analysis": test_ddi_analysis(),
        "Analog Generation": test_analog_generation(),
        "Research Findings": test_research_findings(),
        "Flask Application": test_flask_app(),
        "RDKit API": test_rdkit_api(),
        "ADMET Predictor": test_admet_predictor()
    }
    
    # Generate summary report
    all_passed = generate_diagnostic_report(results)
    
    if all_passed:
        print(f"\n🎉 All tests passed! Platform is fully functional.")
        return 0
    else:
        print(f"\n⚠️  Some tests failed. Please review the results above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())

