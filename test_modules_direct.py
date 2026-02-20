#!/usr/bin/env python3
"""
Manual test of core functionality without external requests
Tests the core PharmaSight modules directly
"""

import sys
import os
sys.path.append('src')

def test_compound_database():
    """Test compound database directly"""
    print("Testing Compound Database...")
    
    try:
        from pharmasight_complete import COMPOUND_DATABASE
        
        total_compounds = len(COMPOUND_DATABASE)
        print(f"  ✅ Compound Database: {total_compounds} compounds loaded")
        
        # Test a few specific compounds
        test_compounds = ["psilocybin", "ketamine", "sertraline", "alprazolam"]
        found = 0
        
        for compound in test_compounds:
            if compound in COMPOUND_DATABASE:
                data = COMPOUND_DATABASE[compound]
                print(f"    - {data['name']}: MW {data['molecular_weight']}g/mol, {data['therapeutic_area']}")
                found += 1
        
        print(f"  Found {found}/{len(test_compounds)} test compounds")
        return True
        
    except Exception as e:
        print(f"  ❌ Compound Database: {e}")
        return False

def test_analog_generation():
    """Test analog generation module"""
    print("Testing Analog Generation Module...")
    
    try:
        from analog_generation_fix import generate_analog_report
        
        result = generate_analog_report("psilocybin", "patent-free")
        
        if 'error' not in result:
            analogs = result.get('analogs', [])
            print(f"  ✅ Analog Generation: {len(analogs)} analogs for psilocybin")
            
            if analogs:
                print(f"    - Top analog: {analogs[0]['name']} (similarity: {analogs[0]['similarity']:.2f})")
            
            return True
        else:
            print(f"  ❌ Analog Generation: {result['error']}")
            return False
            
    except Exception as e:
        print(f"  ❌ Analog Generation: {e}")
        return False

def test_ddi_analysis():
    """Test drug-drug interaction analysis"""
    print("Testing DDI Analysis Module...")
    
    try:
        from ddi_analysis_fix import get_detailed_interaction_info
        
        interaction = get_detailed_interaction_info("sertraline", "tramadol")
        
        if interaction:
            print(f"  ✅ DDI Analysis: {interaction['risk_level']} risk interaction found")
            print(f"    - Mechanism: {interaction['mechanism'][:50]}...")
            return True
        else:
            print(f"  ❌ DDI Analysis: No interaction data found")
            return False
            
    except Exception as e:
        print(f"  ❌ DDI Analysis: {e}")
        return False

def test_research_findings():
    """Test research findings module"""
    print("Testing Research Findings Module...")
    
    try:
        from research_findings_fix import get_research_findings_with_hypotheses
        
        findings = get_research_findings_with_hypotheses()
        
        if findings:
            print(f"  ✅ Research Findings: {len(findings)} findings loaded")
            
            # Check for different types
            regular_findings = [f for f in findings if not f['id'].startswith('HYP')]
            hypotheses = [f for f in findings if f['id'].startswith('HYP')]
            
            print(f"    - Regular findings: {len(regular_findings)}")
            print(f"    - AI hypotheses: {len(hypotheses)}")
            return True
        else:
            print(f"  ❌ Research Findings: No findings loaded")
            return False
            
    except Exception as e:
        print(f"  ❌ Research Findings: {e}")
        return False

def test_receptor_pharmacology():
    """Test receptor pharmacology module"""
    print("Testing Receptor Pharmacology Module...")
    
    try:
        from receptor_pharmacology import get_receptor_database, ReceptorProfile
        
        db = get_receptor_database()
        receptors = db.get_all_receptors()
        
        total_families = len(receptors)
        total_subtypes = sum(len(family['subtypes']) for family in receptors.values())
        
        print(f"  ✅ Receptor Database: {total_families} families, {total_subtypes} subtypes")
        
        # Test profile creation
        profile = ReceptorProfile("test_compound")
        print(f"    - Profile system operational")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Receptor Pharmacology: {e}")
        return False

def test_external_apis():
    """Test external API modules"""
    print("Testing External API Modules...")
    
    try:
        from external_database_apis import PubChemAPI, ChEMBLAPI, FDAOrangeBookAPI
        
        # Test API class instantiation
        pubchem = PubChemAPI()
        chembl = ChEMBLAPI()
        fda = FDAOrangeBookAPI()
        
        print(f"  ✅ External APIs: All API classes instantiated successfully")
        print(f"    - PubChem API: {pubchem.BASE_URL}")
        print(f"    - ChEMBL API: {chembl.BASE_URL}")
        print(f"    - FDA API: {fda.BASE_URL}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ External APIs: {e}")
        return False

def main():
    print("PharmaSight™ Direct Module Test Suite")
    print("=" * 50)
    
    tests = [
        ("Compound Database", test_compound_database),
        ("Analog Generation", test_analog_generation),
        ("DDI Analysis", test_ddi_analysis),
        ("Research Findings", test_research_findings),
        ("Receptor Pharmacology", test_receptor_pharmacology),
        ("External APIs", test_external_apis),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            success = test_func()
            if success:
                passed += 1
            print()
        except Exception as e:
            print(f"  ❌ {test_name}: Exception during test - {e}")
            print()
    
    print("=" * 50)
    print(f"Module Test Results: {passed}/{total} passed")
    print(f"Success Rate: {(passed/total)*100:.1f}%")
    
    if passed == total:
        print("🎉 All core modules are functional!")
    else:
        print("⚠️  Some modules need attention")

if __name__ == "__main__":
    main()