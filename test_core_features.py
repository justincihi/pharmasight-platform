#!/usr/bin/env python3
"""
Test script for core PharmaSight features
Tests compound analysis, analog generation, DDI analysis, etc.
"""

import sys
import json
import requests
import time

# Test configuration
BASE_URL = "http://127.0.0.1:5008"
TEST_COMPOUNDS = ["psilocybin", "ketamine", "sertraline", "alprazolam"]

def test_compound_analysis():
    """Test compound analysis feature"""
    print("Testing Compound Analysis...")
    success_count = 0
    
    for compound in TEST_COMPOUNDS:
        try:
            response = requests.post(
                f"{BASE_URL}/api/analyze_compound",
                json={"compound": compound},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                if 'error' not in data:
                    print(f"  ✅ {compound}: {data.get('name')} - MW: {data.get('molecular_weight')}g/mol")
                    success_count += 1
                else:
                    print(f"  ❌ {compound}: {data['error']}")
            else:
                print(f"  ❌ {compound}: HTTP {response.status_code}")
                
        except Exception as e:
            print(f"  ❌ {compound}: Exception - {e}")
    
    return success_count, len(TEST_COMPOUNDS)

def test_analog_generation():
    """Test analog generation feature"""
    print("Testing Analog Generation...")
    success_count = 0
    test_compounds = ["psilocybin", "ketamine", "mdma"]
    
    for compound in test_compounds:
        try:
            response = requests.post(
                f"{BASE_URL}/api/generate_analogs",
                json={
                    "parent_compound": compound,
                    "target_properties": "patent-free"
                },
                timeout=15
            )
            
            if response.status_code == 200:
                data = response.json()
                if 'error' not in data:
                    analog_count = len(data.get('analogs', []))
                    print(f"  ✅ {compound}: {analog_count} analogs generated")
                    success_count += 1
                else:
                    print(f"  ❌ {compound}: {data['error']}")
            else:
                print(f"  ❌ {compound}: HTTP {response.status_code}")
                
        except Exception as e:
            print(f"  ❌ {compound}: Exception - {e}")
    
    return success_count, len(test_compounds)

def test_ddi_analysis():
    """Test drug-drug interaction analysis"""
    print("Testing DDI Analysis...")
    
    test_cases = [
        (["sertraline", "tramadol"], "High risk interaction"),
        (["alprazolam", "oxycodone"], "CNS depression"),
        (["ketamine", "alprazolam"], "Moderate interaction")
    ]
    
    success_count = 0
    
    for medications, expected_desc in test_cases:
        try:
            response = requests.post(
                f"{BASE_URL}/api/analyze_pkpd",
                json={
                    "medications": medications,
                    "patient_data": {
                        "age_group": "31-50",
                        "liver_disease": False,
                        "kidney_disease": False
                    }
                },
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                interactions = data.get('interactions', [])
                safety_score = data.get('safety_score', 0)
                
                print(f"  ✅ {' + '.join(medications)}: Safety Score {safety_score}%, {len(interactions)} interactions")
                success_count += 1
            else:
                print(f"  ❌ {' + '.join(medications)}: HTTP {response.status_code}")
                
        except Exception as e:
            print(f"  ❌ {' + '.join(medications)}: Exception - {e}")
    
    return success_count, len(test_cases)

def test_research_findings():
    """Test research findings feature"""
    print("Testing Research Findings...")
    
    try:
        response = requests.get(f"{BASE_URL}/api/research_findings", timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            findings = data.get('findings', [])
            analytics = data.get('analytics', {})
            
            print(f"  ✅ Research Findings: {len(findings)} findings loaded")
            if analytics:
                print(f"      - Average Confidence: {analytics.get('average_confidence', 'N/A')}%")
            return 1, 1
        else:
            print(f"  ❌ Research Findings: HTTP {response.status_code}")
            return 0, 1
            
    except Exception as e:
        print(f"  ❌ Research Findings: Exception - {e}")
        return 0, 1

def test_receptor_profiling():
    """Test receptor profiling feature"""
    print("Testing Receptor Profiling...")
    success_count = 0
    test_compounds = ["psilocybin", "ketamine"]
    
    for compound in test_compounds:
        try:
            response = requests.get(
                f"{BASE_URL}/api/receptor_profile/{compound}",
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                if data.get('success'):
                    profile = data.get('profile', {})
                    interactions = profile.get('interactions', [])
                    print(f"  ✅ {compound}: {len(interactions)} receptor interactions")
                    success_count += 1
                else:
                    print(f"  ❌ {compound}: {data.get('error', 'Unknown error')}")
            else:
                print(f"  ❌ {compound}: HTTP {response.status_code}")
                
        except Exception as e:
            print(f"  ❌ {compound}: Exception - {e}")
    
    return success_count, len(test_compounds)

def test_health_check():
    """Test application health"""
    print("Testing Application Health...")
    
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            print(f"  ✅ Health Check: {data.get('status')} - Version {data.get('version')}")
            return True
        else:
            print(f"  ❌ Health Check: HTTP {response.status_code}")
            return False
            
    except Exception as e:
        print(f"  ❌ Health Check: Exception - {e}")
        return False

def main():
    print("PharmaSight™ Core Feature Test Suite")
    print("=" * 50)
    
    # First check if the application is running
    if not test_health_check():
        print("\n❌ Application not accessible. Please ensure it's running on http://127.0.0.1:5008")
        return
    
    print()
    
    # Run all feature tests
    tests = [
        ("Compound Analysis", test_compound_analysis),
        ("Analog Generation", test_analog_generation),
        ("DDI Analysis", test_ddi_analysis),
        ("Research Findings", test_research_findings),
        ("Receptor Profiling", test_receptor_profiling),
    ]
    
    total_passed = 0
    total_tests = 0
    
    for test_name, test_func in tests:
        try:
            passed, total = test_func()
            total_passed += passed
            total_tests += total
            print(f"  Result: {passed}/{total} passed\n")
        except Exception as e:
            print(f"  ❌ Test suite error: {e}\n")
    
    print("=" * 50)
    print(f"Overall Results: {total_passed}/{total_tests} tests passed")
    print(f"Success Rate: {(total_passed/total_tests)*100:.1f}%")
    
    if total_passed == total_tests:
        print("🎉 All core features are functional!")
    else:
        print("⚠️  Some features need attention")

if __name__ == "__main__":
    main()