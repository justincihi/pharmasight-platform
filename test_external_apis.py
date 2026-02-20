#!/usr/bin/env python3
"""
Test script for external API integrations
"""

import sys
import os
sys.path.append('src')

from external_database_apis import PubChemAPI, ChEMBLAPI, FDAOrangeBookAPI

def test_pubchem():
    """Test PubChem API"""
    print("Testing PubChem API...")
    api = PubChemAPI()
    
    # Test compound search
    result = api.search_compound("aspirin")
    if result:
        print(f"✅ PubChem: Found {result.get('name', 'Unknown')} (CID: {result.get('cid')})")
        return True
    else:
        print("❌ PubChem: Search failed")
        return False

def test_chembl():
    """Test ChEMBL API"""
    print("Testing ChEMBL API...")
    api = ChEMBLAPI()
    
    # Test compound search
    result = api.search_compound("aspirin")
    if result:
        print(f"✅ ChEMBL: Found {result.get('name', 'Unknown')} (ID: {result.get('chembl_id')})")
        return True
    else:
        print("❌ ChEMBL: Search failed")
        return False

def test_fda():
    """Test FDA Orange Book API"""
    print("Testing FDA Orange Book API...")
    api = FDAOrangeBookAPI()
    
    # Test drug search
    result = api.search_drug("Advil")
    if result:
        print(f"✅ FDA: Found {result.get('brand_name', 'Unknown')} - {result.get('generic_name')}")
        return True
    else:
        print("❌ FDA: Search failed")
        return False

if __name__ == "__main__":
    print("PharmaSight™ External API Test Suite")
    print("=" * 50)
    
    tests = [
        ("PubChem", test_pubchem),
        ("ChEMBL", test_chembl),
        ("FDA Orange Book", test_fda),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"❌ {name}: Exception - {e}")
            results.append((name, False))
        print()
    
    print("Test Results Summary:")
    print("=" * 30)
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{name}: {status}")
    
    print(f"\nOverall: {passed}/{total} APIs functional")