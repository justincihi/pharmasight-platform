#!/usr/bin/env python3
"""
Test External Database Connectivity for PharmaSight™ Platform
Tests connections to public chemical information repositories
"""

import sys
sys.path.insert(0, 'src')

from external_database_apis import PubChemAPI, ChEMBLAPI, FDAOrangeBookAPI, UnifiedDatabaseSearch

def test_database_connections():
    """Test connections to external databases"""
    
    print("=" * 60)
    print("Testing External Database Connections")
    print("=" * 60)
    
    # Test compound: Aspirin (well-known drug)
    test_compound = "aspirin"
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin SMILES
    
    # 1. Test PubChem
    print("\n1. Testing PubChem API...")
    try:
        pubchem = PubChemAPI()
        result = pubchem.search_compound(test_compound)
        if result:
            print(f"   ✅ PubChem: Connected - Found {test_compound}")
            print(f"      Molecular Formula: {result.get('molecular_formula')}")
            print(f"      Molecular Weight: {result.get('molecular_weight')}")
        else:
            print(f"   ⚠️  PubChem: No data returned for {test_compound}")
    except Exception as e:
        print(f"   ❌ PubChem: Connection failed - {e}")
    
    # 2. Test ChEMBL
    print("\n2. Testing ChEMBL API...")
    try:
        chembl = ChEMBLAPI()
        result = chembl.search_compound(test_compound)
        if result:
            print(f"   ✅ ChEMBL: Connected - Found {test_compound}")
            print(f"      ChEMBL ID: {result.get('chembl_id')}")
            print(f"      Max Phase: {result.get('max_phase')}")
        else:
            print(f"   ⚠️  ChEMBL: No data returned for {test_compound}")
    except Exception as e:
        print(f"   ❌ ChEMBL: Connection failed - {e}")
    
    # 3. Test FDA Orange Book
    print("\n3. Testing FDA Orange Book API...")
    try:
        fda = FDAOrangeBookAPI()
        result = fda.search_drug(test_compound)
        if result:
            print(f"   ✅ FDA: Connected - Found {test_compound}")
            print(f"      Brand Name: {result.get('brand_name')}")
            print(f"      Generic Name: {result.get('generic_name')}")
        else:
            print(f"   ⚠️  FDA: No data returned for {test_compound}")
    except Exception as e:
        print(f"   ❌ FDA: Connection failed - {e}")
    
    # 4. Test Unified Search
    print("\n4. Testing Unified Database Search...")
    try:
        unified = UnifiedDatabaseSearch()
        results = unified.search_compound_all_databases(test_compound)
        
        active_databases = sum(1 for v in results.values() if v is not None)
        print(f"   📊 Active databases: {active_databases}/4")
        
        for db_name, data in results.items():
            if data:
                print(f"      ✅ {db_name.upper()}: Data retrieved")
            else:
                print(f"      ⚠️  {db_name.upper()}: No data")
                
    except Exception as e:
        print(f"   ❌ Unified Search: Failed - {e}")
    
    # 5. Test Similar Compound Search
    print("\n5. Testing Similar Compound Search...")
    try:
        pubchem = PubChemAPI()
        similar = pubchem.search_similar_compounds(test_smiles, threshold=0.9)
        if similar:
            print(f"   ✅ Similar Compounds: Found {len(similar)} compounds")
            for i, comp in enumerate(similar[:3], 1):
                print(f"      {i}. CID: {comp.get('cid')}")
        else:
            print(f"   ⚠️  Similar Compounds: No results")
    except Exception as e:
        print(f"   ❌ Similar Compounds: Search failed - {e}")
    
    print("\n" + "=" * 60)
    print("Database Connection Test Complete")
    print("=" * 60)

if __name__ == "__main__":
    test_database_connections()