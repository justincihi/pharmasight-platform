#!/usr/bin/env python3
"""
Test molecular visualization and notation features
"""

import requests
import json
import base64

# Test API endpoints
BASE_URL = "http://127.0.0.1:5000"

def test_2d_visualization():
    """Test 2D visualization with multiple notations"""
    print("\n1. Testing 2D Visualization API...")
    
    # Test SMILES input
    data = {
        "notation": "smiles",
        "value": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "size": [400, 400]
    }
    
    response = requests.post(f"{BASE_URL}/api/visualize/2d", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ 2D Visualization: Success")
        print(f"      - SVG generated: {len(result.get('svg', '')) > 100}")
        print(f"      - PNG base64: {result.get('png_base64', '')[:50]}...")
        print(f"      - SMILES: {result.get('notations', {}).get('smiles')}")
        print(f"      - InChI: {result.get('notations', {}).get('inchi')[:50]}...")
        print(f"      - InChI Key: {result.get('notations', {}).get('inchi_key')}")
        print(f"      - Molecular Weight: {result.get('properties', {}).get('molecular_weight')}")
    else:
        print(f"   ❌ 2D Visualization: Failed ({response.status_code})")

def test_notation_conversion():
    """Test chemical notation conversion"""
    print("\n2. Testing Notation Conversion API...")
    
    data = {
        "from_notation": "smiles",
        "to_notation": "inchi",
        "value": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
    }
    
    response = requests.post(f"{BASE_URL}/api/notation/convert", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ Notation Conversion: Success")
        print(f"      - Requested (InChI): {result.get('requested')[:60]}...")
        print(f"      - SMILES: {result.get('all_formats', {}).get('smiles')}")
        print(f"      - InChI Key: {result.get('all_formats', {}).get('inchi_key')}")
        print(f"      - Molecular Formula: {result.get('all_formats', {}).get('molecular_formula')}")
    else:
        print(f"   ❌ Notation Conversion: Failed ({response.status_code})")

def test_molecule_comparison():
    """Test molecule comparison visualization"""
    print("\n3. Testing Molecule Comparison API...")
    
    data = {
        "smiles1": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "smiles2": "CC(=O)Oc1ccccc1C(=O)OC",  # Aspirin methyl ester
        "labels": ["Aspirin", "Aspirin Methyl Ester"]
    }
    
    response = requests.post(f"{BASE_URL}/api/visualize/compare", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ Molecule Comparison: Success")
        print(f"      - Comparison image generated: {result.get('comparison_image', '')[:50]}...")
        print(f"      - Similarity: {result.get('similarity')}")
        print(f"      - Differences: {result.get('differences')}")
        print(f"      - Molecule 1 atoms: {result.get('molecule1', {}).get('num_atoms')}")
        print(f"      - Molecule 2 atoms: {result.get('molecule2', {}).get('num_atoms')}")
    else:
        print(f"   ❌ Molecule Comparison: Failed ({response.status_code})")

def test_quantum_optimization():
    """Test quantum optimization with visualization"""
    print("\n4. Testing Quantum Optimization API...")
    
    data = {
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "type": "lead_optimization"
    }
    
    response = requests.post(f"{BASE_URL}/api/quantum/optimize", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ Quantum Optimization: Success")
        
        if 'optimized_structures' in result:
            for i, struct in enumerate(result['optimized_structures'][:2]):
                print(f"\n      Optimized Structure {i+1}:")
                print(f"      - Modification: {struct.get('modification')}")
                
                if 'enhanced_visualization' in struct:
                    viz = struct['enhanced_visualization']
                    print(f"      - Similarity Score: {viz.get('similarity_score')}")
                    print(f"      - Image generated: {viz.get('comparison_image', '')[:50]}...")
                    
                    changes = viz.get('structural_changes', {})
                    print(f"      - Atoms changed: {changes.get('atoms_added', 0)}")
                    print(f"      - Bonds changed: {changes.get('bonds_changed', 0)}")
                    
                    notations = viz.get('notations', {})
                    if 'optimized' in notations:
                        print(f"      - New SMILES: {notations['optimized'].get('smiles')}")
                        print(f"      - New InChI Key: {notations['optimized'].get('inchi_key')}")
    else:
        print(f"   ❌ Quantum Optimization: Failed ({response.status_code})")
        print(f"      Error: {response.json()}")

def test_lead_optimization():
    """Test enhanced lead optimization with visualization"""
    print("\n5. Testing Enhanced Lead Optimization API...")
    
    data = {
        "smiles": "CCc1cnn(c1)C2CCNCC2",  # Example compound
        "target_profile": {
            "binding_affinity": "increase",
            "selectivity": "increase"
        },
        "include_visualization": True
    }
    
    response = requests.post(f"{BASE_URL}/api/lead_opt/optimize", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ Lead Optimization: Success")
        
        if 'modified_structures' in result:
            for i, struct in enumerate(result['modified_structures'][:2]):
                print(f"\n      Modified Structure {i+1}:")
                print(f"      - SMILES: {struct.get('smiles')}")
                print(f"      - Score: {struct.get('score')}")
                
                if 'visualization' in struct:
                    viz = struct['visualization']
                    print(f"      - Similarity: {viz.get('similarity')}")
                    print(f"      - Image: {viz.get('comparison_image', '')[:50]}...")
                    
                    changes = viz.get('chemical_changes', {})
                    print(f"      - Atom diff: {changes.get('atom_diff', 0)}")
                    print(f"      - Bond diff: {changes.get('bond_diff', 0)}")
                    print(f"      - Modification type: {changes.get('modification_type')}")
                    
                    if 'notations' in viz:
                        print(f"      - InChI Key: {viz['notations'].get('inchi_key')}")
    else:
        print(f"   ❌ Lead Optimization: Failed ({response.status_code})")

def test_3d_visualization():
    """Test 3D visualization data generation"""
    print("\n6. Testing 3D Visualization API...")
    
    data = {
        "smiles": "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F"  # Celecoxib
    }
    
    response = requests.post(f"{BASE_URL}/api/visualize/3d", json=data)
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ 3D Visualization: Success")
        print(f"      - SDF block generated: {len(result.get('sdf', '')) > 100}")
        print(f"      - Number of atoms: {result.get('num_atoms')}")
        print(f"      - Number of bonds: {result.get('num_bonds')}")
    else:
        print(f"   ❌ 3D Visualization: Failed ({response.status_code})")

if __name__ == "__main__":
    print("=" * 60)
    print("Testing Molecular Visualization & Notation Features")
    print("=" * 60)
    
    test_2d_visualization()
    test_notation_conversion()
    test_molecule_comparison()
    test_quantum_optimization()
    test_lead_optimization()
    test_3d_visualization()
    
    print("\n" + "=" * 60)
    print("Visualization Testing Complete!")
    print("=" * 60)