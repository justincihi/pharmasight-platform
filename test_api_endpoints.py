#!/usr/bin/env python3
"""Test script to verify all API endpoints are working"""

import requests
import json

BASE_URL = "http://localhost:5000"

def test_endpoint(name, method, endpoint, data=None, expected_type=dict):
    """Test a single endpoint with validation"""
    print(f"\nTesting {name}...")
    try:
        if method == "POST":
            response = requests.post(f"{BASE_URL}{endpoint}", json=data)
        else:
            response = requests.get(f"{BASE_URL}{endpoint}")
        
        if response.status_code == 200:
            print(f"✅ {name} - SUCCESS (Status: 200)")
            result = response.json()
            
            # Handle both dict and list responses
            if isinstance(result, list):
                print(f"   Response type: list with {len(result)} items")
                if result:
                    print(f"   First item keys: {list(result[0].keys())[:5]}...")
            elif isinstance(result, dict):
                if 'error' not in result:
                    print(f"   Response keys: {list(result.keys())[:5]}...")
                else:
                    print(f"   Error: {result['error']}")
            
            # Validate response type
            if expected_type == list and not isinstance(result, list):
                print(f"   ⚠️  Warning: Expected list but got {type(result).__name__}")
            elif expected_type == dict and not isinstance(result, dict):
                print(f"   ⚠️  Warning: Expected dict but got {type(result).__name__}")
                
            return True, result
        else:
            print(f"❌ {name} - FAILED (Status: {response.status_code})")
            print(f"   Response: {response.text[:200]}")
            return False, None
    except Exception as e:
        print(f"❌ {name} - ERROR: {str(e)}")
        return False, None

def validate_response_structure(name, response, required_fields):
    """Validate that response contains required fields"""
    print(f"   Validating structure...")
    missing_fields = []
    
    if isinstance(response, dict):
        for field in required_fields:
            if field not in response:
                missing_fields.append(field)
    elif isinstance(response, list) and response:
        for field in required_fields:
            if field not in response[0]:
                missing_fields.append(field)
    
    if missing_fields:
        print(f"   ⚠️  Missing required fields: {missing_fields}")
        return False
    else:
        print(f"   ✅ All required fields present")
        return True

# Test data
test_smiles = "CC(C)c1ccc(cc1)C(C)C(O)=O"  # Ibuprofen
test_compounds = [
    {"smiles": "CC(C)c1ccc(cc1)C(C)C(O)=O", "activity": 0.8},
    {"smiles": "CC(C)Cc1ccc(cc1)C(C)C(O)=O", "activity": 0.6},
    {"smiles": "CCCc1ccc(cc1)C(C)C(O)=O", "activity": 0.7}
]

# Test all new endpoints with validation
print("=" * 60)
print("Testing PharmaSight API Endpoints with Validation")
print("=" * 60)

# 1. Virtual Screening - Single Compound
success, response = test_endpoint(
    "Virtual Screening - Single Compound",
    "POST",
    "/api/vhts/screen_compound",
    {"smiles": test_smiles, "compound_name": "Ibuprofen"}
)
if success:
    validate_response_structure("Virtual Screening", response, 
        ["compound_id", "receptor_hits", "polypharmacology_score"])

# 2. Virtual Screening - Batch
success, response = test_endpoint(
    "Virtual Screening - Batch",
    "POST",
    "/api/vhts/batch_screen",
    {"compounds": test_compounds[:2]},
    expected_type=list  # Expect a list response
)
if success:
    validate_response_structure("Batch Screening", response,
        ["compound_id", "receptor_hits"])

# 3. Lead Optimization
success, response = test_endpoint(
    "Lead Optimization",
    "POST",
    "/api/lead_opt/optimize",
    {"smiles": test_smiles, "target_profile": {"logp": 3.5, "mw": 300}}
)
if success:
    validate_response_structure("Lead Optimization", response,
        ["original_smiles", "optimization_strategies", "modified_structures"])

# 4. Off-Target Prediction
success, response = test_endpoint(
    "Off-Target Prediction",
    "POST",
    "/api/off_target/predict",
    {"smiles": test_smiles, "primary_target": "COX-2"}
)
if success:
    validate_response_structure("Off-Target Prediction", response,
        ["off_target_hits", "risk_score", "safety_classification"])

# 5. SAR Analysis
success, response = test_endpoint(
    "SAR Analysis",
    "POST",
    "/api/sar/analyze",
    {"compound_series": test_compounds}
)
if success:
    validate_response_structure("SAR Analysis", response,
        ["sar_matrix", "activity_cliffs", "key_features"])

# Test validation with invalid inputs
print("\n" + "=" * 60)
print("Testing Input Validation")
print("=" * 60)

# Test with missing SMILES
test_endpoint(
    "Virtual Screening - Missing SMILES",
    "POST",
    "/api/vhts/screen_compound",
    {"compound_name": "No SMILES"}
)

# Test with invalid SMILES
test_endpoint(
    "Lead Optimization - Invalid SMILES",
    "POST",
    "/api/lead_opt/optimize",
    {"smiles": "INVALID_SMILES_STRING"}
)

# Also test some existing endpoints to ensure they still work
print("\n" + "=" * 60)
print("Testing Existing Endpoints")
print("=" * 60)

# 6. Analyze Compound
success, response = test_endpoint(
    "Analyze Compound",
    "POST",
    "/api/analyze_compound",
    {"smiles": test_smiles}
)

# 7. Health Check
success, response = test_endpoint(
    "Health Check",
    "GET",
    "/health"
)

print("\n" + "=" * 60)
print("API Endpoint Testing Complete!")
print("=" * 60)