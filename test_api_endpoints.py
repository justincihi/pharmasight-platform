#!/usr/bin/env python3
"""Test script to verify all API endpoints are working"""

import requests
import json

BASE_URL = "http://localhost:5000"

def test_endpoint(name, method, endpoint, data=None):
    """Test a single endpoint"""
    print(f"\nTesting {name}...")
    try:
        if method == "POST":
            response = requests.post(f"{BASE_URL}{endpoint}", json=data)
        else:
            response = requests.get(f"{BASE_URL}{endpoint}")
        
        if response.status_code == 200:
            print(f"✅ {name} - SUCCESS")
            result = response.json()
            if 'error' not in result:
                print(f"   Response keys: {list(result.keys())[:5]}...")
            else:
                print(f"   Error: {result['error']}")
        else:
            print(f"❌ {name} - FAILED (Status: {response.status_code})")
            print(f"   Response: {response.text[:200]}")
    except Exception as e:
        print(f"❌ {name} - ERROR: {str(e)}")

# Test data
test_smiles = "CC(C)c1ccc(cc1)C(C)C(O)=O"  # Ibuprofen
test_compounds = [
    {"smiles": "CC(C)c1ccc(cc1)C(C)C(O)=O", "activity": 0.8},
    {"smiles": "CC(C)Cc1ccc(cc1)C(C)C(O)=O", "activity": 0.6},
    {"smiles": "CCCc1ccc(cc1)C(C)C(O)=O", "activity": 0.7}
]

# Test all new endpoints
print("=" * 60)
print("Testing PharmaSight API Endpoints")
print("=" * 60)

# 1. Virtual Screening - Single Compound
test_endpoint(
    "Virtual Screening - Single Compound",
    "POST",
    "/api/vhts/screen_compound",
    {"smiles": test_smiles, "compound_name": "Ibuprofen"}
)

# 2. Virtual Screening - Batch
test_endpoint(
    "Virtual Screening - Batch",
    "POST",
    "/api/vhts/batch_screen",
    {"compounds": test_compounds[:2]}
)

# 3. Lead Optimization
test_endpoint(
    "Lead Optimization",
    "POST",
    "/api/lead_opt/optimize",
    {"smiles": test_smiles, "target_profile": {"logp": 3.5, "mw": 300}}
)

# 4. Off-Target Prediction
test_endpoint(
    "Off-Target Prediction",
    "POST",
    "/api/off_target/predict",
    {"smiles": test_smiles, "primary_target": "COX-2"}
)

# 5. SAR Analysis
test_endpoint(
    "SAR Analysis",
    "POST",
    "/api/sar/analyze",
    {"compound_series": test_compounds}
)

# Also test some existing endpoints to ensure they still work
print("\n" + "=" * 60)
print("Testing Existing Endpoints")
print("=" * 60)

# 6. Analyze Compound
test_endpoint(
    "Analyze Compound",
    "POST",
    "/api/analyze_compound",
    {"smiles": test_smiles}
)

# 7. Health Check
test_endpoint(
    "Health Check",
    "GET",
    "/health"
)

print("\n" + "=" * 60)
print("API Endpoint Testing Complete!")
print("=" * 60)