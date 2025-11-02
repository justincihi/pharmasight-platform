#!/usr/bin/env python3
"""
Test script for PharmaSight Microservices Architecture
Tests all services through the API Gateway
"""

import requests
import json
import sys
from typing import Dict, Any

# Configuration
GATEWAY_URL = "http://localhost:8080"
AUTH_URL = f"{GATEWAY_URL}/auth-service"
COMPOUND_URL = f"{GATEWAY_URL}/compound-analysis"
ANALOG_URL = f"{GATEWAY_URL}/analog-generation"
ML_URL = f"{GATEWAY_URL}/ml-models"
QUANTUM_URL = f"{GATEWAY_URL}/quantum-calculator"

# Test credentials
TEST_USER = "researcher"
TEST_PASSWORD = "password"

# ANSI color codes
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
RESET = "\033[0m"

def print_test(name: str):
    """Print test name"""
    print(f"\n{BLUE}Testing: {name}{RESET}")

def print_success(message: str):
    """Print success message"""
    print(f"{GREEN}✓ {message}{RESET}")

def print_error(message: str):
    """Print error message"""
    print(f"{RED}✗ {message}{RESET}")

def print_warning(message: str):
    """Print warning message"""
    print(f"{YELLOW}⚠ {message}{RESET}")

def test_gateway_health() -> bool:
    """Test API Gateway health check"""
    print_test("API Gateway Health")
    try:
        response = requests.get(f"{GATEWAY_URL}/health", timeout=10)
        if response.status_code == 200:
            data = response.json()
            print_success(f"Gateway status: {data.get('gateway')}")
            print_success(f"Overall status: {data.get('overall_status')}")
            
            # Check each service
            services = data.get('services', {})
            all_healthy = True
            for service_name, service_status in services.items():
                status = service_status.get('status')
                if status == 'healthy':
                    print_success(f"  {service_name}: {status}")
                else:
                    print_error(f"  {service_name}: {status}")
                    all_healthy = False
            
            return all_healthy
        else:
            print_error(f"Health check failed: {response.status_code}")
            return False
    except Exception as e:
        print_error(f"Health check error: {e}")
        return False

def test_authentication() -> str:
    """Test authentication service and return access token"""
    print_test("Authentication Service")
    try:
        response = requests.post(
            f"{AUTH_URL}/token",
            data={
                "username": TEST_USER,
                "password": TEST_PASSWORD
            },
            timeout=10
        )
        
        if response.status_code == 200:
            data = response.json()
            token = data.get('access_token')
            print_success(f"Login successful, token received")
            
            # Test getting user info
            headers = {"Authorization": f"Bearer {token}"}
            user_response = requests.get(f"{AUTH_URL}/users/me", headers=headers, timeout=10)
            if user_response.status_code == 200:
                user_data = user_response.json()
                print_success(f"User verified: {user_data.get('username')} ({user_data.get('role')})")
            
            return token
        else:
            print_error(f"Authentication failed: {response.status_code}")
            return None
    except Exception as e:
        print_error(f"Authentication error: {e}")
        return None

def test_compound_analysis(token: str = None) -> bool:
    """Test compound analysis service"""
    print_test("Compound Analysis Service")
    try:
        headers = {"Content-Type": "application/json"}
        if token:
            headers["Authorization"] = f"Bearer {token}"
        
        response = requests.post(
            f"{COMPOUND_URL}/analyze",
            headers=headers,
            json={"compound": "psilocybin", "analysis_type": "full"},
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            print_success("Compound analysis successful")
            
            # Check response structure
            if 'base_properties' in data:
                props = data['base_properties']
                print_success(f"  Name: {props.get('name')}")
                print_success(f"  MW: {props.get('molecular_weight')}")
            
            if 'rdkit_analysis' in data:
                rdkit = data['rdkit_analysis']
                print_success(f"  LogP: {rdkit.get('logp')}")
                print_success(f"  TPSA: {rdkit.get('tpsa')}")
            
            return True
        else:
            print_error(f"Compound analysis failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print_error(f"Compound analysis error: {e}")
        return False

def test_analog_generation(token: str = None) -> bool:
    """Test analog generation service"""
    print_test("Analog Generation Service")
    try:
        headers = {"Content-Type": "application/json"}
        if token:
            headers["Authorization"] = f"Bearer {token}"
        
        response = requests.post(
            f"{ANALOG_URL}/generate",
            headers=headers,
            json={"parent_compound": "psilocybin", "target_properties": "all"},
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            print_success("Analog generation successful")
            print_success(f"  Generated {data.get('count')} analogs")
            
            analogs = data.get('analogs', [])
            for i, analog in enumerate(analogs[:3], 1):
                print_success(f"  {i}. {analog.get('name')} (similarity: {analog.get('similarity'):.2f})")
            
            return True
        else:
            print_error(f"Analog generation failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print_error(f"Analog generation error: {e}")
        return False

def test_ml_predictions(token: str = None) -> bool:
    """Test ML models service"""
    print_test("ML Models Service")
    try:
        headers = {"Content-Type": "application/json"}
        if token:
            headers["Authorization"] = f"Bearer {token}"
        
        # Test ADMET prediction
        response = requests.post(
            f"{ML_URL}/predict/admet",
            headers=headers,
            json={"smiles_list": ["CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]},
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            print_success("ADMET prediction successful")
            
            if 'warning' in data:
                print_warning(f"  {data['warning']}")
            
            predictions = data.get('predictions', {})
            for smiles, pred in list(predictions.items())[:2]:
                print_success(f"  {smiles[:20]}... - Absorption: {pred.get('absorption', 'N/A')}, Toxicity: {pred.get('toxicity', 'N/A')}")
            
            return True
        else:
            print_error(f"ML prediction failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print_error(f"ML prediction error: {e}")
        return False

def test_quantum_calculator(token: str = None) -> bool:
    """Test quantum calculator service"""
    print_test("Quantum Calculator Service")
    try:
        headers = {"Content-Type": "application/json"}
        if token:
            headers["Authorization"] = f"Bearer {token}"
        
        # Test with water molecule
        response = requests.post(
            f"{QUANTUM_URL}/calculate",
            headers=headers,
            json={
                "mol_geometry": "O 0 0 0; H 0 1 0; H 0 0 1",
                "basis": "sto-3g",
                "xc_functional": "b3lyp"
            },
            timeout=60
        )
        
        if response.status_code == 200:
            data = response.json()
            print_success("Quantum calculation successful")
            print_success(f"  Total energy: {data.get('total_energy_hartree', 'N/A')} Hartree")
            print_success(f"  HOMO-LUMO gap: {data.get('homo_lumo_gap_ev', 'N/A')} eV")
            return True
        else:
            print_error(f"Quantum calculation failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print_error(f"Quantum calculation error: {e}")
        return False

def main():
    """Run all tests"""
    print(f"\n{BLUE}{'='*60}")
    print("PharmaSight Microservices Test Suite")
    print(f"{'='*60}{RESET}\n")
    
    results = {}
    
    # Test 1: Gateway health
    results['gateway_health'] = test_gateway_health()
    
    # Test 2: Authentication
    token = test_authentication()
    results['authentication'] = token is not None
    
    # Test 3: Compound analysis
    results['compound_analysis'] = test_compound_analysis(token)
    
    # Test 4: Analog generation
    results['analog_generation'] = test_analog_generation(token)
    
    # Test 5: ML predictions
    results['ml_predictions'] = test_ml_predictions(token)
    
    # Test 6: Quantum calculator (may be slow)
    print_warning("Note: Quantum calculations may take 30-60 seconds")
    results['quantum_calculator'] = test_quantum_calculator(token)
    
    # Summary
    print(f"\n{BLUE}{'='*60}")
    print("Test Summary")
    print(f"{'='*60}{RESET}\n")
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for test_name, result in results.items():
        status = f"{GREEN}PASS{RESET}" if result else f"{RED}FAIL{RESET}"
        print(f"{test_name.replace('_', ' ').title()}: {status}")
    
    print(f"\n{BLUE}Results: {passed}/{total} tests passed{RESET}")
    
    if passed == total:
        print(f"{GREEN}All tests passed! ✓{RESET}\n")
        return 0
    else:
        print(f"{RED}Some tests failed. Please check the output above.{RESET}\n")
        return 1

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(f"\n{YELLOW}Tests interrupted by user{RESET}")
        sys.exit(1)
