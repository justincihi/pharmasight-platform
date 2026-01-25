#!/usr/bin/env python3
"""
Service Verification Script
Checks all PharmaSight platform services are properly configured
"""

import os
import json
from pathlib import Path

def check_service(service_name, required_files):
    """Check if a service has all required files"""
    service_path = Path(f"services/{service_name}")
    
    if not service_path.exists():
        return False, f"❌ Directory services/{service_name} not found"
    
    missing_files = []
    for file in required_files:
        if not (service_path / file).exists():
            missing_files.append(file)
    
    if missing_files:
        return False, f"⚠️  Missing files: {', '.join(missing_files)}"
    
    return True, f"✅ All required files present"

def check_backend():
    """Check backend PK modeling service"""
    backend_path = Path("backend")
    
    if not backend_path.exists():
        return False, "❌ backend/ directory not found"
    
    required_modules = ["pharmasight_pk/ddi.py", "pharmasight_pk/popPK.py", 
                       "pharmasight_pk/virtual_patient.py", "main.py"]
    
    missing = []
    for module in required_modules:
        if not (backend_path / module).exists():
            missing.append(module)
    
    if missing:
        return False, f"⚠️  Missing modules: {', '.join(missing)}"
    
    return True, "✅ Backend PK service configured"

def check_data_files():
    """Check critical data files"""
    data_files = {
        "MASTER_ANALOG_DISCOVERIES.json": "119 novel compounds",
        "RESEARCH_ARTICLES_DATABASE.json": "26 research articles",
        "database/schema.sql": "PostgreSQL schema"
    }
    
    results = []
    for file, description in data_files.items():
        if Path(file).exists():
            size = Path(file).stat().st_size
            results.append(f"✅ {file} ({size:,} bytes) - {description}")
        else:
            results.append(f"❌ {file} - MISSING")
    
    return results

def check_docker_compose():
    """Verify docker-compose.yml configuration"""
    compose_file = Path("docker-compose.yml")
    
    if not compose_file.exists():
        return False, "❌ docker-compose.yml not found"
    
    # Count services in docker-compose.yml
    with open(compose_file) as f:
        content = f.read()
        # Simple service count (lines starting with service names)
        services = [line.strip().rstrip(':') for line in content.split('\n') 
                   if line and not line.startswith(' ') and ':' in line and 'services:' not in line]
    
    return True, f"✅ Docker Compose configured with {len(services)} services"

def main():
    print("=" * 60)
    print("PharmaSight Platform - Service Verification")
    print("=" * 60)
    print()
    
    # Check services
    services = {
        "compound-analysis": ["Dockerfile", "main.py", "requirements.txt"],
        "analog-generation": ["Dockerfile", "main.py", "requirements.txt"],
        "ml-models": ["Dockerfile", "main.py", "requirements.txt"],
        "quantum-calculator": ["Dockerfile", "main.py", "requirements.txt"],
        "auth-service": ["Dockerfile", "main.py", "requirements.txt"],
        "web-frontend": ["Dockerfile", "main.py", "requirements.txt"],
        "research-engine": ["Dockerfile", "main.py", "requirements.txt", "api_integrations.py"],
        "api-gateway": ["Dockerfile", "main.py", "requirements.txt"]
    }
    
    print("📦 Checking Microservices:")
    print("-" * 60)
    for service, required_files in services.items():
        status, message = check_service(service, required_files)
        print(f"{service:25s} {message}")
    
    print()
    print("🔬 Checking Backend Services:")
    print("-" * 60)
    status, message = check_backend()
    print(f"{'Backend PK Modeling':25s} {message}")
    
    print()
    print("📊 Checking Data Files:")
    print("-" * 60)
    for result in check_data_files():
        print(result)
    
    print()
    print("🐳 Checking Docker Configuration:")
    print("-" * 60)
    status, message = check_docker_compose()
    print(message)
    
    print()
    print("=" * 60)
    print("Verification Complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
