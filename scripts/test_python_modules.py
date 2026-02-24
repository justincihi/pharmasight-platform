#!/usr/bin/env python3
"""
Test Python modules for basic functionality
"""

import sys
import importlib.util
from pathlib import Path

def test_import(module_path, module_name):
    """Test if a Python module can be imported"""
    try:
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None:
            return False, "❌ Failed to load spec"
        
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        
        return True, f"✅ Loaded successfully"
    except Exception as e:
        return False, f"❌ Error: {str(e)[:50]}"

def check_json_validity(file_path):
    """Check if JSON file is valid"""
    try:
        import json
        with open(file_path) as f:
            data = json.load(f)
        return True, f"✅ Valid JSON ({len(data)} items)"
    except Exception as e:
        return False, f"❌ Invalid JSON: {str(e)[:50]}"

def main():
    print("=" * 60)
    print("PharmaSight Platform - Python Module Tests")
    print("=" * 60)
    print()
    
    # Test backend modules
    print("🔬 Testing Backend PK Modules:")
    print("-" * 60)
    
    backend_modules = [
        ("backend/pharmasight_pk/ddi.py", "ddi"),
        ("backend/pharmasight_pk/popPK.py", "popPK"),
        ("backend/pharmasight_pk/virtual_patient.py", "virtual_patient"),
    ]
    
    for path, name in backend_modules:
        if Path(path).exists():
            status, message = test_import(path, name)
            print(f"{name:25s} {message}")
        else:
            print(f"{name:25s} ❌ File not found")
    
    print()
    print("📊 Testing Data Files:")
    print("-" * 60)
    
    data_files = [
        "MASTER_ANALOG_DISCOVERIES.json",
        "RESEARCH_ARTICLES_DATABASE.json"
    ]
    
    for file in data_files:
        if Path(file).exists():
            status, message = check_json_validity(file)
            print(f"{Path(file).name:40s} {message}")
        else:
            print(f"{Path(file).name:40s} ❌ Not found")
    
    print()
    print("=" * 60)
    print("Module Testing Complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
