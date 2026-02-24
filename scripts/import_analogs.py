#!/usr/bin/env python3
"""
Import analogs from PharmaSight Platform ANALOG_GENERATION_DATABASE
Converts Python dict to JSON and imports via Platform API
"""

import sys
import os
import json
import requests

# Add pharmasight-platform src to path
sys.path.insert(0, '/home/ubuntu/pharmasight-platform-latest/src')

from analog_generation_fix import ANALOG_GENERATION_DATABASE

# Dashboard API configuration
API_BASE_URL = "http://localhost:3000"
API_KEY = os.getenv("PLATFORM_API_KEY", "temp-dev-key-12345")

def convert_analog_to_discovery(parent_compound: str, analog: dict, index: int) -> dict:
    """Convert analog dict to discovery format for API"""
    compound_id = f"{parent_compound.upper()}-{analog['name'].replace(' ', '-')}-A{index:03d}"
    
    return {
        "compoundId": compound_id,
        "compoundName": analog["name"],
        "smiles": analog["smiles"],
        "parentCompound": parent_compound.title(),
        "mechanismOfAction": f"Similar mechanism to {parent_compound}",
        "keyDifferences": f"Structural analog with {analog['similarity']*100:.0f}% similarity",
        "confidence": int(analog.get("novelty_score", 85)),
        "similarity": int(analog["similarity"] * 100),
        "safetyScore": analog.get("safety_score", 80),
        "efficacyScore": analog.get("efficacy_score", 80),
        "drugLikenessScore": analog.get("drug_likeness", 85),
        "patentStatus": analog["patent_status"].lower().replace(" ", "-") if "patent" in analog["patent_status"].lower() else "patent-free",
        "marketValue": analog.get("estimated_value", "$0"),
        "discoveryMethod": "analog-generation-database"
    }

def import_all_analogs():
    """Import all analogs from ANALOG_GENERATION_DATABASE"""
    all_discoveries = []
    
    for parent_compound, data in ANALOG_GENERATION_DATABASE.items():
        for idx, analog in enumerate(data["analogs"], 1):
            discovery = convert_analog_to_discovery(parent_compound, analog, idx)
            all_discoveries.append(discovery)
    
    print(f"Prepared {len(all_discoveries)} analogs for import")
    
    # Import via Platform API
    url = f"{API_BASE_URL}/api/platform/discoveries/import"
    payload = {
        "apiKey": API_KEY,
        "discoveries": all_discoveries
    }
    
    try:
        response = requests.post(url, json=payload, timeout=30)
        result = response.json()
        
        print(f"\n✅ Import complete:")
        print(f"  - Imported: {result.get('imported', 0)}")
        print(f"  - Skipped: {result.get('skipped', 0)}")
        print(f"  - Errors: {result.get('errors', 0)}")
        
        if result.get('errors', 0) > 0:
            print(f"\nError details:")
            for error in result.get('details', {}).get('errors', []):
                print(f"  - {error}")
        
        return result
    
    except requests.exceptions.ConnectionError:
        print("❌ Error: Could not connect to dashboard API")
        print("   Make sure the development server is running on http://localhost:3000")
        return None
    except Exception as e:
        print(f"❌ Error: {e}")
        return None

def export_to_json(output_file: str = "master_analogs.json"):
    """Export all analogs to JSON file"""
    all_analogs = []
    
    for parent_compound, data in ANALOG_GENERATION_DATABASE.items():
        for idx, analog in enumerate(data["analogs"], 1):
            discovery = convert_analog_to_discovery(parent_compound, analog, idx)
            all_analogs.append(discovery)
    
    with open(output_file, 'w') as f:
        json.dump(all_analogs, f, indent=2)
    
    print(f"✅ Exported {len(all_analogs)} analogs to {output_file}")
    return all_analogs

if __name__ == "__main__":
    print("PharmaSight™ Analog Import Tool")
    print("=" * 50)
    
    # Export to JSON first
    output_path = "/home/ubuntu/pharmasight-admin-dashboard/master_analogs.json"
    analogs = export_to_json(output_path)
    
    # Import to database
    print("\nImporting to database...")
    import_all_analogs()
