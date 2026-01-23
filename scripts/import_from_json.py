#!/usr/bin/env python3
"""
Import all analogs from master_analogs.json into the database
"""

import json
import requests
import os

API_BASE_URL = "http://localhost:3000"
API_KEY = os.getenv("PLATFORM_API_KEY", "")
MASTER_FILE = "/home/ubuntu/pharmasight-admin-dashboard/master_analogs.json"

def import_analogs_from_json():
    """Import all analogs from master JSON file"""
    
    if not API_KEY:
        print("❌ Error: PLATFORM_API_KEY environment variable not set")
        return False
    
    # Load analogs from JSON
    with open(MASTER_FILE, 'r') as f:
        analogs = json.load(f)
    
    print(f"Loaded {len(analogs)} analogs from {MASTER_FILE}")
    
    # Import in batches of 50 to avoid timeout
    batch_size = 50
    total_imported = 0
    total_skipped = 0
    total_errors = 0
    
    for i in range(0, len(analogs), batch_size):
        batch = analogs[i:i+batch_size]
        print(f"\nImporting batch {i//batch_size + 1} ({len(batch)} analogs)...")
        
        url = f"{API_BASE_URL}/api/platform/discoveries/import"
        payload = {
            "apiKey": API_KEY,
            "discoveries": batch
        }
        
        try:
            response = requests.post(url, json=payload, timeout=60)
            result = response.json()
            
            total_imported += result.get('imported', 0)
            total_skipped += result.get('skipped', 0)
            total_errors += result.get('errors', 0)
            
            print(f"  ✅ Imported: {result.get('imported', 0)}")
            print(f"  ⏭️  Skipped: {result.get('skipped', 0)}")
            if result.get('errors', 0) > 0:
                print(f"  ❌ Errors: {result.get('errors', 0)}")
        
        except Exception as e:
            print(f"  ❌ Error: {e}")
            total_errors += len(batch)
    
    print(f"\n{'='*50}")
    print(f"Import Summary:")
    print(f"  Total analogs: {len(analogs)}")
    print(f"  ✅ Imported: {total_imported}")
    print(f"  ⏭️  Skipped: {total_skipped}")
    print(f"  ❌ Errors: {total_errors}")
    print(f"{'='*50}")
    
    return total_errors == 0

if __name__ == "__main__":
    print("PharmaSight™ Master Analog Import")
    print("=" * 50)
    success = import_analogs_from_json()
    exit(0 if success else 1)
