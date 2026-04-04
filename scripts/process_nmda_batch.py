#!/usr/bin/env python3
"""
Process batch-downloaded NMDA receptor structures from RCSB PDB
Converts PDB to PDBQT format and prepares for docking
"""

import os
import subprocess
import json
from pathlib import Path

# NMDA structure metadata
NMDA_STRUCTURES = {
    '8xlk': {
        'name': 'NMDA GluN1-GluN2A-GluN2B (Tri-heteromeric)',
        'species': 'rat',
        'description': 'Native tri-heteromeric NMDA receptor from rat cortex and hippocampus',
        'resolution': 2.85,
        'year': 2023,
        'ligands': ['glutamate', 'glycine'],
        'pdb_id': '8XLK',
        'family': 'nmda',
        'subtype': 'glun2a-glun2b',
        'modulation_mode': 'ion-channel',
    },
    '9jnn': {
        'name': 'NMDA GluN1-GluN2B (Di-heteromeric)',
        'species': 'rat',
        'description': 'Native di-heteromeric NMDA receptor from rat cortex and hippocampus',
        'resolution': 2.95,
        'year': 2024,
        'ligands': ['glutamate', 'glycine'],
        'pdb_id': '9JNN',
        'family': 'nmda',
        'subtype': 'glun2b',
        'modulation_mode': 'ion-channel',
    },
}

def prepare_structure(pdb_id, input_dir, output_dir):
    """Prepare PDB structure for docking"""
    pdb_file = os.path.join(input_dir, f'{pdb_id}.pdb')
    pdbqt_file = os.path.join(output_dir, f'{pdb_id}.pdbqt')
    
    if not os.path.exists(pdb_file):
        print(f"❌ PDB file not found: {pdb_file}")
        return None
    
    print(f"\n📦 Processing {pdb_id}...")
    
    try:
        # Use meeko to prepare the structure
        cmd = [
            'python3', '-m', 'meeko',
            '--input', pdb_file,
            '--output', pdbqt_file,
            '--hydrate', 'yes',
            '--add_hydrogens',
            '--add_atom_types',
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            print(f"⚠️  Meeko preparation failed, trying obabel...")
            # Fallback to obabel
            cmd = [
                'obabel',
                pdb_file,
                '-O', pdbqt_file,
                '-xh',  # Add hydrogens
                '-xr',  # Add rigid root
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if result.returncode != 0:
                print(f"❌ Failed to prepare {pdb_id}")
                print(f"Error: {result.stderr}")
                return None
        
        # Verify output file
        if os.path.exists(pdbqt_file):
            file_size = os.path.getsize(pdbqt_file) / 1024  # KB
            print(f"✅ {pdb_id} prepared successfully ({file_size:.1f} KB)")
            return pdbqt_file
        else:
            print(f"❌ PDBQT file not created for {pdb_id}")
            return None
            
    except subprocess.TimeoutExpired:
        print(f"❌ Timeout processing {pdb_id}")
        return None
    except Exception as e:
        print(f"❌ Error processing {pdb_id}: {e}")
        return None

def main():
    # Setup directories
    input_dir = '/home/ubuntu/nmda_structures'
    output_dir = '/home/ubuntu/pharmasight-admin-dashboard/data/nmda_receptors'
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("NMDA Receptor Batch Processing")
    print("=" * 60)
    
    results = {}
    
    for pdb_id, metadata in NMDA_STRUCTURES.items():
        pdbqt_file = prepare_structure(pdb_id, input_dir, output_dir)
        
        if pdbqt_file:
            results[pdb_id] = {
                'pdbqt_file': pdbqt_file,
                'metadata': metadata,
                'status': 'success',
            }
        else:
            results[pdb_id] = {
                'metadata': metadata,
                'status': 'failed',
            }
    
    # Save results to JSON
    results_file = os.path.join(output_dir, 'nmda_structures_manifest.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Processing Summary")
    print("=" * 60)
    
    successful = sum(1 for r in results.values() if r['status'] == 'success')
    failed = sum(1 for r in results.values() if r['status'] == 'failed')
    
    print(f"✅ Successful: {successful}/{len(results)}")
    print(f"❌ Failed: {failed}/{len(results)}")
    print(f"\nResults saved to: {results_file}")
    
    # List prepared files
    print("\nPrepared PDBQT files:")
    for pdb_id, result in results.items():
        if result['status'] == 'success':
            print(f"  - {pdb_id}: {result['pdbqt_file']}")

if __name__ == '__main__':
    main()
