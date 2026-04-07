#!/usr/bin/env python3

"""
Download and prepare PDBQT receptor files for psychiatric drug targets
Uses RCSB PDB API to fetch structures and converts them to PDBQT format
"""

import os
import json
import subprocess
import requests
from pathlib import Path

# Psychiatric drug targets with their PDB IDs and descriptions
PSYCHIATRIC_TARGETS = {
    'NMDA': {
        'pdb_id': '6NQH',  # NMDA receptor with GluN1/GluN2B subunits
        'description': 'NMDA receptor (N-methyl-D-aspartate receptor)',
        'chain': 'A',
        'ligand_chain': 'C',
    },
    '5HT2A': {
        'pdb_id': '9AS1',  # 5-HT2A receptor bound to DMT
        'description': 'Serotonin 5-HT2A receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
    'DOPAMINE_D2': {
        'pdb_id': '6A93',  # Dopamine D2 receptor
        'description': 'Dopamine D2 receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
    'GABA_A': {
        'pdb_id': '6D6T',  # GABA-A receptor alpha1-beta2-gamma2
        'description': 'GABA-A receptor (α1β2γ2)',
        'chain': 'A',
        'ligand_chain': None,
    },
    'GABA_B': {
        'pdb_id': '7EB2',  # GABA-B receptor
        'description': 'GABA-B receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
    'SEROTONIN_5HT1A': {
        'pdb_id': '7E2X',  # 5-HT1A receptor
        'description': 'Serotonin 5-HT1A receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
    'DOPAMINE_D3': {
        'pdb_id': '3PBL',  # Dopamine D3 receptor
        'description': 'Dopamine D3 receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
    'MUSCARINIC_M1': {
        'pdb_id': '5CXV',  # Muscarinic M1 receptor
        'description': 'Muscarinic M1 receptor',
        'chain': 'A',
        'ligand_chain': None,
    },
}

def create_receptor_dir():
    """Create receptors directory if it doesn't exist"""
    receptor_dir = Path(__file__).parent.parent / 'receptors'
    receptor_dir.mkdir(exist_ok=True)
    return receptor_dir

def download_pdb(pdb_id: str, output_path: Path) -> bool:
    """Download PDB file from RCSB PDB"""
    try:
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        print(f'  Downloading {pdb_id} from {url}...')
        
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(output_path, 'w') as f:
                f.write(response.text)
            print(f'  ✅ Downloaded: {output_path}')
            return True
        else:
            print(f'  ❌ Failed to download {pdb_id}: HTTP {response.status_code}')
            return False
    except Exception as e:
        print(f'  ❌ Error downloading {pdb_id}: {str(e)}')
        return False

def convert_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert PDB to PDBQT using Open Babel"""
    try:
        print(f'  Converting to PDBQT...')
        
        # Use obabel to convert PDB to PDBQT
        cmd = [
            'obabel',
            str(pdb_path),
            '-O', str(pdbqt_path),
            '-xr',  # Remove water
            '-p', '7.4',  # pH for protonation
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0 and pdbqt_path.exists():
            print(f'  ✅ Converted: {pdbqt_path}')
            return True
        else:
            print(f'  ❌ Conversion failed: {result.stderr}')
            return False
    except Exception as e:
        print(f'  ❌ Error converting to PDBQT: {str(e)}')
        return False

def create_metadata(target_name: str, target_info: dict, pdb_path: Path, pdbqt_path: Path) -> dict:
    """Create metadata for the receptor"""
    return {
        'name': target_name,
        'pdb_id': target_info['pdb_id'],
        'description': target_info['description'],
        'pdb_file': str(pdb_path.name),
        'pdbqt_file': str(pdbqt_path.name),
        'chain': target_info['chain'],
        'ligand_chain': target_info['ligand_chain'],
        'source': 'RCSB PDB',
        'downloaded_at': str(Path(pdb_path).stat().st_mtime),
    }

def main():
    print('🧬 Downloading Psychiatric Drug Target Receptors\n')
    print('=' * 60)
    
    receptor_dir = create_receptor_dir()
    print(f'Receptor directory: {receptor_dir}\n')
    
    metadata_list = []
    success_count = 0
    
    for target_name, target_info in PSYCHIATRIC_TARGETS.items():
        print(f'\n📥 Processing {target_name} (PDB: {target_info["pdb_id"]})')
        print(f'   Description: {target_info["description"]}')
        
        pdb_id = target_info['pdb_id']
        pdb_path = receptor_dir / f'{target_name}_{pdb_id}.pdb'
        pdbqt_path = receptor_dir / f'{target_name}_{pdb_id}.pdbqt'
        
        # Download PDB
        if not pdb_path.exists():
            if not download_pdb(pdb_id, pdb_path):
                continue
        else:
            print(f'  ℹ️  PDB file already exists: {pdb_path}')
        
        # Convert to PDBQT
        if not pdbqt_path.exists():
            if not convert_to_pdbqt(pdb_path, pdbqt_path):
                continue
        else:
            print(f'  ℹ️  PDBQT file already exists: {pdbqt_path}')
        
        # Create metadata
        metadata = create_metadata(target_name, target_info, pdb_path, pdbqt_path)
        metadata_list.append(metadata)
        success_count += 1
    
    # Save metadata
    metadata_file = receptor_dir / 'metadata.json'
    with open(metadata_file, 'w') as f:
        json.dump(metadata_list, f, indent=2)
    
    print(f'\n' + '=' * 60)
    print(f'✅ Downloaded {success_count}/{len(PSYCHIATRIC_TARGETS)} receptors')
    print(f'📄 Metadata saved to: {metadata_file}\n')
    
    # Print summary
    print('📋 Receptor Summary:')
    for metadata in metadata_list:
        print(f'  • {metadata["name"]}: {metadata["description"]}')
        print(f'    PDB: {metadata["pdb_file"]} → PDBQT: {metadata["pdbqt_file"]}')

if __name__ == '__main__':
    main()
