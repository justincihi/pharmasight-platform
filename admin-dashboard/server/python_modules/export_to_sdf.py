#!/usr/bin/env python3
"""
Export analogs to SDF format with 3D coordinates
"""

import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem

def export_to_sdf(analogs_json: str, output_path: str):
    """Export analogs to SDF file"""
    try:
        analogs = json.loads(analogs_json)
        writer = Chem.SDWriter(output_path)
        
        for analog in analogs:
            smiles = analog.get('smiles', '')
            if not smiles:
                continue
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == 0:
                AllChem.UFFOptimizeMolecule(mol)
            
            # Set properties
            mol.SetProp('COMPOUND_ID', str(analog.get('compoundId', '')))
            mol.SetProp('COMPOUND_NAME', str(analog.get('compoundName', '')))
            mol.SetProp('PARENT_COMPOUND', str(analog.get('parentCompound', '')))
            mol.SetProp('SMILES', smiles)
            mol.SetProp('CONFIDENCE_SCORE', str(analog.get('confidenceScore', 0)))
            mol.SetProp('SAFETY_SCORE', str(analog.get('safetyScore', 0)))
            mol.SetProp('EFFICACY_SCORE', str(analog.get('efficacyScore', 0)))
            mol.SetProp('DRUG_LIKENESS_SCORE', str(analog.get('drugLikenessScore', 0)))
            mol.SetProp('PATENT_STATUS', str(analog.get('patentStatus', '')))
            mol.SetProp('THERAPEUTIC_POTENTIAL', str(analog.get('therapeuticPotential', '')))
            mol.SetProp('MECHANISM_OF_ACTION', str(analog.get('mechanismOfAction', '')))
            mol.SetProp('MOLECULAR_WEIGHT', str(analog.get('molecularWeight', '')))
            mol.SetProp('LOG_P', str(analog.get('logP', '')))
            mol.SetProp('H_BOND_DONORS', str(analog.get('hBondDonors', 0)))
            mol.SetProp('H_BOND_ACCEPTORS', str(analog.get('hBondAcceptors', 0)))
            mol.SetProp('DISCOVERED_AT', str(analog.get('discoveredAt', '')))
            mol.SetProp('APPROVAL_STATUS', str(analog.get('approvalStatus', '')))
            
            writer.write(mol)
        
        writer.close()
        print(json.dumps({'success': True, 'filePath': output_path}))
        
    except Exception as e:
        print(json.dumps({'success': False, 'error': str(e)}))
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(json.dumps({'success': False, 'error': 'Usage: export_to_sdf.py <analogs_json> <output_path>'}))
        sys.exit(1)
    
    export_to_sdf(sys.argv[1], sys.argv[2])
