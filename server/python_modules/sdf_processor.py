#!/usr/bin/env python3
"""
SDF Processor for PharmaSight™
Parses SDF files, generates 3D coordinates, computes molecular descriptors
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem import rdPartialCharges
from typing import Dict, List, Optional, Tuple
import json
import os


class SDFProcessor:
    """Process SDF files and extract molecular properties"""
    
    def __init__(self):
        self.molecules = []
    
    def parse_sdf(self, sdf_path: str) -> List[Dict]:
        """
        Parse SDF file and extract all molecules with properties
        
        Args:
            sdf_path: Path to SDF file
        
        Returns:
            List of molecule dictionaries with properties
        """
        if not os.path.exists(sdf_path):
            raise FileNotFoundError(f"SDF file not found: {sdf_path}")
        
        supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
        molecules = []
        
        for idx, mol in enumerate(supplier):
            # Extract properties even if mol is None (empty structure)
            props = {}
            if mol is not None:
                props = mol.GetPropsAsDict()
            
            # Get SMILES from properties first (more reliable for empty structures)
            smiles = props.get('SMILES') or props.get('smiles')
            
            # If no SMILES in properties and mol exists, generate from structure
            if not smiles and mol is not None and mol.GetNumAtoms() > 0:
                smiles = Chem.MolToSmiles(mol)
            
            if not smiles:
                print(f"Warning: No SMILES found for molecule {idx}, skipping")
                continue
            
            # ONLY rebuild from SMILES if mol is None or invalid (preserves conformers!)
            if mol is None or mol.GetNumAtoms() == 0:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"Warning: Could not rebuild molecule {idx} from SMILES: {smiles}")
                    continue
            
            # Check if molecule has 3D coordinates
            has_3d_coords = self._has_3d_coordinates(mol)
            
            mol_data = {
                'index': idx,
                'smiles': smiles,
                'mol': mol,
                'has_3d_coords': has_3d_coords,
                'properties': props,
                'compound_name': props.get('COMPOUND_NAME') or props.get('compound_name') or props.get('name') or f'Compound_{idx}',
            }
            
            molecules.append(mol_data)
        
        self.molecules = molecules
        return molecules
    
    def _has_3d_coordinates(self, mol: Chem.Mol) -> bool:
        """Check if molecule has 3D coordinates"""
        try:
            # First check if molecule has any conformers
            if mol.GetNumConformers() == 0:
                return False
            
            conf = mol.GetConformer()
            # Check if Z coordinates have non-zero variance (indicates 3D)
            positions = conf.GetPositions()
            z_coords = positions[:, 2]
            return not all(abs(z) < 0.001 for z in z_coords)
        except:
            return False
    
    def generate_3d_coordinates(self, mol: Chem.Mol, optimize: bool = True) -> Chem.Mol:
        """
        Generate 3D coordinates for a molecule
        
        Args:
            mol: RDKit molecule object
            optimize: Whether to optimize geometry with MMFF
        
        Returns:
            Molecule with 3D coordinates
        """
        # Add hydrogens if not present
        mol_h = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
        
        if result == -1:
            # If embedding fails, try with random coordinates
            AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
        
        # Optimize geometry with MMFF force field
        if optimize:
            try:
                AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)
            except:
                # If MMFF fails, try UFF
                try:
                    AllChem.UFFOptimizeMolecule(mol_h, maxIters=200)
                except:
                    print("Warning: Could not optimize molecule geometry")
        
        return mol_h
    
    def compute_descriptors(self, mol: Chem.Mol) -> Dict:
        """
        Compute comprehensive molecular descriptors
        
        Args:
            mol: RDKit molecule object
        
        Returns:
            Dictionary of molecular descriptors
        """
        descriptors = {}
        
        try:
            # Basic properties
            descriptors['molecular_weight'] = round(Descriptors.MolWt(mol), 2)
            descriptors['exact_mass'] = round(Descriptors.ExactMolWt(mol), 2)
            
            # Lipophilicity
            descriptors['clogp'] = round(Descriptors.MolLogP(mol), 2)
            descriptors['molar_refractivity'] = round(Descriptors.MolMR(mol), 2)
            
            # Hydrogen bonding
            descriptors['hbd'] = rdMolDescriptors.CalcNumHBD(mol)
            descriptors['hba'] = rdMolDescriptors.CalcNumHBA(mol)
            
            # Polar surface area
            descriptors['tpsa'] = round(rdMolDescriptors.CalcTPSA(mol), 2)
            
            # Rotatable bonds (flexibility)
            descriptors['rotatable_bonds'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
            
            # Ring information
            descriptors['num_rings'] = rdMolDescriptors.CalcNumRings(mol)
            descriptors['num_aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
            descriptors['num_aliphatic_rings'] = rdMolDescriptors.CalcNumAliphaticRings(mol)
            
            # Atom counts
            descriptors['num_atoms'] = mol.GetNumAtoms()
            descriptors['num_heavy_atoms'] = mol.GetNumHeavyAtoms()
            descriptors['num_heteroatoms'] = rdMolDescriptors.CalcNumHeteroatoms(mol)
            
            # Charge-related (important for NMDA antagonists)
            descriptors['formal_charge'] = Chem.GetFormalCharge(mol)
            
            # Calculate partial charges for basicity estimation
            try:
                AllChem.ComputeGasteigerCharges(mol)
                charges = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') 
                          for i in range(mol.GetNumAtoms())]
                descriptors['max_positive_charge'] = round(max(charges), 3)
                descriptors['max_negative_charge'] = round(min(charges), 3)
            except:
                descriptors['max_positive_charge'] = None
                descriptors['max_negative_charge'] = None
            
            # Basicity indicator (count basic nitrogens)
            basic_n_count = 0
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'N':
                    # Check if nitrogen is likely basic (sp3, not aromatic, not amide)
                    if atom.GetHybridization() == Chem.HybridizationType.SP3:
                        basic_n_count += 1
            descriptors['basic_nitrogen_count'] = basic_n_count
            
            # Lipinski's Rule of Five violations
            descriptors['lipinski_violations'] = sum([
                descriptors['molecular_weight'] > 500,
                descriptors['clogp'] > 5,
                descriptors['hbd'] > 5,
                descriptors['hba'] > 10,
            ])
            
            # Drug-likeness score
            descriptors['drug_likeness_score'] = max(0, 100 - (descriptors['lipinski_violations'] * 20))
            
        except Exception as e:
            print(f"Error computing descriptors: {e}")
        
        return descriptors
    
    def write_sdf(self, molecules: List[Tuple[Chem.Mol, Dict]], output_path: str):
        """
        Write molecules to SDF file with properties
        
        Args:
            molecules: List of (mol, properties) tuples
            output_path: Output SDF file path
        """
        writer = Chem.SDWriter(output_path)
        
        for mol, props in molecules:
            # Set properties on molecule
            for key, value in props.items():
                if value is not None:
                    mol.SetProp(key, str(value))
            
            writer.write(mol)
        
        writer.close()
        print(f"Wrote {len(molecules)} molecules to {output_path}")
    
    def process_analog_sdf(self, sdf_path: str, output_dir: str = None) -> List[Dict]:
        """
        Complete pipeline: parse SDF, generate 3D coords, compute descriptors
        
        Args:
            sdf_path: Input SDF file path
            output_dir: Optional output directory for processed SDF
        
        Returns:
            List of processed analog dictionaries ready for database import
        """
        # Parse SDF
        molecules = self.parse_sdf(sdf_path)
        print(f"Parsed {len(molecules)} molecules from {sdf_path}")
        
        processed_analogs = []
        molecules_for_sdf = []
        
        for mol_data in molecules:
            mol = mol_data['mol']
            props = mol_data['properties']
            
            # Generate 3D coordinates if needed
            if not mol_data['has_3d_coords']:
                print(f"Generating 3D coordinates for {mol_data['compound_name']}")
                mol = self.generate_3d_coordinates(mol)
            
            # Compute descriptors
            descriptors = self.compute_descriptors(mol)
            
            # Extract scores from SDF properties
            efficacy_score = self._parse_numeric(props.get('efficacy_score') or props.get('efficacy'))
            safety_score = self._parse_numeric(props.get('safety_score') or props.get('safety'))
            confidence_score = self._parse_numeric(props.get('confidence_score') or props.get('confidence'))
            
            # Extract patent and value info
            patent_free = self._parse_boolean(props.get('patent-free') or props.get('patent_free'))
            ip_value = self._parse_value(props.get('IP_value_estimate') or props.get('market_value'))
            
            # Build analog dictionary
            analog = {
                'compound_name': mol_data['compound_name'],
                'smiles': mol_data['smiles'],
                'description': props.get('description') or '',
                'mechanism': props.get('mechanism') or 'NMDA receptor antagonist',
                'therapeutic_area': props.get('therapeutic_area') or 'CNS Disorders',
                
                # Scores
                'efficacy_score': efficacy_score or 85,
                'safety_score': safety_score or 85,
                'confidence_score': confidence_score or 85,
                
                # Patent info
                'patent_status': 'patent-free' if patent_free else 'patent-opportunity',
                'market_value': ip_value or 0,
                
                # Molecular descriptors
                **descriptors,
            }
            
            processed_analogs.append(analog)
            molecules_for_sdf.append((mol, {**props, **descriptors}))
        
        # Write processed SDF if output directory specified
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            base_name = os.path.splitext(os.path.basename(sdf_path))[0]
            output_path = os.path.join(output_dir, f"{base_name}_processed.sdf")
            self.write_sdf(molecules_for_sdf, output_path)
        
        return processed_analogs
    
    def _parse_numeric(self, value) -> Optional[float]:
        """Parse numeric value from string"""
        if value is None:
            return None
        try:
            return float(str(value).strip())
        except:
            return None
    
    def _parse_boolean(self, value) -> bool:
        """Parse boolean value from string"""
        if value is None:
            return False
        value_str = str(value).lower().strip()
        return value_str in ['true', '1', 'yes', 'patent-free']
    
    def _parse_value(self, value) -> float:
        """Parse monetary value from string (e.g., '$65M')"""
        if value is None:
            return 0.0
        try:
            value_str = str(value).strip().upper()
            # Remove currency symbols
            value_str = value_str.replace('$', '').replace('£', '').replace('€', '')
            
            # Handle M (millions), K (thousands), B (billions)
            multiplier = 1
            if 'M' in value_str:
                multiplier = 1_000_000
                value_str = value_str.replace('M', '')
            elif 'K' in value_str:
                multiplier = 1_000
                value_str = value_str.replace('K', '')
            elif 'B' in value_str:
                multiplier = 1_000_000_000
                value_str = value_str.replace('B', '')
            
            return float(value_str.strip()) * multiplier
        except:
            return 0.0


def process_sdf_file(sdf_path: str, output_dir: str = None) -> str:
    """
    Convenience function to process SDF file and return JSON
    
    Args:
        sdf_path: Input SDF file path
        output_dir: Optional output directory
    
    Returns:
        JSON string of processed analogs
    """
    processor = SDFProcessor()
    analogs = processor.process_analog_sdf(sdf_path, output_dir)
    return json.dumps(analogs, indent=2)


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python sdf_processor.py <sdf_file> [output_dir]")
        sys.exit(1)
    
    sdf_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    result = process_sdf_file(sdf_path, output_dir)
    print("\n=== Processed Analogs ===")
    print(result)
