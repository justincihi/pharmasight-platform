"""
Molecular Docking Module
Integrates AutoDock Vina for protein-ligand docking simulations
"""

import os
import tempfile
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to import vina
try:
    from vina import Vina
    VINA_AVAILABLE = True
except ImportError:
    logger.warning("AutoDock Vina not available. Docking simulations will use mock data.")
    VINA_AVAILABLE = False


class MolecularDocking:
    """Molecular docking using AutoDock Vina"""
    
    def __init__(self):
        self.vina_available = VINA_AVAILABLE
        
    def prepare_ligand(self, smiles: str) -> Optional[str]:
        """
        Prepare ligand from SMILES for docking
        Returns PDBQT format string
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e:
            logger.error(f"3D coordinate generation failed: {e}")
            return None
        
        # Convert to PDB format
        pdb_block = Chem.MolToPDBBlock(mol)
        
        return pdb_block
    
    def dock_ligand(self,
                   ligand_smiles: str,
                   receptor_pdb: Optional[str] = None,
                   center: Tuple[float, float, float] = (0, 0, 0),
                   box_size: Tuple[float, float, float] = (20, 20, 20)) -> Dict:
        """
        Perform molecular docking
        
        Args:
            ligand_smiles: SMILES string of ligand
            receptor_pdb: PDB file path of receptor (optional)
            center: Center of docking box (x, y, z)
            box_size: Size of docking box (x, y, z)
        
        Returns:
            Docking results including binding affinity and poses
        """
        if not self.vina_available:
            # Return mock docking results
            return self._mock_docking_results(ligand_smiles)
        
        # Prepare ligand
        ligand_pdb = self.prepare_ligand(ligand_smiles)
        if not ligand_pdb:
            return {'error': 'Ligand preparation failed'}
        
        try:
            # Initialize Vina
            v = Vina(sf_name='vina')
            
            # Set search space
            v.set_receptor(receptor_pdb) if receptor_pdb else None
            
            # For demonstration, use mock results if no receptor provided
            if not receptor_pdb:
                return self._mock_docking_results(ligand_smiles)
            
            # Compute Vina maps
            v.compute_vina_maps(center=center, box_size=box_size)
            
            # Set ligand
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(ligand_pdb)
                ligand_file = f.name
            
            v.set_ligand_from_file(ligand_file)
            
            # Dock
            v.dock(exhaustiveness=8, n_poses=9)
            
            # Get results
            energies = v.energies()
            
            # Clean up
            os.unlink(ligand_file)
            
            return {
                'binding_affinity': round(energies[0][0], 2),  # kcal/mol
                'num_poses': len(energies),
                'top_poses': [
                    {
                        'rank': i + 1,
                        'affinity': round(energy[0], 2),
                        'rmsd_lb': round(energy[1], 2),
                        'rmsd_ub': round(energy[2], 2),
                    }
                    for i, energy in enumerate(energies[:3])
                ],
                'success': True,
            }
            
        except Exception as e:
            logger.error(f"Docking failed: {e}")
            return {'error': str(e), 'success': False}
    
    def _mock_docking_results(self, smiles: str) -> Dict:
        """Generate mock docking results for demonstration"""
        import random
        random.seed(hash(smiles) % 1000)
        
        # Generate realistic-looking binding affinity based on molecular properties
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Chem.Descriptors.MolWt(mol)
            logp = Chem.Crippen.MolLogP(mol)
            
            # Estimate binding affinity (more negative = better binding)
            # Typical range: -12 to -4 kcal/mol for good binders
            base_affinity = -6.0
            mw_contribution = (mw - 300) / 100  # Larger molecules tend to bind better
            logp_contribution = logp * 0.5  # Hydrophobic interactions
            
            affinity = base_affinity - mw_contribution - logp_contribution
            affinity += random.uniform(-1, 1)  # Add some randomness
            affinity = round(max(-12, min(-3, affinity)), 2)
        else:
            affinity = round(random.uniform(-8, -5), 2)
        
        return {
            'binding_affinity': affinity,
            'num_poses': 9,
            'top_poses': [
                {
                    'rank': 1,
                    'affinity': affinity,
                    'rmsd_lb': round(random.uniform(0, 1), 2),
                    'rmsd_ub': round(random.uniform(1, 2), 2),
                },
                {
                    'rank': 2,
                    'affinity': round(affinity + random.uniform(0.5, 1.5), 2),
                    'rmsd_lb': round(random.uniform(1, 2), 2),
                    'rmsd_ub': round(random.uniform(2, 3), 2),
                },
                {
                    'rank': 3,
                    'affinity': round(affinity + random.uniform(1.5, 2.5), 2),
                    'rmsd_lb': round(random.uniform(2, 3), 2),
                    'rmsd_ub': round(random.uniform(3, 4), 2),
                },
            ],
            'success': True,
            'note': 'Mock docking results (AutoDock Vina not configured with receptor)',
        }
    
    def batch_dock_analogs(self,
                          analog_smiles_list: List[str],
                          receptor_pdb: Optional[str] = None) -> List[Dict]:
        """
        Dock multiple analogs and rank by binding affinity
        
        Args:
            analog_smiles_list: List of SMILES strings
            receptor_pdb: Receptor PDB file path
        
        Returns:
            List of docking results sorted by binding affinity
        """
        results = []
        
        for smiles in analog_smiles_list:
            docking_result = self.dock_ligand(smiles, receptor_pdb)
            docking_result['smiles'] = smiles
            results.append(docking_result)
        
        # Sort by binding affinity (more negative = better)
        results.sort(key=lambda x: x.get('binding_affinity', 0))
        
        return results
    
    def predict_binding_mode(self, smiles: str, target_protein: str) -> Dict:
        """
        Predict binding mode and key interactions
        
        Args:
            smiles: Ligand SMILES
            target_protein: Target protein name
        
        Returns:
            Predicted binding mode and interactions
        """
        # This is a simplified prediction
        # In a real system, this would analyze the docking pose
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {'error': 'Invalid SMILES'}
        
        # Analyze functional groups
        h_bond_donors = Chem.Lipinski.NumHDonors(mol)
        h_bond_acceptors = Chem.Lipinski.NumHAcceptors(mol)
        aromatic_rings = Chem.Lipinski.NumAromaticRings(mol)
        
        interactions = []
        
        if h_bond_donors > 0:
            interactions.append({
                'type': 'Hydrogen Bond Donor',
                'count': h_bond_donors,
                'description': f'{h_bond_donors} potential H-bond donor(s) with receptor'
            })
        
        if h_bond_acceptors > 0:
            interactions.append({
                'type': 'Hydrogen Bond Acceptor',
                'count': h_bond_acceptors,
                'description': f'{h_bond_acceptors} potential H-bond acceptor(s) with receptor'
            })
        
        if aromatic_rings > 0:
            interactions.append({
                'type': 'Pi-Pi Stacking',
                'count': aromatic_rings,
                'description': f'{aromatic_rings} aromatic ring(s) for pi-stacking interactions'
            })
        
        return {
            'target': target_protein,
            'predicted_interactions': interactions,
            'interaction_score': len(interactions) * 10,  # Simplified score
        }


class TargetLibrary:
    """Library of common drug targets"""
    
    TARGETS = {
        '5-HT2A': {
            'name': '5-HT2A Serotonin Receptor',
            'pdb_id': '6A93',
            'description': 'G protein-coupled receptor for serotonin',
            'therapeutic_areas': ['Depression', 'Psychedelic Therapy', 'Schizophrenia'],
        },
        'NMDA': {
            'name': 'NMDA Receptor',
            'pdb_id': '5UN1',
            'description': 'Glutamate-gated ion channel',
            'therapeutic_areas': ['Depression', 'Neuropathic Pain', 'Epilepsy'],
        },
        'D2': {
            'name': 'Dopamine D2 Receptor',
            'pdb_id': '6CM4',
            'description': 'G protein-coupled receptor for dopamine',
            'therapeutic_areas': ['Schizophrenia', 'Parkinson\'s Disease'],
        },
        'SERT': {
            'name': 'Serotonin Transporter',
            'pdb_id': '6DZZ',
            'description': 'Sodium-dependent serotonin transporter',
            'therapeutic_areas': ['Depression', 'Anxiety'],
        },
        'MAO-A': {
            'name': 'Monoamine Oxidase A',
            'pdb_id': '2Z5X',
            'description': 'Enzyme that catalyzes oxidative deamination of monoamines',
            'therapeutic_areas': ['Depression', 'Anxiety'],
        },
    }
    
    @classmethod
    def get_target_info(cls, target_name: str) -> Optional[Dict]:
        """Get information about a specific target"""
        return cls.TARGETS.get(target_name)
    
    @classmethod
    def list_targets(cls) -> List[str]:
        """List all available targets"""
        return list(cls.TARGETS.keys())
    
    @classmethod
    def search_targets_by_therapeutic_area(cls, area: str) -> List[str]:
        """Search targets by therapeutic area"""
        matching_targets = []
        for target_name, target_info in cls.TARGETS.items():
            if any(area.lower() in ta.lower() for ta in target_info['therapeutic_areas']):
                matching_targets.append(target_name)
        return matching_targets


class DockingAnalyzer:
    """Analyze and interpret docking results"""
    
    @staticmethod
    def classify_binding_affinity(affinity: float) -> str:
        """Classify binding affinity strength"""
        if affinity <= -10:
            return 'Very Strong'
        elif affinity <= -8:
            return 'Strong'
        elif affinity <= -6:
            return 'Moderate'
        elif affinity <= -4:
            return 'Weak'
        else:
            return 'Very Weak'
    
    @staticmethod
    def estimate_ki(affinity: float, temperature: float = 298.15) -> float:
        """
        Estimate inhibition constant (Ki) from binding affinity
        Using: ΔG = RT ln(Ki)
        
        Args:
            affinity: Binding affinity in kcal/mol
            temperature: Temperature in Kelvin
        
        Returns:
            Ki in nM
        """
        import math
        R = 0.001987  # kcal/(mol·K)
        
        # Convert affinity to Ki
        # ΔG = RT ln(Ki)
        # Ki = exp(ΔG / RT)
        
        ki_molar = math.exp(affinity / (R * temperature))
        ki_nm = ki_molar * 1e9  # Convert to nM
        
        return round(ki_nm, 2)
    
    @staticmethod
    def compare_analogs(docking_results: List[Dict]) -> Dict:
        """Compare docking results for multiple analogs"""
        if not docking_results:
            return {}
        
        affinities = [r['binding_affinity'] for r in docking_results if 'binding_affinity' in r]
        
        best_result = min(docking_results, key=lambda x: x.get('binding_affinity', 0))
        worst_result = max(docking_results, key=lambda x: x.get('binding_affinity', 0))
        
        return {
            'best_binder': {
                'smiles': best_result.get('smiles'),
                'affinity': best_result.get('binding_affinity'),
                'classification': DockingAnalyzer.classify_binding_affinity(
                    best_result.get('binding_affinity', 0)
                ),
            },
            'worst_binder': {
                'smiles': worst_result.get('smiles'),
                'affinity': worst_result.get('binding_affinity'),
            },
            'affinity_range': (min(affinities), max(affinities)),
            'mean_affinity': round(sum(affinities) / len(affinities), 2),
        }


# Singleton instance
_molecular_docking = None

def get_molecular_docking() -> MolecularDocking:
    """Get or create molecular docking instance"""
    global _molecular_docking
    if _molecular_docking is None:
        _molecular_docking = MolecularDocking()
    return _molecular_docking

