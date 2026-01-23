#!/usr/bin/env python3
"""
AutoDock Vina Integration for Molecular Docking
Simulates drug-protein interactions
"""

import os
import json
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class AutoDockSimulator:
    """Integrate AutoDock Vina for molecular docking simulations"""
    
    def __init__(self):
        self.vina_executable = "vina"  # Assuming vina is installed
        self.protein_library = self.load_protein_targets()
        
    def load_protein_targets(self) -> Dict:
        """Load comprehensive drug target proteins from receptor database"""
        from src.receptor_database import RECEPTOR_DATABASE
        
        # Map receptors with PDB structures to docking targets
        targets = {}
        
        for receptor_name, receptor_data in RECEPTOR_DATABASE.items():
            if 'pdb_ids' in receptor_data and receptor_data['pdb_ids']:
                # Use first PDB ID as primary structure
                pdb_id = receptor_data['pdb_ids'][0] if receptor_data['pdb_ids'] else None
                
                if pdb_id:
                    # Estimate binding site (would be receptor-specific in production)
                    center = [30.0, 35.0, 40.0]  # Default center
                    size = [20, 20, 20]  # Default box size
                    
                    # Adjust for specific receptor types
                    if 'Ion channel' in receptor_data.get('type', ''):
                        size = [22, 22, 25]  # Larger for channels
                    elif 'Nuclear' in receptor_data.get('type', ''):
                        size = [18, 18, 18]  # Smaller for nuclear receptors
                    
                    targets[receptor_name] = {
                        "pdb_id": pdb_id,
                        "name": f"{receptor_data.get('family', '')} {receptor_name} receptor",
                        "description": receptor_data.get('therapeutic_relevance', 'Drug target'),
                        "binding_site": {"center": center, "size": size},
                        "receptor_type": receptor_data.get('type', 'Unknown'),
                        "family": receptor_data.get('family', 'Unknown'),
                        "known_ligands": {
                            "agonists": receptor_data.get('agonists', []),
                            "antagonists": receptor_data.get('antagonists', []),
                            "modulators": receptor_data.get('positive_modulators', [])
                        }
                    }
        
        # Keep some essential targets even without full receptor data
        if "5-HT2A" not in targets:
            targets["5-HT2A"] = {
                "pdb_id": "6A93",
                "name": "Serotonin 5-HT2A receptor",
                "description": "Target for psychedelics and antipsychotics",
                "binding_site": {"center": [25.0, 30.0, 40.0], "size": [20, 20, 20]}
            }
        
        return targets
    
    def prepare_ligand(self, smiles: str, output_path: str) -> bool:
        """Prepare ligand from SMILES for docking"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False
            
            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Write as PDBQT format (required by Vina)
            # For now, save as PDB and note that conversion to PDBQT is needed
            pdb_path = output_path.replace('.pdbqt', '.pdb')
            Chem.MolToPDBFile(mol, pdb_path)
            
            # Note: In production, use OpenBabel or ADFRsuite to convert PDB to PDBQT
            # For simulation, we'll create a mock PDBQT
            self._create_mock_pdbqt(pdb_path, output_path)
            
            return True
        except Exception as e:
            print(f"Error preparing ligand: {e}")
            return False
    
    def _create_mock_pdbqt(self, pdb_path: str, pdbqt_path: str):
        """Create a mock PDBQT file for demonstration"""
        # In production, use: obabel input.pdb -O output.pdbqt
        with open(pdbqt_path, 'w') as f:
            f.write("REMARK  Mock PDBQT file for demonstration\n")
            if os.path.exists(pdb_path):
                with open(pdb_path, 'r') as pdb:
                    for line in pdb:
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            f.write(line)
            f.write("TORSDOF 0\n")
    
    def simulate_docking(self, smiles: str, target_name: str) -> Dict:
        """Simulate molecular docking with a protein target"""
        if target_name not in self.protein_library:
            return {"error": f"Unknown target: {target_name}"}
        
        target = self.protein_library[target_name]
        
        # Since we don't have actual Vina installed, simulate results
        # In production, this would run actual Vina calculations
        
        # Simulate binding affinity based on molecular properties
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        # Simple scoring function (mock)
        mw = Chem.Descriptors.MolWt(mol)
        logp = Chem.Crippen.MolLogP(mol)
        hbd = Chem.Descriptors.NumHDonors(mol)
        hba = Chem.Descriptors.NumHAcceptors(mol)
        
        # Mock binding affinity calculation
        base_affinity = -6.0  # kcal/mol
        mw_penalty = max(0, (mw - 500) / 100) * 0.5
        logp_bonus = min(2.0, abs(logp - 2.5)) * 0.3
        hbond_bonus = min(hbd + hba, 10) * 0.1
        
        binding_affinity = base_affinity - mw_penalty - logp_bonus - hbond_bonus
        binding_affinity = round(binding_affinity + np.random.normal(0, 0.5), 2)
        
        # Classify binding strength
        if binding_affinity < -9.0:
            binding_class = "Excellent"
            confidence = 95
        elif binding_affinity < -7.0:
            binding_class = "Good"
            confidence = 85
        elif binding_affinity < -5.0:
            binding_class = "Moderate"
            confidence = 70
        else:
            binding_class = "Weak"
            confidence = 50
        
        # Generate mock binding poses
        poses = []
        for i in range(3):
            pose_affinity = binding_affinity + np.random.normal(0, 0.3)
            poses.append({
                "rank": i + 1,
                "affinity": round(pose_affinity, 2),
                "rmsd_lower": round(np.random.uniform(0, 2), 2),
                "rmsd_upper": round(np.random.uniform(2, 5), 2)
            })
        
        # Identify key interactions
        interactions = self._predict_interactions(mol, target_name)
        
        return {
            "success": True,
            "target": target,
            "ligand_smiles": smiles,
            "binding_affinity": binding_affinity,
            "binding_class": binding_class,
            "confidence": confidence,
            "poses": poses,
            "interactions": interactions,
            "recommendations": self._generate_recommendations(binding_affinity, interactions)
        }
    
    def _predict_interactions(self, mol, target_name: str) -> List[Dict]:
        """Predict molecular interactions with target"""
        interactions = []
        
        # Check for common pharmacophores
        smarts_patterns = {
            "H-bond donor": "[NH,NH2,NH3,OH]",
            "H-bond acceptor": "[O,N;!H]",
            "Aromatic": "a",
            "Positive charge": "[+,NH3+,NH2+,NH+]",
            "Negative charge": "[-,O-,COO-]",
            "Hydrophobic": "[CH3,CH2,CH,C]"
        }
        
        for interaction_type, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                interactions.append({
                    "type": interaction_type,
                    "count": len(matches),
                    "strength": "Strong" if len(matches) > 2 else "Moderate"
                })
        
        # Add target-specific interactions
        if "5-HT" in target_name:
            interactions.append({
                "type": "π-π stacking",
                "count": 1,
                "strength": "Strong",
                "residue": "Trp336"
            })
        
        return interactions
    
    def _generate_recommendations(self, affinity: float, interactions: List[Dict]) -> List[str]:
        """Generate optimization recommendations"""
        recommendations = []
        
        if affinity > -7.0:
            recommendations.append("Consider adding hydrogen bond donors/acceptors")
            recommendations.append("Optimize molecular weight (ideal: 300-500 Da)")
        
        if not any(i["type"] == "Aromatic" for i in interactions):
            recommendations.append("Add aromatic rings for π-π interactions")
        
        if affinity < -8.0:
            recommendations.append("Excellent binding predicted - proceed to synthesis")
            recommendations.append("Consider ADMET profiling for drug-likeness")
        
        return recommendations
    
    def batch_docking(self, smiles_list: List[str], target_name: str) -> List[Dict]:
        """Perform docking for multiple compounds"""
        results = []
        for smiles in smiles_list:
            result = self.simulate_docking(smiles, target_name)
            results.append(result)
        
        # Sort by binding affinity
        results.sort(key=lambda x: x.get('binding_affinity', 0))
        
        return results
    
    def virtual_screening(self, smiles_list: List[str], 
                         target_names: List[str] = None) -> Dict:
        """Screen compounds against multiple targets"""
        if not target_names:
            target_names = list(self.protein_library.keys())
        
        screening_results = {}
        
        for target in target_names:
            results = self.batch_docking(smiles_list, target)
            
            # Get top hits
            top_hits = [r for r in results if r.get('binding_affinity', 0) < -7.0]
            
            screening_results[target] = {
                "target_info": self.protein_library[target],
                "compounds_screened": len(smiles_list),
                "hits": len(top_hits),
                "top_compounds": top_hits[:5] if top_hits else []
            }
        
        return screening_results
    
    def generate_docking_report(self, docking_result: Dict) -> str:
        """Generate detailed docking report"""
        report = f"""
        ================== MOLECULAR DOCKING REPORT ==================
        
        Target Protein: {docking_result['target']['name']}
        PDB ID: {docking_result['target']['pdb_id']}
        Description: {docking_result['target']['description']}
        
        Ligand SMILES: {docking_result['ligand_smiles']}
        
        === BINDING RESULTS ===
        Binding Affinity: {docking_result['binding_affinity']} kcal/mol
        Binding Class: {docking_result['binding_class']}
        Confidence: {docking_result['confidence']}%
        
        === TOP BINDING POSES ===
        """
        
        for pose in docking_result.get('poses', []):
            report += f"""
        Pose {pose['rank']}:
          - Affinity: {pose['affinity']} kcal/mol
          - RMSD lower bound: {pose['rmsd_lower']} Å
          - RMSD upper bound: {pose['rmsd_upper']} Å
            """
        
        report += "\n=== MOLECULAR INTERACTIONS ===\n"
        for interaction in docking_result.get('interactions', []):
            report += f"- {interaction['type']}: {interaction['count']} sites ({interaction['strength']})\n"
        
        report += "\n=== RECOMMENDATIONS ===\n"
        for rec in docking_result.get('recommendations', []):
            report += f"• {rec}\n"
        
        report += "\n" + "="*60
        
        return report