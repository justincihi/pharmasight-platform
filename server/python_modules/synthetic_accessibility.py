"""
Synthetic Accessibility (SA) Score Calculator
Estimates how difficult it is to synthesize a molecule
Based on fragment contributions and complexity penalties
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski
import math
from typing import Dict

class SyntheticAccessibilityScorer:
    """
    Calculate synthetic accessibility scores for molecules
    Score range: 1 (easy to synthesize) to 10 (very difficult)
    """
    
    def __init__(self):
        # Complexity penalties
        self.size_penalty_threshold = 45  # Number of atoms
        self.stereo_penalty = 1.0
        self.spiro_penalty = 0.5
        self.bridge_penalty = 0.5
        self.macrocycle_penalty = 1.5
        
    def calculate_complexity_penalty(self, mol: Chem.Mol) -> float:
        """
        Calculate complexity penalty based on molecular features
        """
        penalty = 0.0
        
        # Size penalty
        num_atoms = mol.GetNumHeavyAtoms()
        if num_atoms > self.size_penalty_threshold:
            penalty += math.log10(num_atoms - self.size_penalty_threshold + 1)
        
        # Stereochemistry penalty
        num_stereo_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        num_unspecified_stereo = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        if num_stereo_centers > 0:
            penalty += self.stereo_penalty * math.log10(num_stereo_centers + 1)
        
        # Ring complexity penalties
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        
        # Spiro atoms
        num_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        if num_spiro > 0:
            penalty += self.spiro_penalty * num_spiro
        
        # Bridgehead atoms
        num_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        if num_bridgehead > 0:
            penalty += self.bridge_penalty * num_bridgehead
        
        # Macrocycles (rings with >8 atoms)
        for ring in ring_info.AtomRings():
            if len(ring) > 8:
                penalty += self.macrocycle_penalty
        
        # Fused ring systems
        if num_rings > 2:
            penalty += 0.3 * (num_rings - 2)
        
        return penalty
    
    def calculate_fragment_score(self, mol: Chem.Mol) -> float:
        """
        Calculate score based on common fragments
        (Simplified - in production, use pre-computed fragment library)
        """
        # Common functional groups (easier to synthesize)
        common_patterns = [
            'c1ccccc1',  # Benzene
            'C(=O)O',    # Carboxylic acid
            'C(=O)N',    # Amide
            'c1ccc(O)cc1',  # Phenol
            'C(=O)C',    # Ketone
            'COC',       # Ether
            'c1ccc(N)cc1',  # Aniline
        ]
        
        # Uncommon/reactive patterns (harder to synthesize)
        uncommon_patterns = [
            'C1OC1',     # Epoxide
            'C1NC1',     # Aziridine
            'N=[N+]=[N-]',  # Azide
            'C#N',       # Nitrile
            'C(F)(F)F',  # Trifluoromethyl
            '[N+](=O)[O-]',  # Nitro
        ]
        
        score = 0.0
        
        # Bonus for common fragments
        for pattern in common_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                score -= 0.2
        
        # Penalty for uncommon/reactive fragments
        for pattern in uncommon_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                score += 0.5
        
        return score
    
    def calculate_sa_score(self, smiles: str) -> Dict:
        """
        Calculate synthetic accessibility score
        Returns score from 1 (easy) to 10 (very difficult)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Base score
        base_score = 5.0
        
        # Complexity penalty
        complexity_penalty = self.calculate_complexity_penalty(mol)
        
        # Fragment score
        fragment_score = self.calculate_fragment_score(mol)
        
        # Additional factors
        num_atoms = mol.GetNumHeavyAtoms()
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        num_heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
        
        # Molecular complexity indicators
        if num_rotatable > 10:
            complexity_penalty += 0.5
        
        if num_heteroatoms / max(num_atoms, 1) > 0.5:
            complexity_penalty += 0.3
        
        # Calculate final score
        sa_score = base_score + complexity_penalty + fragment_score
        
        # Clamp to 1-10 range
        sa_score = max(1.0, min(10.0, sa_score))
        
        # Classify difficulty
        if sa_score <= 3:
            difficulty = 'Easy'
            synthesis_time = '1-2 steps'
            recommendation = 'Straightforward synthesis'
        elif sa_score <= 5:
            difficulty = 'Moderate'
            synthesis_time = '3-5 steps'
            recommendation = 'Achievable with standard methods'
        elif sa_score <= 7:
            difficulty = 'Challenging'
            synthesis_time = '6-10 steps'
            recommendation = 'Requires experienced synthetic chemist'
        else:
            difficulty = 'Very Difficult'
            synthesis_time = '>10 steps'
            recommendation = 'Consider structural simplification'
        
        return {
            'sa_score': round(sa_score, 2),
            'difficulty': difficulty,
            'estimated_steps': synthesis_time,
            'recommendation': recommendation,
            'num_atoms': num_atoms,
            'num_rings': num_rings,
            'num_stereocenters': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
            'num_rotatable_bonds': num_rotatable,
            'complexity_penalty': round(complexity_penalty, 2),
            'fragment_score': round(fragment_score, 2)
        }
    
    def compare_analogs(self, smiles_list: list) -> list:
        """
        Compare synthetic accessibility of multiple analogs
        """
        results = []
        for smiles in smiles_list:
            result = self.calculate_sa_score(smiles)
            result['smiles'] = smiles
            results.append(result)
        
        # Sort by SA score (easiest first)
        results.sort(key=lambda x: x.get('sa_score', 10))
        return results

def calculate_sa_score(smiles: str) -> Dict:
    """
    Main function to calculate SA score
    """
    scorer = SyntheticAccessibilityScorer()
    return scorer.calculate_sa_score(smiles)

def compare_analog_sa_scores(smiles_list: list) -> list:
    """
    Compare SA scores for multiple analogs
    """
    scorer = SyntheticAccessibilityScorer()
    return scorer.compare_analogs(smiles_list)

if __name__ == '__main__':
    # Test with different molecules
    test_molecules = [
        ('Aspirin', 'CC(=O)Oc1ccccc1C(=O)O'),
        ('Ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'),
        ('Taxol (complex)', 'CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C'),
    ]
    
    print('Synthetic Accessibility Scores:\n')
    for name, smiles in test_molecules:
        result = calculate_sa_score(smiles)
        if 'error' not in result:
            print(f"{name}:")
            print(f"  SA Score: {result['sa_score']} ({result['difficulty']})")
            print(f"  Estimated Steps: {result['estimated_steps']}")
            print(f"  Recommendation: {result['recommendation']}")
            print()
