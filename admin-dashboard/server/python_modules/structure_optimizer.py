"""
AI-Driven Structure Optimization Engine
Suggests molecular modifications to improve ADMET properties while maintaining activity
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, Lipinski
from typing import Dict, List, Tuple, Optional
import json

class StructureOptimizer:
    """
    Suggests structural modifications to optimize drug-like properties
    """
    
    def __init__(self):
        self.optimization_rules = self._load_optimization_rules()
        
    def _load_optimization_rules(self) -> Dict:
        """
        Load structure optimization rules based on medicinal chemistry principles
        """
        return {
            'reduce_lipophilicity': [
                {
                    'name': 'Add hydroxyl to aromatic ring',
                    'pattern': '[cH:1]',
                    'replacement': '[c:1]O',
                    'rationale': 'Increases polarity, reduces LogP',
                    'impact': {'logP': -0.5, 'solubility': +0.3}
                },
                {
                    'name': 'Replace phenyl with pyridine',
                    'pattern': 'c1ccccc1',
                    'replacement': 'c1ccncc1',
                    'rationale': 'Nitrogen reduces lipophilicity',
                    'impact': {'logP': -0.7, 'hERG': -0.2}
                },
                {
                    'name': 'Add carboxylic acid',
                    'pattern': '[CH3:1]',
                    'replacement': 'C(=O)O',
                    'rationale': 'Increases polarity and solubility',
                    'impact': {'logP': -1.5, 'solubility': +0.5}
                },
            ],
            'reduce_molecular_weight': [
                {
                    'name': 'Remove methyl group',
                    'pattern': '[C:1][CH3:2]',
                    'replacement': '[C:1]',
                    'rationale': 'Reduces MW by 14 Da',
                    'impact': {'mw': -14, 'size': -1}
                },
                {
                    'name': 'Replace ethyl with methyl',
                    'pattern': '[C:1][CH2:2][CH3:3]',
                    'replacement': '[C:1][CH3:2]',
                    'rationale': 'Reduces MW by 14 Da',
                    'impact': {'mw': -14, 'size': -1}
                },
            ],
            'improve_solubility': [
                {
                    'name': 'Add amino group',
                    'pattern': '[cH:1]',
                    'replacement': '[c:1]N',
                    'rationale': 'Increases HBA, improves solubility',
                    'impact': {'solubility': +0.4, 'hbd': +1}
                },
                {
                    'name': 'Add hydroxyl group',
                    'pattern': '[CH2:1]',
                    'replacement': '[CH:1]O',
                    'rationale': 'Increases polarity',
                    'impact': {'solubility': +0.3, 'hbd': +1}
                },
                {
                    'name': 'Replace alkyl with ether',
                    'pattern': '[CH2:1][CH2:2]',
                    'replacement': '[CH2:1]O[CH2:2]',
                    'rationale': 'Oxygen increases polarity',
                    'impact': {'solubility': +0.2, 'logP': -0.3}
                },
            ],
            'reduce_hERG_risk': [
                {
                    'name': 'Remove basic nitrogen',
                    'pattern': '[N:1]([CH3:2])([CH3:3])',
                    'replacement': '[N:1]([CH3:2])',
                    'rationale': 'Reduces basicity, lowers hERG risk',
                    'impact': {'hERG': -0.3}
                },
                {
                    'name': 'Add polar substituent',
                    'pattern': '[c:1][CH3:2]',
                    'replacement': '[c:1]C(=O)N',
                    'rationale': 'Reduces lipophilicity',
                    'impact': {'hERG': -0.2, 'logP': -0.5}
                },
            ],
            'improve_metabolic_stability': [
                {
                    'name': 'Block metabolic site with fluorine',
                    'pattern': '[CH3:1]',
                    'replacement': '[CH2:1]F',
                    'rationale': 'C-F bond resists metabolism',
                    'impact': {'stability': +0.3}
                },
                {
                    'name': 'Replace ester with amide',
                    'pattern': '[C:1](=O)[O:2][C:3]',
                    'replacement': '[C:1](=O)[N:2][C:3]',
                    'rationale': 'Amides more stable than esters',
                    'impact': {'stability': +0.4}
                },
            ],
            'reduce_toxicity': [
                {
                    'name': 'Remove aromatic amine',
                    'pattern': '[c:1][NH2:2]',
                    'replacement': '[c:1][NH:2]C(=O)C',
                    'rationale': 'Acetylation reduces mutagenic risk',
                    'impact': {'mutagenicity': -0.5}
                },
                {
                    'name': 'Remove nitro group',
                    'pattern': '[c:1][N+:2](=O)[O-]',
                    'replacement': '[c:1][NH2:2]',
                    'rationale': 'Nitro groups are mutagenic',
                    'impact': {'mutagenicity': -0.6}
                },
            ]
        }
    
    def analyze_properties(self, smiles: str) -> Dict:
        """
        Analyze current molecular properties
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        return {
            'mw': Descriptors.MolWt(mol),
            'logP': Crippen.MolLogP(mol),
            'hbd': rdMolDescriptors.CalcNumHBD(mol),
            'hba': rdMolDescriptors.CalcNumHBA(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
            'num_atoms': mol.GetNumHeavyAtoms()
        }
    
    def identify_issues(self, properties: Dict) -> List[str]:
        """
        Identify property issues that need optimization
        """
        issues = []
        
        if properties['logP'] > 5:
            issues.append('high_lipophilicity')
        if properties['mw'] > 500:
            issues.append('high_molecular_weight')
        if properties['logP'] < 0:
            issues.append('low_lipophilicity')
        if properties['tpsa'] > 140:
            issues.append('high_polarity')
        if properties['rotatable_bonds'] > 10:
            issues.append('high_flexibility')
        if properties['hbd'] > 5 or properties['hba'] > 10:
            issues.append('poor_solubility')
        
        return issues
    
    def generate_suggestions(self, smiles: str, target_property: Optional[str] = None) -> List[Dict]:
        """
        Generate optimization suggestions for a molecule
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [{'error': 'Invalid SMILES'}]
        
        properties = self.analyze_properties(smiles)
        issues = self.identify_issues(properties)
        
        suggestions = []
        
        # If target property specified, focus on that
        if target_property:
            rule_categories = [target_property] if target_property in self.optimization_rules else []
        else:
            # Auto-detect what needs optimization
            rule_categories = []
            if 'high_lipophilicity' in issues:
                rule_categories.append('reduce_lipophilicity')
            if 'high_molecular_weight' in issues:
                rule_categories.append('reduce_molecular_weight')
            if 'poor_solubility' in issues:
                rule_categories.append('improve_solubility')
            
            # Always check these
            rule_categories.extend(['reduce_hERG_risk', 'improve_metabolic_stability', 'reduce_toxicity'])
        
        # Apply rules
        for category in rule_categories:
            if category not in self.optimization_rules:
                continue
                
            for rule in self.optimization_rules[category]:
                pattern = Chem.MolFromSmarts(rule['pattern'])
                if pattern and mol.HasSubstructMatch(pattern):
                    # Try to apply the transformation
                    try:
                        rxn = AllChem.ReactionFromSmarts(f"{rule['pattern']}>>{rule['replacement']}")
                        products = rxn.RunReactants((mol,))
                        
                        if products and len(products) > 0:
                            product_mol = products[0][0]
                            Chem.SanitizeMol(product_mol)
                            product_smiles = Chem.MolToSmiles(product_mol)
                            
                            # Calculate new properties
                            new_properties = self.analyze_properties(product_smiles)
                            
                            suggestions.append({
                                'modification': rule['name'],
                                'category': category,
                                'rationale': rule['rationale'],
                                'original_smiles': smiles,
                                'optimized_smiles': product_smiles,
                                'property_changes': {
                                    'mw_change': round(new_properties['mw'] - properties['mw'], 1),
                                    'logP_change': round(new_properties['logP'] - properties['logP'], 2),
                                    'tpsa_change': round(new_properties['tpsa'] - properties['tpsa'], 1),
                                },
                                'original_properties': properties,
                                'new_properties': new_properties,
                                'priority': self._calculate_priority(properties, new_properties, issues)
                            })
                    except Exception as e:
                        # Skip invalid transformations
                        continue
        
        # Sort by priority
        suggestions.sort(key=lambda x: x['priority'], reverse=True)
        
        # Limit to top 10 suggestions
        return suggestions[:10]
    
    def _calculate_priority(self, old_props: Dict, new_props: Dict, issues: List[str]) -> float:
        """
        Calculate priority score for a suggestion
        Higher score = more beneficial modification
        """
        score = 0.0
        
        # Reward improvements in problematic areas
        if 'high_lipophilicity' in issues and new_props['logP'] < old_props['logP']:
            score += 2.0
        if 'high_molecular_weight' in issues and new_props['mw'] < old_props['mw']:
            score += 1.5
        if 'poor_solubility' in issues and new_props['tpsa'] > old_props['tpsa']:
            score += 1.5
        
        # General improvements
        if new_props['logP'] < old_props['logP'] and old_props['logP'] > 3:
            score += 1.0
        if new_props['mw'] < old_props['mw'] and old_props['mw'] > 400:
            score += 0.5
        
        # Penalize if making things worse
        if new_props['mw'] > 600:
            score -= 2.0
        if new_props['logP'] > 6:
            score -= 2.0
        if new_props['logP'] < -2:
            score -= 1.5
        
        return score
    
    def optimize_for_target(self, smiles: str, target_property: str, iterations: int = 3) -> List[Dict]:
        """
        Iteratively optimize a molecule for a specific property
        """
        current_smiles = smiles
        optimization_history = []
        
        for i in range(iterations):
            suggestions = self.generate_suggestions(current_smiles, target_property)
            
            if not suggestions or len(suggestions) == 0:
                break
            
            # Take the best suggestion
            best = suggestions[0]
            optimization_history.append({
                'iteration': i + 1,
                'modification': best['modification'],
                'smiles': best['optimized_smiles'],
                'properties': best['new_properties'],
                'rationale': best['rationale']
            })
            
            current_smiles = best['optimized_smiles']
        
        return optimization_history

def optimize_structure(smiles: str, target_property: Optional[str] = None) -> List[Dict]:
    """
    Main function to get optimization suggestions
    """
    optimizer = StructureOptimizer()
    return optimizer.generate_suggestions(smiles, target_property)

def iterative_optimization(smiles: str, target_property: str, iterations: int = 3) -> List[Dict]:
    """
    Perform iterative optimization
    """
    optimizer = StructureOptimizer()
    return optimizer.optimize_for_target(smiles, target_property, iterations)

if __name__ == '__main__':
    # Test with a sample molecule
    test_smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'  # Ibuprofen
    print('Optimization Suggestions for Ibuprofen:\n')
    
    optimizer = StructureOptimizer()
    suggestions = optimizer.generate_suggestions(test_smiles)
    
    for i, sug in enumerate(suggestions[:5], 1):
        print(f"{i}. {sug['modification']}")
        print(f"   Category: {sug['category']}")
        print(f"   Rationale: {sug['rationale']}")
        print(f"   LogP change: {sug['property_changes']['logP_change']}")
        print(f"   MW change: {sug['property_changes']['mw_change']}")
        print(f"   Priority: {sug['priority']:.2f}")
        print()
