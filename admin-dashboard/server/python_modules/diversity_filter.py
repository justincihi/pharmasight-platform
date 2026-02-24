#!/usr/bin/env python3
"""
Diversity Filter for PharmaSight™
Filters discoveries based on structural diversity against existing database analogs
"""

import sys
import os

# Add current directory to path
src_dir = os.path.dirname(os.path.abspath(__file__))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from rdkit_analog_generator import RDKitAnalogGenerator
from typing import List, Dict

class DiversityFilter:
    """Filter discoveries based on structural diversity"""
    
    def __init__(self, min_diversity_threshold: float = 0.3):
        """
        Initialize diversity filter
        
        Args:
            min_diversity_threshold: Minimum Tanimoto distance required (default 0.3)
                                     0.3 means 70% similarity maximum
        """
        self.generator = RDKitAnalogGenerator()
        self.min_diversity_threshold = min_diversity_threshold
    
    def filter_discoveries(self, new_discoveries: List[Dict], 
                          existing_smiles: List[str]) -> List[Dict]:
        """
        Filter new discoveries to keep only structurally diverse ones
        
        Args:
            new_discoveries: List of discovery dicts with 'compound_smiles' field
            existing_smiles: List of SMILES strings already in database
        
        Returns:
            Filtered list of discoveries that meet diversity threshold
        """
        filtered = []
        
        for discovery in new_discoveries:
            smiles = discovery.get('compound_smiles', '')
            if not smiles:
                continue
            
            diversity_score, is_diverse = self.generator.calculate_diversity_score(
                smiles, 
                existing_smiles,
                self.min_diversity_threshold
            )
            
            # Add diversity metadata
            discovery['diversity_score'] = diversity_score
            discovery['is_structurally_diverse'] = is_diverse
            
            if is_diverse:
                filtered.append(discovery)
                # Add to existing list to check against future discoveries in this batch
                existing_smiles.append(smiles)
        
        return filtered


def filter_discoveries_from_json(discoveries_json: List[Dict], 
                                 existing_smiles_json: List[str],
                                 min_threshold: float = 0.3) -> List[Dict]:
    """
    Convenience function for filtering discoveries from JSON data
    
    Args:
        discoveries_json: List of discovery dictionaries
        existing_smiles_json: List of existing SMILES strings
        min_threshold: Minimum diversity threshold
    
    Returns:
        Filtered discoveries list
    """
    filter_engine = DiversityFilter(min_diversity_threshold=min_threshold)
    return filter_engine.filter_discoveries(discoveries_json, existing_smiles_json)


if __name__ == '__main__':
    # Test the diversity filter
    test_discoveries = [
        {'compound_name': 'Test-1', 'compound_smiles': 'CC(C)NCC(O)c1ccc(O)c(CO)c1'},
        {'compound_name': 'Test-2', 'compound_smiles': 'CN1C(=O)N(C)c2ncn(C)c2C1=O'},
    ]
    
    test_existing = ['CC(C)NCC(O)c1ccc(O)c(CO)c1']  # Same as Test-1
    
    filter_engine = DiversityFilter(min_diversity_threshold=0.3)
    filtered = filter_engine.filter_discoveries(test_discoveries, test_existing)
    
    print(f"Original discoveries: {len(test_discoveries)}")
    print(f"Filtered discoveries: {len(filtered)}")
    for disc in filtered:
        print(f"  {disc['compound_name']}: diversity={disc['diversity_score']}")
