#!/usr/bin/env python3
"""
Compound Name Resolution System
Converts compound names (with or without salts) to SMILES strings
"""

import re
import json

class CompoundNameResolver:
    def __init__(self):
        # Comprehensive compound database with names, synonyms, and SMILES
        self.compound_database = {
            # Psychedelics
            'psilocybin': {
                'smiles': 'CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12',
                'molecular_weight': 284.25,
                'synonyms': ['4-phosphoryloxy-N,N-dimethyltryptamine', 'psilocybin phosphate']
            },
            'psilocin': {
                'smiles': 'CN(C)CCc1c[nH]c2ccc(O)cc12',
                'molecular_weight': 204.27,
                'synonyms': ['4-hydroxy-N,N-dimethyltryptamine', '4-HO-DMT']
            },
            'lsd': {
                'smiles': 'CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C',
                'molecular_weight': 323.43,
                'synonyms': ['lysergic acid diethylamide', 'acid', 'lucy']
            },
            'dmt': {
                'smiles': 'CN(C)CCc1c[nH]c2ccccc12',
                'molecular_weight': 188.27,
                'synonyms': ['N,N-dimethyltryptamine', 'dimethyltryptamine']
            },
            'mdma': {
                'smiles': 'CC(Cc1ccc2c(c1)OCO2)NC',
                'molecular_weight': 193.25,
                'synonyms': ['3,4-methylenedioxymethamphetamine', 'ecstasy', 'molly']
            },
            '2c-b': {
                'smiles': 'COc1cc(CCN)cc(Br)c1OC',
                'molecular_weight': 260.13,
                'synonyms': ['4-bromo-2,5-dimethoxyphenethylamine', '2C-B']
            },
            'mescaline': {
                'smiles': 'COc1cc(CCN)cc(OC)c1OC',
                'molecular_weight': 211.26,
                'synonyms': ['3,4,5-trimethoxyphenethylamine']
            },
            
            # Ketamine and analogs
            'ketamine': {
                'smiles': 'CNC1(c2ccccc2Cl)CCCCC1=O',
                'molecular_weight': 237.73,
                'synonyms': ['2-(2-chlorophenyl)-2-(methylamino)cyclohexanone']
            },
            'esketamine': {
                'smiles': 'CN[C@]1(CCCCC1=O)c2ccccc2Cl',
                'molecular_weight': 237.73,
                'synonyms': ['S-ketamine', '(S)-ketamine']
            },
            'arketamine': {
                'smiles': 'CN[C@@]1(CCCCC1=O)c2ccccc2Cl',
                'molecular_weight': 237.73,
                'synonyms': ['R-ketamine', '(R)-ketamine']
            },
            'norketamine': {
                'smiles': 'N[C@]1(CCCCC1=O)c2ccccc2Cl',
                'molecular_weight': 223.70,
                'synonyms': ['desmethylketamine']
            },
            
            # Opioids
            'morphine': {
                'smiles': 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5',
                'molecular_weight': 285.34,
                'synonyms': ['morphine sulfate', 'morphine hydrochloride']
            },
            'oxycodone': {
                'smiles': 'COc1ccc2c3c1O[C@H]1[C@@H](O)C=C[C@H]4[C@@H](C2)N(C)CC[C@@]341',
                'molecular_weight': 315.36,
                'synonyms': ['oxycodone hydrochloride', 'oxycontin']
            },
            'fentanyl': {
                'smiles': 'CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1',
                'molecular_weight': 336.47,
                'synonyms': ['fentanyl citrate', 'fentanyl hydrochloride']
            },
            'buprenorphine': {
                'smiles': 'COc1ccc2c3c1O[C@H]1[C@@H](O)C=C[C@H]4[C@@H](C2)N(CC1)CC[C@@]43C(C)(C)C',
                'molecular_weight': 467.64,
                'synonyms': ['buprenorphine hydrochloride', 'suboxone']
            },
            'tramadol': {
                'smiles': 'COc1cccc(C2(O)CCCCC2CN(C)C)c1',
                'molecular_weight': 263.38,
                'synonyms': ['tramadol hydrochloride', 'ultram']
            },
            
            # Benzodiazepines and analogs
            'alprazolam': {
                'smiles': 'Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2',
                'molecular_weight': 308.76,
                'synonyms': ['xanax', 'alprazolam tablets']
            },
            'diazepam': {
                'smiles': 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21',
                'molecular_weight': 284.74,
                'synonyms': ['valium', 'diazepam tablets']
            },
            'lorazepam': {
                'smiles': 'O=C1CN=C(c2ccc(Cl)cc2Cl)c2cc(Cl)ccc2N1O',
                'molecular_weight': 321.16,
                'synonyms': ['ativan', 'lorazepam tablets']
            },
            'clonazepam': {
                'smiles': 'O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1',
                'molecular_weight': 315.71,
                'synonyms': ['klonopin', 'clonazepam tablets']
            },
            'etizolam': {
                'smiles': 'CCc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2',
                'molecular_weight': 342.83,
                'synonyms': ['etizolam tablets']
            },
            'flualprazolam': {
                'smiles': 'Cc1nnc2n1-c1ccc(F)cc1C(c1ccccc1)=NC2',
                'molecular_weight': 326.75,
                'synonyms': ['flualprazolam powder']
            },
            'bromazolam': {
                'smiles': 'Cc1nnc2n1-c1ccc(Br)cc1C(c1ccccc1)=NC2',
                'molecular_weight': 353.21,
                'synonyms': ['bromazolam powder']
            },
            
            # Antidepressants
            'sertraline': {
                'smiles': 'CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21',
                'molecular_weight': 306.23,
                'synonyms': ['sertraline hydrochloride', 'zoloft']
            },
            'fluoxetine': {
                'smiles': 'CNCCC(c1ccc(C(F)(F)F)cc1)c1ccccc1',
                'molecular_weight': 309.33,
                'synonyms': ['fluoxetine hydrochloride', 'prozac']
            },
            'paroxetine': {
                'smiles': 'Fc1ccc([C@@H]2CCNCC2COc2ccc3c(c2)OCO3)cc1',
                'molecular_weight': 329.37,
                'synonyms': ['paroxetine hydrochloride', 'paxil']
            },
            'venlafaxine': {
                'smiles': 'COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1',
                'molecular_weight': 277.40,
                'synonyms': ['venlafaxine hydrochloride', 'effexor']
            },
            
            # Antipsychotics
            'aripiprazole': {
                'smiles': 'O=C1CCc2c(cccc2Cl)N1CCCCN1CCN(c2cccc(Cl)c2Cl)CC1',
                'molecular_weight': 448.39,
                'synonyms': ['abilify', 'aripiprazole tablets']
            },
            'risperidone': {
                'smiles': 'CC1=C(CCN2CCC(c3noc4cc(F)ccc34)CC2)C(=O)N2CCCCC2=N1',
                'molecular_weight': 410.49,
                'synonyms': ['risperdal', 'risperidone tablets']
            },
            'quetiapine': {
                'smiles': 'OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1',
                'molecular_weight': 383.51,
                'synonyms': ['quetiapine fumarate', 'seroquel']
            },
            'olanzapine': {
                'smiles': 'CN1CCN(C2=Nc3cc(C)ccc3Nc3ccccc32)CC1',
                'molecular_weight': 312.44,
                'synonyms': ['zyprexa', 'olanzapine tablets']
            },
            
            # Stimulants
            'amphetamine': {
                'smiles': 'CC(N)Cc1ccccc1',
                'molecular_weight': 135.21,
                'synonyms': ['amphetamine sulfate', 'adderall']
            },
            'methylphenidate': {
                'smiles': 'COC(=O)[C@H](c1ccccc1)[C@H]1CCCCN1',
                'molecular_weight': 233.31,
                'synonyms': ['methylphenidate hydrochloride', 'ritalin', 'concerta']
            },
            'modafinil': {
                'smiles': 'NC(=O)C[S@](=O)c1ccc(C(F)(F)F)cc1',
                'molecular_weight': 273.35,
                'synonyms': ['provigil', 'modafinil tablets']
            },
            
            # Novel compounds
            'dmnpc': {
                'smiles': 'CN(C)CCc1c[nH]c2ccc(OC(=O)N3CCCC3)cc12',
                'molecular_weight': 315.37,
                'synonyms': ['4-(pyrrolidin-1-ylcarbonyloxy)-N,N-dimethyltryptamine']
            },
            'ethylbromazelam': {
                'smiles': 'CCc1nnc2n1-c1ccc(Br)cc1C(c1ccccc1)=NC2',
                'molecular_weight': 367.24,
                'synonyms': ['ethyl-bromazolam']
            }
        }
        
        # Salt and form suffixes to remove
        self.salt_suffixes = [
            'hcl', 'hydrochloride', 'hydrochloric acid',
            'sulfate', 'sulphate', 'so4',
            'phosphate', 'po4',
            'citrate', 'tartrate', 'fumarate', 'maleate',
            'acetate', 'succinate', 'lactate',
            'free base', 'freebase', 'base',
            'salt', 'monohydrate', 'dihydrate',
            'tablets', 'capsules', 'powder', 'solution'
        ]
    
    def normalize_compound_name(self, compound_name):
        """Normalize compound name by removing salts and common suffixes"""
        if not compound_name:
            return ""
        
        # Convert to lowercase and strip whitespace
        normalized = compound_name.lower().strip()
        
        # Remove common punctuation and extra spaces
        normalized = re.sub(r'[^\w\s-]', ' ', normalized)
        normalized = re.sub(r'\s+', ' ', normalized).strip()
        
        # Remove salt suffixes
        for suffix in self.salt_suffixes:
            # Remove suffix at the end
            if normalized.endswith(' ' + suffix):
                normalized = normalized[:-len(' ' + suffix)].strip()
            # Remove suffix in the middle (e.g., "compound hcl tablets")
            normalized = re.sub(r'\s+' + re.escape(suffix) + r'\s+', ' ', normalized)
            normalized = re.sub(r'\s+' + re.escape(suffix) + r'$', '', normalized)
        
        # Handle common abbreviations
        abbreviation_map = {
            'lsd-25': 'lsd',
            '25i-nbome': '25i',
            '25b-nbome': '25b',
            '25c-nbome': '25c',
            '2cb': '2c-b',
            '2ci': '2c-i',
            '4-aco-dmt': '4-acetoxy-dmt',
            '4-ho-dmt': 'psilocin',
            '5-meo-dmt': '5-methoxy-dmt'
        }
        
        if normalized in abbreviation_map:
            normalized = abbreviation_map[normalized]
        
        return normalized
    
    def find_compound_by_name(self, compound_name):
        """Find compound data by name or synonym"""
        normalized_name = self.normalize_compound_name(compound_name)
        
        # Direct match
        if normalized_name in self.compound_database:
            return self.compound_database[normalized_name]
        
        # Search synonyms
        for compound_key, compound_data in self.compound_database.items():
            # Check if normalized name matches compound key
            if normalized_name == compound_key:
                return compound_data
            
            # Check synonyms
            for synonym in compound_data.get('synonyms', []):
                if normalized_name == self.normalize_compound_name(synonym):
                    return compound_data
        
        return None
    
    def resolve_compound_name(self, compound_name):
        """Resolve compound name to SMILES string and properties"""
        compound_data = self.find_compound_by_name(compound_name)
        
        if compound_data:
            return {
                'found': True,
                'smiles': compound_data['smiles'],
                'molecular_weight': compound_data['molecular_weight'],
                'synonyms': compound_data.get('synonyms', []),
                'normalized_name': self.normalize_compound_name(compound_name)
            }
        else:
            return {
                'found': False,
                'error': f'Compound "{compound_name}" not found in database',
                'suggestions': self.get_suggestions(compound_name),
                'normalized_name': self.normalize_compound_name(compound_name)
            }
    
    def get_suggestions(self, compound_name, max_suggestions=5):
        """Get suggestions for similar compound names"""
        normalized_input = self.normalize_compound_name(compound_name)
        suggestions = []
        
        for compound_key, compound_data in self.compound_database.items():
            # Check similarity with compound key
            if self.calculate_similarity(normalized_input, compound_key) > 0.6:
                suggestions.append(compound_key.title())
            
            # Check similarity with synonyms
            for synonym in compound_data.get('synonyms', []):
                normalized_synonym = self.normalize_compound_name(synonym)
                if self.calculate_similarity(normalized_input, normalized_synonym) > 0.6:
                    suggestions.append(synonym.title())
        
        # Remove duplicates and limit results
        suggestions = list(set(suggestions))[:max_suggestions]
        return suggestions
    
    def calculate_similarity(self, str1, str2):
        """Calculate simple string similarity"""
        if not str1 or not str2:
            return 0.0
        
        # Simple Jaccard similarity based on character n-grams
        def get_ngrams(s, n=2):
            return set(s[i:i+n] for i in range(len(s)-n+1))
        
        ngrams1 = get_ngrams(str1)
        ngrams2 = get_ngrams(str2)
        
        if not ngrams1 and not ngrams2:
            return 1.0
        if not ngrams1 or not ngrams2:
            return 0.0
        
        intersection = len(ngrams1.intersection(ngrams2))
        union = len(ngrams1.union(ngrams2))
        
        return intersection / union if union > 0 else 0.0
    
    def add_compound(self, name, smiles, molecular_weight, synonyms=None):
        """Add a new compound to the database"""
        normalized_name = self.normalize_compound_name(name)
        self.compound_database[normalized_name] = {
            'smiles': smiles,
            'molecular_weight': molecular_weight,
            'synonyms': synonyms or []
        }
    
    def get_all_compounds(self):
        """Get list of all available compounds"""
        compounds = []
        for compound_key, compound_data in self.compound_database.items():
            compounds.append({
                'name': compound_key.title(),
                'smiles': compound_data['smiles'],
                'molecular_weight': compound_data['molecular_weight'],
                'synonyms': compound_data.get('synonyms', [])
            })
        return compounds

# Global instance
compound_resolver = CompoundNameResolver()

def resolve_compound_name_to_smiles(compound_name):
    """Convenience function to resolve compound name to SMILES"""
    result = compound_resolver.resolve_compound_name(compound_name)
    if result['found']:
        return result['smiles']
    else:
        return None

def get_compound_suggestions(compound_name):
    """Convenience function to get compound suggestions"""
    result = compound_resolver.resolve_compound_name(compound_name)
    if not result['found']:
        return result.get('suggestions', [])
    return []

