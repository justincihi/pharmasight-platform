#!/usr/bin/env python3
"""
PharmaSight Cheminformatics Workflow Integration
Executes similarity-based screening, patent checks, and analog generation
"""

import sys
import json
import time
import requests
from typing import List, Dict, Optional, Set, Tuple

try:
    import pubchempy as pcp
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs, BRICS, Descriptors
except ImportError as e:
    print(json.dumps({"error": f"Missing dependency: {e}"}))
    sys.exit(1)


# ============================================================================
# WORKFLOW 1: Similarity-Based PubChem Screening
# ============================================================================

def confirm_and_fetch_similars(
    name_or_smiles: str,
    threshold: float = 0.70,
    max_hits: int = 25
) -> Dict:
    """
    Step 1: Confirm/canonicalize SMILES via PubChem
    Step 2: Run PubChem similarity search (Tanimoto >= threshold)
    Step 3: Return up to max_hits screened, confirmed compounds
    """
    try:
        # --- Confirm parent compound ---
        results = pcp.get_compounds(name_or_smiles, 'name')
        if not results:
            results = pcp.get_compounds(name_or_smiles, 'smiles')
        
        if not results:
            return {
                "success": False,
                "error": f"Compound '{name_or_smiles}' not found in PubChem."
            }
        
        parent = results[0]
        parent_smiles = parent.canonical_smiles
        parent_mol = Chem.MolFromSmiles(parent_smiles)
        
        if parent_mol is None:
            return {
                "success": False,
                "error": f"Invalid SMILES for parent compound: {parent_smiles}"
            }
        
        parent_fp = AllChem.GetMorganFingerprintAsBitVect(parent_mol, radius=2, nBits=2048)
        
        print(f"Confirmed: {parent.iupac_name or name_or_smiles}")
        print(f"Canonical SMILES: {parent_smiles}")
        print(f"PubChem CID: {parent.cid}")
        
        # --- PubChem similarity search (server-side Tanimoto) ---
        pct = int(threshold * 100)
        similar_hits = pcp.get_compounds(
            parent_smiles,
            namespace='smiles',
            searchtype='similarity',
            Threshold=pct,
            MaxRecords=max_hits
        )
        
        # --- Local RDKit re-scoring + filtering ---
        screened = []
        for hit in similar_hits:
            hit_mol = Chem.MolFromSmiles(hit.canonical_smiles)
            if hit_mol is None:
                continue
            
            hit_fp = AllChem.GetMorganFingerprintAsBitVect(hit_mol, radius=2, nBits=2048)
            tanimoto = DataStructs.TanimotoSimilarity(parent_fp, hit_fp)
            
            if tanimoto >= threshold:
                screened.append({
                    'cid': hit.cid,
                    'name': hit.iupac_name or 'Unknown',
                    'smiles': hit.canonical_smiles,
                    'mw': hit.molecular_weight,
                    'tanimoto': round(tanimoto, 4)
                })
            
            time.sleep(0.3)  # Rate limit
        
        screened = sorted(screened, key=lambda x: x['tanimoto'], reverse=True)
        
        return {
            "success": True,
            "parent_smiles": parent_smiles,
            "parent_cid": parent.cid,
            "parent_name": parent.iupac_name or name_or_smiles,
            "hits": screened[:max_hits],
            "total_hits": len(screened)
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


# ============================================================================
# WORKFLOW 2: Patent Status Check
# ============================================================================

def check_patent_status(cid: int) -> Dict:
    """
    Returns list of patent IDs linked to a PubChem CID.
    Empty list = no patent associations found (patent-free indicator).
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PatentID/JSON"
        resp = requests.get(url, timeout=10)
        
        if resp.status_code == 404:
            return {
                "success": True,
                "cid": cid,
                "patents": [],
                "patent_free": True
            }
        
        data = resp.json()
        patents = data['InformationList']['Information'][0].get('PatentID', [])
        
        return {
            "success": True,
            "cid": cid,
            "patents": patents,
            "patent_free": len(patents) == 0
        }
    
    except Exception as e:
        return {
            "success": False,
            "cid": cid,
            "error": str(e)
        }


def screen_and_flag_for_masterlist(
    hits: List[Dict],
    max_hits: int = 25
) -> Dict:
    """
    Cross-references each hit against PubChem patent data.
    Returns patent-free compounds ready for master list entry.
    """
    try:
        master_candidates = []
        
        for hit in hits:
            time.sleep(0.3)  # Respect PubChem rate limits
            
            patents_result = check_patent_status(hit['cid'])
            
            if patents_result['success']:
                hit['patents'] = patents_result['patents']
                hit['patent_free'] = patents_result['patent_free']
                hit['flag'] = '✅ CLEAR' if hit['patent_free'] else f"⚠️ PATENTED ({len(patents_result['patents'])} patents)"
                
                if hit['patent_free']:
                    master_candidates.append(hit)
                
                print(f"CID {hit['cid']} | {hit['tanimoto']:.2%} similar | {hit['flag']}")
        
        print(f"\n{len(master_candidates)} patent-free compounds flagged for master list.")
        
        return {
            "success": True,
            "master_candidates": master_candidates[:max_hits],
            "total_candidates": len(master_candidates)
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


# ============================================================================
# WORKFLOW 3: RDKit-Driven Analog Generation (BRICS Fragmentation)
# ============================================================================

def generate_brics_analogs(
    smiles: str,
    n: int = 25
) -> Dict:
    """
    Fragments the input mol via BRICS rules, then reassembles fragment
    combinations to produce novel analogs.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "success": False,
                "error": "Invalid SMILES."
            }
        
        # Break into BRICS fragments
        fragments = list(BRICS.BRICSDecompose(mol))
        print(f"Fragments: {fragments}")
        
        # Build new molecules from fragment pool
        fragment_mols = [Chem.MolFromSmiles(f) for f in fragments if Chem.MolFromSmiles(f)]
        analogs_raw = list(BRICS.BRICSBuild(fragment_mols))
        
        # Canonicalize, deduplicate, and score similarity to parent
        parent_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        analog_list = []
        seen: Set[str] = set()
        
        for analog in analogs_raw:
            if analog is None:
                continue
            
            smi = Chem.MolToSmiles(analog)
            if smi in seen or smi == Chem.MolToSmiles(mol):
                continue
            
            seen.add(smi)
            a_fp = AllChem.GetMorganFingerprintAsBitVect(analog, radius=2, nBits=2048)
            tanimoto = DataStructs.TanimotoSimilarity(parent_fp, a_fp)
            
            analog_list.append({
                'smiles': smi,
                'tanimoto': round(tanimoto, 4),
                'mw': Descriptors.MolWt(analog)
            })
        
        analog_list = sorted(analog_list, key=lambda x: x['tanimoto'], reverse=True)
        
        return {
            "success": True,
            "parent_smiles": smiles,
            "fragments": fragments,
            "analogs": analog_list[:n],
            "total_generated": len(analog_list)
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


# ============================================================================
# WORKFLOW 4: RDKit-Driven Analog Generation (R-Group Substitution)
# ============================================================================

SUBSTITUENTS = {
    'methoxy': 'OC',
    'fluoro': 'F',
    'chloro': 'Cl',
    'methyl': 'C',
    'trifluoromethyl': 'C(F)(F)F',
    'hydroxy': 'O',
    'amino': 'N',
    'ethyl': 'CC',
    'acetyl': 'C(=O)C',
    'cyano': 'C#N',
}


def enumerate_substituent_analogs(
    base_smiles: str,
    attachment_idx: int,
    n: int = 25
) -> Dict:
    """
    Manually swap substituents at a specified atom index.
    """
    try:
        parent_mol = Chem.MolFromSmiles(base_smiles)
        if parent_mol is None:
            return {
                "success": False,
                "error": "Invalid base SMILES"
            }
        
        parent_fp = AllChem.GetMorganFingerprintAsBitVect(parent_mol, radius=2, nBits=2048)
        analogs = []
        
        for name, sub_smiles in SUBSTITUENTS.items():
            try:
                sub_mol = Chem.MolFromSmiles(sub_smiles)
                if sub_mol:
                    fp = AllChem.GetMorganFingerprintAsBitVect(sub_mol, radius=2, nBits=2048)
                    sim = DataStructs.TanimotoSimilarity(parent_fp, fp)
                    
                    analogs.append({
                        'variant': name,
                        'smiles': sub_smiles,
                        'tanimoto': round(sim, 4),
                        'mw': Descriptors.MolWt(sub_mol)
                    })
            except Exception:
                continue
        
        analogs = sorted(analogs, key=lambda x: x['tanimoto'], reverse=True)
        
        return {
            "success": True,
            "parent_smiles": base_smiles,
            "analogs": analogs[:n],
            "total_generated": len(analogs)
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


# ============================================================================
# WORKFLOW 5: Full Integrated Analog Pipeline
# ============================================================================

def full_analog_pipeline(
    input_smiles: str,
    threshold: float = 0.70,
    max_hits: int = 25
) -> Dict:
    """
    Full workflow: manual SMILES → confirm → generate → PubChem screen → patent check
    """
    try:
        print("=" * 60)
        print(f"INPUT SMILES: {input_smiles}")
        
        # Step 1: Canonicalize and confirm
        mol = Chem.MolFromSmiles(input_smiles)
        if mol is None:
            return {
                "success": False,
                "error": "Invalid SMILES string entered."
            }
        
        canonical = Chem.MolToSmiles(mol)
        print(f"Canonical SMILES: {canonical}")
        
        # Step 2: Check if it already exists in PubChem
        pubchem_results = pcp.get_compounds(canonical, 'smiles')
        if pubchem_results:
            print(f"Confirmed in PubChem as CID {pubchem_results[0].cid}")
            in_pubchem = True
            parent_cid = pubchem_results[0].cid
        else:
            print("Not found in PubChem — novel compound, proceeding to analog generation.")
            in_pubchem = False
            parent_cid = None
        
        # Step 3: Generate BRICS analogs
        brics_result = generate_brics_analogs(canonical, n=max_hits * 2)
        
        if not brics_result['success']:
            return brics_result
        
        analogs = brics_result['analogs']
        print(f"\nGenerated {len(analogs)} raw analogs")
        
        # Step 4: Filter by similarity threshold
        filtered = [a for a in analogs if a['tanimoto'] >= threshold][:max_hits]
        print(f"Filtered to {len(filtered)} analogs at {threshold:.0%}+ similarity\n")
        
        # Step 5: PubChem screening for each analog
        master_list = []
        for a in filtered:
            time.sleep(0.3)  # Rate limit
            
            try:
                hits = pcp.get_compounds(a['smiles'], 'smiles')
                a['in_pubchem'] = bool(hits)
                a['cid'] = hits[0].cid if hits else None
                
                if a['cid']:
                    patent_result = check_patent_status(a['cid'])
                    a['patents'] = patent_result.get('patents', [])
                    a['patent_free'] = patent_result.get('patent_free', False)
                else:
                    a['patents'] = []
                    a['patent_free'] = True  # Novel compounds are patent-free
                
                a['flag'] = '✅ CLEAR' if a['patent_free'] else f"⚠️ {len(a['patents'])} patents"
                master_list.append(a)
            
            except Exception as e:
                print(f"Error screening {a['smiles']}: {e}")
                continue
        
        return {
            "success": True,
            "input_smiles": input_smiles,
            "canonical_smiles": canonical,
            "in_pubchem": in_pubchem,
            "parent_cid": parent_cid,
            "master_list": master_list,
            "total_candidates": len(master_list),
            "patent_free_count": sum(1 for a in master_list if a['patent_free'])
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def validate_smiles(smiles: str) -> Dict:
    """
    Validate a SMILES string
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "success": False,
                "valid": False,
                "error": "Invalid SMILES string"
            }
        
        canonical = Chem.MolToSmiles(mol)
        return {
            "success": True,
            "valid": True,
            "canonical_smiles": canonical,
            "mw": Descriptors.MolWt(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds()
        }
    
    except Exception as e:
        return {
            "success": False,
            "valid": False,
            "error": str(e)
        }


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No workflow specified"}))
        sys.exit(1)
    
    workflow = sys.argv[1]
    
    try:
        if workflow == 'confirm_and_fetch_similars':
            name_or_smiles = sys.argv[2]
            threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.70
            max_hits = int(sys.argv[4]) if len(sys.argv) > 4 else 25
            result = confirm_and_fetch_similars(name_or_smiles, threshold, max_hits)
        
        elif workflow == 'check_patent_status':
            cid = int(sys.argv[2])
            result = check_patent_status(cid)
        
        elif workflow == 'screen_and_flag_for_masterlist':
            hits_json = json.loads(sys.argv[2])
            max_hits = int(sys.argv[3]) if len(sys.argv) > 3 else 25
            result = screen_and_flag_for_masterlist(hits_json, max_hits)
        
        elif workflow == 'generate_brics_analogs':
            smiles = sys.argv[2]
            n = int(sys.argv[3]) if len(sys.argv) > 3 else 25
            result = generate_brics_analogs(smiles, n)
        
        elif workflow == 'enumerate_substituent_analogs':
            base_smiles = sys.argv[2]
            attachment_idx = int(sys.argv[3])
            n = int(sys.argv[4]) if len(sys.argv) > 4 else 25
            result = enumerate_substituent_analogs(base_smiles, attachment_idx, n)
        
        elif workflow == 'full_analog_pipeline':
            input_smiles = sys.argv[2]
            threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.70
            max_hits = int(sys.argv[4]) if len(sys.argv) > 4 else 25
            result = full_analog_pipeline(input_smiles, threshold, max_hits)
        
        elif workflow == 'validate_smiles':
            smiles = sys.argv[2]
            result = validate_smiles(smiles)
        
        else:
            result = {"error": f"Unknown workflow: {workflow}"}
        
        print(json.dumps(result))
    
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)


if __name__ == '__main__':
    main()
