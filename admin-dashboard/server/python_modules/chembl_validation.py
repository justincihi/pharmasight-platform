#!/usr/bin/env python3
"""
ChEMBL Validation Module
Compares PharmaSight predictions against experimental binding data from ChEMBL
"""

import requests
import json
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import sqlite3
import os

class ChEMBLValidator:
    """Validates predictions against ChEMBL experimental data"""
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    # Map receptor names to ChEMBL target IDs
    RECEPTOR_TO_CHEMBL = {
        # Serotonin receptors
        "5-HT1A": "CHEMBL214",
        "5-HT2A": "CHEMBL224",
        "5-HT2B": "CHEMBL225",
        "5-HT2C": "CHEMBL226",
        "5-HT3": "CHEMBL1899",
        # Dopamine receptors
        "D1": "CHEMBL2056",
        "D2": "CHEMBL217",
        "D3": "CHEMBL234",
        "D4": "CHEMBL219",
        # GABA receptors
        "GABA-A": "CHEMBL2093872",
        "GABA-A_alpha1": "CHEMBL4296078",
        # Adrenergic receptors
        "Alpha1A": "CHEMBL229",
        "Alpha2A": "CHEMBL1867",
        "Beta1": "CHEMBL213",
        "Beta2": "CHEMBL210",
        # Muscarinic receptors
        "M1": "CHEMBL216",
        "M2": "CHEMBL211",
        "M3": "CHEMBL245",
        # Opioid receptors
        "MOR": "CHEMBL233",
        "DOR": "CHEMBL236",
        "KOR": "CHEMBL237",
        # Histamine receptors
        "H1": "CHEMBL231",
        "H2": "CHEMBL1782361",
        "H3": "CHEMBL264",
        # Glutamate receptors
        "NMDA": "CHEMBL1907601",
        # Cannabinoid receptors
        "CB1": "CHEMBL218",
        "CB2": "CHEMBL253",
    }
    
    def __init__(self, cache_db: str = "chembl_cache.db"):
        self.cache_db = cache_db
        self._init_cache()
        
    def _init_cache(self):
        """Initialize SQLite cache for ChEMBL data"""
        conn = sqlite3.connect(self.cache_db)
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS chembl_cache (
                smiles TEXT,
                receptor TEXT,
                experimental_ki REAL,
                experimental_ic50 REAL,
                activity_type TEXT,
                assay_description TEXT,
                chembl_id TEXT,
                timestamp TEXT,
                PRIMARY KEY (smiles, receptor)
            )
        """)
        conn.commit()
        conn.close()
        
    def _get_cached_data(self, smiles: str, receptor: str) -> Optional[Dict]:
        """Get cached experimental data"""
        conn = sqlite3.connect(self.cache_db)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT * FROM chembl_cache WHERE smiles = ? AND receptor = ?",
            (smiles, receptor)
        )
        row = cursor.fetchone()
        conn.close()
        
        if row:
            return {
                "smiles": row[0],
                "receptor": row[1],
                "experimental_ki": row[2],
                "experimental_ic50": row[3],
                "activity_type": row[4],
                "assay_description": row[5],
                "chembl_id": row[6],
                "timestamp": row[7]
            }
        return None
        
    def _cache_data(self, data: Dict):
        """Cache experimental data"""
        conn = sqlite3.connect(self.cache_db)
        cursor = conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO chembl_cache 
            (smiles, receptor, experimental_ki, experimental_ic50, activity_type, 
             assay_description, chembl_id, timestamp)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            data.get("smiles"),
            data.get("receptor"),
            data.get("experimental_ki"),
            data.get("experimental_ic50"),
            data.get("activity_type"),
            data.get("assay_description"),
            data.get("chembl_id"),
            datetime.now().isoformat()
        ))
        conn.commit()
        conn.close()
    
    def search_compound(self, smiles: str) -> Optional[str]:
        """Search for a compound in ChEMBL by SMILES, return ChEMBL ID"""
        try:
            url = f"{self.BASE_URL}/molecule/search"
            params = {
                "q": smiles,
                "format": "json"
            }
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data.get("molecules"):
                    return data["molecules"][0].get("molecule_chembl_id")
        except Exception as e:
            print(f"ChEMBL search error: {e}")
        return None
    
    def get_compound_by_smiles(self, smiles: str) -> Optional[Dict]:
        """Get compound info from ChEMBL by SMILES"""
        try:
            # Use similarity search with exact match
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            canonical_smiles = Chem.MolToSmiles(mol)
            
            url = f"{self.BASE_URL}/molecule.json"
            params = {
                "molecule_structures__canonical_smiles__flexmatch": canonical_smiles,
                "limit": 1
            }
            response = requests.get(url, params=params, timeout=15)
            if response.status_code == 200:
                data = response.json()
                if data.get("molecules"):
                    return data["molecules"][0]
        except Exception as e:
            print(f"ChEMBL compound lookup error: {e}")
        return None
    
    def get_activities_for_compound(self, chembl_id: str, receptor: str = None) -> List[Dict]:
        """Get bioactivity data for a compound"""
        activities = []
        try:
            url = f"{self.BASE_URL}/activity.json"
            params = {
                "molecule_chembl_id": chembl_id,
                "pchembl_value__isnull": False,
                "limit": 100,
                "format": "json"
            }
            
            # Filter by target if receptor specified
            if receptor and receptor in self.RECEPTOR_TO_CHEMBL:
                params["target_chembl_id"] = self.RECEPTOR_TO_CHEMBL[receptor]
            
            response = requests.get(url, params=params, timeout=15)
            if response.status_code == 200:
                data = response.json()
                for activity in data.get("activities", []):
                    act_data = {
                        "target_chembl_id": activity.get("target_chembl_id"),
                        "target_name": activity.get("target_pref_name"),
                        "activity_type": activity.get("standard_type"),
                        "activity_value": activity.get("standard_value"),
                        "activity_units": activity.get("standard_units"),
                        "pchembl_value": activity.get("pchembl_value"),
                        "assay_description": activity.get("assay_description"),
                        "assay_type": activity.get("assay_type")
                    }
                    activities.append(act_data)
        except Exception as e:
            print(f"ChEMBL activity lookup error: {e}")
            
        return activities
    
    def validate_prediction(self, smiles: str, receptor: str, predicted_ki: str, 
                          predicted_score: float) -> Dict:
        """Validate a single prediction against ChEMBL data"""
        
        # Check cache first
        cached = self._get_cached_data(smiles, receptor)
        if cached:
            return self._compare_prediction(cached, predicted_ki, predicted_score)
        
        # Look up compound
        compound = self.get_compound_by_smiles(smiles)
        if not compound:
            return {
                "status": "not_found",
                "message": "Compound not found in ChEMBL",
                "experimental_data": None,
                "prediction_accuracy": None
            }
        
        chembl_id = compound.get("molecule_chembl_id")
        
        # Get activities for this receptor
        activities = self.get_activities_for_compound(chembl_id, receptor)
        
        if not activities:
            return {
                "status": "no_data",
                "message": f"No experimental data for {receptor} in ChEMBL",
                "chembl_id": chembl_id,
                "experimental_data": None,
                "prediction_accuracy": None
            }
        
        # Find Ki or IC50 data
        ki_data = None
        ic50_data = None
        
        for activity in activities:
            if activity["activity_type"] == "Ki":
                ki_data = activity
            elif activity["activity_type"] == "IC50" and not ki_data:
                ic50_data = activity
        
        experimental = ki_data or ic50_data
        if not experimental:
            return {
                "status": "no_binding_data",
                "message": "No Ki/IC50 data available",
                "chembl_id": chembl_id,
                "experimental_data": activities[:3],
                "prediction_accuracy": None
            }
        
        # Cache the data
        cache_entry = {
            "smiles": smiles,
            "receptor": receptor,
            "experimental_ki": float(experimental["activity_value"]) if experimental["activity_type"] == "Ki" else None,
            "experimental_ic50": float(experimental["activity_value"]) if experimental["activity_type"] == "IC50" else None,
            "activity_type": experimental["activity_type"],
            "assay_description": experimental.get("assay_description", ""),
            "chembl_id": chembl_id
        }
        self._cache_data(cache_entry)
        
        return self._compare_prediction(cache_entry, predicted_ki, predicted_score)
    
    def _compare_prediction(self, experimental: Dict, predicted_ki: str, 
                           predicted_score: float) -> Dict:
        """Compare predicted vs experimental values"""
        
        exp_value = experimental.get("experimental_ki") or experimental.get("experimental_ic50")
        exp_type = "Ki" if experimental.get("experimental_ki") else "IC50"
        
        if not exp_value:
            return {
                "status": "no_value",
                "experimental_data": experimental,
                "prediction_accuracy": None
            }
        
        # Parse predicted Ki range
        pred_range = self._parse_ki_range(predicted_ki)
        
        # Calculate accuracy
        accuracy = self._calculate_accuracy(exp_value, pred_range, predicted_score)
        
        return {
            "status": "validated",
            "experimental_value_nM": exp_value,
            "experimental_type": exp_type,
            "predicted_range": predicted_ki,
            "prediction_accuracy": accuracy,
            "accuracy_category": self._categorize_accuracy(accuracy),
            "chembl_id": experimental.get("chembl_id"),
            "assay_description": experimental.get("assay_description", "")[:100]
        }
    
    def _parse_ki_range(self, ki_str: str) -> Tuple[float, float]:
        """Parse Ki range string to min/max nM values"""
        ki_str = ki_str.strip()
        
        if ki_str.startswith("<"):
            val = float(ki_str[1:].replace("nM", "").replace("µM", "").strip())
            if "µM" in ki_str:
                val *= 1000
            return (0, val)
        elif ki_str.startswith(">"):
            val = float(ki_str[1:].replace("nM", "").replace("µM", "").strip())
            if "µM" in ki_str:
                val *= 1000
            return (val, val * 10)
        elif "-" in ki_str:
            parts = ki_str.replace("nM", "").replace("µM", "").split("-")
            low = float(parts[0].strip())
            high = float(parts[1].strip())
            if "µM" in ki_str:
                low *= 1000
                high *= 1000
            return (low, high)
        else:
            try:
                val = float(ki_str.replace("nM", "").replace("µM", "").strip())
                if "µM" in ki_str:
                    val *= 1000
                return (val * 0.5, val * 2)
            except:
                return (100, 10000)  # Default range
    
    def _calculate_accuracy(self, experimental: float, predicted_range: Tuple[float, float],
                           score: float) -> float:
        """Calculate prediction accuracy (0-1)"""
        pred_low, pred_high = predicted_range
        
        # Check if experimental falls within predicted range
        if pred_low <= experimental <= pred_high:
            return 1.0
        
        # Calculate how far off we are (log scale)
        import math
        pred_mid = (pred_low + pred_high) / 2
        
        log_diff = abs(math.log10(experimental) - math.log10(pred_mid))
        
        # Accuracy decreases with log difference
        # Within 1 log unit = good, 2 log units = moderate, >2 = poor
        accuracy = max(0, 1 - (log_diff / 3))
        
        return round(accuracy, 3)
    
    def _categorize_accuracy(self, accuracy: float) -> str:
        """Categorize accuracy score"""
        if accuracy >= 0.9:
            return "Excellent"
        elif accuracy >= 0.7:
            return "Good"
        elif accuracy >= 0.5:
            return "Moderate"
        elif accuracy >= 0.3:
            return "Fair"
        else:
            return "Poor"
    
    def validate_screening_results(self, screening_result: Dict) -> Dict:
        """Validate all hits from a screening result"""
        smiles = screening_result.get("smiles")
        hits = screening_result.get("receptor_hits", [])
        
        validations = []
        validated_count = 0
        accuracy_sum = 0
        
        for hit in hits[:10]:  # Limit to top 10 hits
            receptor = hit.get("receptor")
            if receptor in self.RECEPTOR_TO_CHEMBL:
                validation = self.validate_prediction(
                    smiles,
                    receptor,
                    hit.get("predicted_ki", ""),
                    hit.get("binding_score", 0)
                )
                validation["receptor"] = receptor
                validations.append(validation)
                
                if validation.get("prediction_accuracy") is not None:
                    validated_count += 1
                    accuracy_sum += validation["prediction_accuracy"]
        
        overall_accuracy = accuracy_sum / validated_count if validated_count > 0 else None
        
        return {
            "smiles": smiles,
            "total_hits_checked": len(validations),
            "validated_hits": validated_count,
            "overall_accuracy": round(overall_accuracy, 3) if overall_accuracy else None,
            "accuracy_category": self._categorize_accuracy(overall_accuracy) if overall_accuracy else "No data",
            "validations": validations
        }

    def get_known_ligands(self, receptor: str, limit: int = 10) -> List[Dict]:
        """Get known ligands for a receptor from ChEMBL"""
        if receptor not in self.RECEPTOR_TO_CHEMBL:
            return []
        
        target_id = self.RECEPTOR_TO_CHEMBL[receptor]
        
        try:
            url = f"{self.BASE_URL}/activity.json"
            params = {
                "target_chembl_id": target_id,
                "pchembl_value__gte": 7,  # High affinity (Ki < 100 nM)
                "limit": limit,
                "format": "json"
            }
            
            response = requests.get(url, params=params, timeout=15)
            if response.status_code == 200:
                data = response.json()
                ligands = []
                seen_molecules = set()
                
                for activity in data.get("activities", []):
                    mol_id = activity.get("molecule_chembl_id")
                    if mol_id not in seen_molecules:
                        seen_molecules.add(mol_id)
                        ligands.append({
                            "chembl_id": mol_id,
                            "name": activity.get("molecule_pref_name") or mol_id,
                            "activity_type": activity.get("standard_type"),
                            "activity_value": activity.get("standard_value"),
                            "pchembl": activity.get("pchembl_value"),
                            "smiles": activity.get("canonical_smiles")
                        })
                
                return ligands
        except Exception as e:
            print(f"ChEMBL ligand lookup error: {e}")
        
        return []


# Helper function for API endpoints
def validate_compound(smiles: str, screening_result: Dict = None) -> Dict:
    """Validate a compound's predictions against ChEMBL"""
    validator = ChEMBLValidator()
    
    if screening_result:
        return validator.validate_screening_results(screening_result)
    else:
        return {
            "status": "error",
            "message": "Screening result required for validation"
        }
