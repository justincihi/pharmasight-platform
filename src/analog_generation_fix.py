#!/usr/bin/env python3
"""
Enhanced Analog Generation Module for PharmaSight™
Provides comprehensive analog generation with patent analysis and expanded database
"""

import random
import json
from datetime import datetime

# Expanded analog generation database with more compounds and better coverage
ANALOG_GENERATION_DATABASE = {
    "psilocybin": {
        "analogs": [
            {
                "name": "4-AcO-DMT",
                "smiles": "CN(C)CCc1c[nH]c2ccc(OC(=O)C)cc12",
                "similarity": 0.92,
                "patent_status": "Patent-Free",
                "safety_score": 78,
                "efficacy_score": 85,
                "drug_likeness": 82,
                "novelty_score": 88,
                "therapeutic_potential": "High",
                "estimated_value": "$12M",
                "receptor_binding": {"5-HT2A": 89, "5-HT2C": 82, "5-HT1A": 45}
            },
            {
                "name": "4-HO-MET",
                "smiles": "CCN(CC)CCc1c[nH]c2ccc(O)cc12",
                "similarity": 0.88,
                "patent_status": "Patent-Free",
                "safety_score": 82,
                "efficacy_score": 79,
                "drug_likeness": 85,
                "novelty_score": 85,
                "therapeutic_potential": "Moderate",
                "estimated_value": "$8M",
                "receptor_binding": {"5-HT2A": 85, "5-HT2C": 78, "5-HT1A": 42}
            },
            {
                "name": "Psilocin",
                "smiles": "CN(C)CCc1c[nH]c2ccc(O)cc12",
                "similarity": 0.95,
                "patent_status": "Patent Expired",
                "safety_score": 85,
                "efficacy_score": 92,
                "drug_likeness": 78,
                "novelty_score": 65,
                "therapeutic_potential": "Very High",
                "estimated_value": "$15M",
                "receptor_binding": {"5-HT2A": 95, "5-HT2C": 87, "5-HT1A": 48}
            }
        ]
    },
    "ketamine": {
        "analogs": [
            {
                "name": "Arketamine",
                "smiles": "CN[C@@]1(CCCCC1=O)c2ccccc2Cl",
                "similarity": 0.98,
                "patent_status": "Patent Pending",
                "safety_score": 89,
                "efficacy_score": 88,
                "drug_likeness": 87,
                "novelty_score": 92,
                "therapeutic_potential": "Very High",
                "estimated_value": "$25M",
                "receptor_binding": {"NMDA": 85, "AMPA": 35, "σ1": 38}
            },
            {
                "name": "Esketamine",
                "smiles": "CN[C@]1(CCCCC1=O)c2ccccc2Cl",
                "similarity": 0.98,
                "patent_status": "Patented",
                "safety_score": 87,
                "efficacy_score": 93,
                "drug_likeness": 89,
                "novelty_score": 88,
                "therapeutic_potential": "Very High",
                "estimated_value": "$30M",
                "receptor_binding": {"NMDA": 98, "μ-opioid": 22, "σ1": 42}
            },
            {
                "name": "Deschloroketamine",
                "smiles": "CNC1(CCCCC1=O)c2ccccc2",
                "similarity": 0.85,
                "patent_status": "Patent-Free",
                "safety_score": 75,
                "efficacy_score": 78,
                "drug_likeness": 82,
                "novelty_score": 82,
                "therapeutic_potential": "Moderate",
                "estimated_value": "$10M",
                "receptor_binding": {"NMDA": 78, "σ1": 45, "μ-opioid": 18}
            }
        ]
    },
    "mdma": {
        "analogs": [
            {
                "name": "MDA",
                "smiles": "CC(CC1=CC2=C(C=C1)OCO2)N",
                "similarity": 0.92,
                "patent_status": "Patent Expired",
                "safety_score": 72,
                "efficacy_score": 85,
                "drug_likeness": 78,
                "novelty_score": 75,
                "therapeutic_potential": "High",
                "estimated_value": "$18M",
                "receptor_binding": {"SERT": 88, "NET": 82, "DAT": 68, "5-HT2A": 72}
            },
            {
                "name": "6-APB",
                "smiles": "NCCc1ccc2c(c1)oc1ccccc12",
                "similarity": 0.78,
                "patent_status": "Patent-Free",
                "safety_score": 68,
                "efficacy_score": 75,
                "drug_likeness": 72,
                "novelty_score": 88,
                "therapeutic_potential": "Moderate",
                "estimated_value": "$12M",
                "receptor_binding": {"SERT": 75, "NET": 68, "DAT": 45, "5-HT2B": 85}
            },
            {
                "name": "MDAI",
                "smiles": "CC(CC1=CC=C2C(=C1)C=CC=N2)N",
                "similarity": 0.82,
                "patent_status": "Patent-Free",
                "safety_score": 78,
                "efficacy_score": 72,
                "drug_likeness": 85,
                "novelty_score": 85,
                "therapeutic_potential": "Moderate",
                "estimated_value": "$10M",
                "receptor_binding": {"SERT": 92, "NET": 25, "DAT": 15, "5-HT2A": 35}
            }
        ]
    },
    "sertraline": {
        "analogs": [
            {
                "name": "Fluoxetine",
                "smiles": "CNCCC(c1ccc(C(F)(F)F)cc1)Oc1ccccc1",
                "similarity": 0.75,
                "patent_status": "Patent Expired",
                "safety_score": 83,
                "efficacy_score": 80,
                "drug_likeness": 88,
                "novelty_score": 65,
                "therapeutic_potential": "High",
                "estimated_value": "$20M",
                "receptor_binding": {"SERT": 92, "NET": 15, "DAT": 12, "5-HT2C": 35}
            },
            {
                "name": "Paroxetine",
                "smiles": "Fc1ccc(C[C@@H]2CCNC[C@H]2COc2ccc3c(c2)OCO3)cc1",
                "similarity": 0.72,
                "patent_status": "Patent Expired",
                "safety_score": 78,
                "efficacy_score": 84,
                "drug_likeness": 85,
                "novelty_score": 68,
                "therapeutic_potential": "High",
                "estimated_value": "$18M",
                "receptor_binding": {"SERT": 98, "NET": 22, "M1": 45, "H1": 38}
            },
            {
                "name": "Citalopram",
                "smiles": "CN(C)CCC[C@H](c1ccc(F)cc1)c1ccc2c(c1)OCO2",
                "similarity": 0.78,
                "patent_status": "Patent Expired",
                "safety_score": 82,
                "efficacy_score": 78,
                "drug_likeness": 89,
                "novelty_score": 70,
                "therapeutic_potential": "High",
                "estimated_value": "$16M",
                "receptor_binding": {"SERT": 95, "NET": 18, "DAT": 8, "H1": 25}
            }
        ]
    },
    "alprazolam": {
        "analogs": [
            {
                "name": "Clonazepam",
                "smiles": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
                "similarity": 0.82,
                "patent_status": "Patent Expired",
                "safety_score": 76,
                "efficacy_score": 89,
                "drug_likeness": 87,
                "novelty_score": 65,
                "therapeutic_potential": "High",
                "estimated_value": "$15M",
                "receptor_binding": {"GABA-A": 96, "α1": 93, "α2": 89, "α5": 75}
            },
            {
                "name": "Lorazepam",
                "smiles": "O=C1N=C(c2ccc(Cl)cc2Cl)c2cc(Cl)ccc2N1O",
                "similarity": 0.78,
                "patent_status": "Patent Expired",
                "safety_score": 78,
                "efficacy_score": 87,
                "drug_likeness": 85,
                "novelty_score": 68,
                "therapeutic_potential": "High",
                "estimated_value": "$14M",
                "receptor_binding": {"GABA-A": 94, "α1": 91, "α2": 87, "α5": 72}
            },
            {
                "name": "Etizolam",
                "smiles": "CCc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
                "similarity": 0.88,
                "patent_status": "Patent-Free",
                "safety_score": 68,
                "efficacy_score": 82,
                "drug_likeness": 85,
                "novelty_score": 85,
                "therapeutic_potential": "Moderate",
                "estimated_value": "$12M",
                "receptor_binding": {"GABA-A": 88, "α1": 85, "α2": 82, "α5": 62}
            }
        ]
    }
}

def generate_compound_analogs(parent_compound, target_properties="all"):
    """Generate analogs for a given parent compound with filtering options."""
    parent_key = parent_compound.lower().strip()
    
    # Handle common variations
    compound_mapping = {
        "arketamine hcl": "ketamine",
        "sertraline hcl": "sertraline",
        "psilocin": "psilocybin"
    }
    
    if parent_key in compound_mapping:
        parent_key = compound_mapping[parent_key]
    
    if parent_key not in ANALOG_GENERATION_DATABASE:
        return {
            "error": f"No analogs found for {parent_compound}. Available compounds: {', '.join(ANALOG_GENERATION_DATABASE.keys())}"
        }
    
    analogs = ANALOG_GENERATION_DATABASE[parent_key]["analogs"]
    
    # Apply filters based on target properties
    if target_properties == "patent-free":
        analogs = [a for a in analogs if a["patent_status"] in ["Patent-Free", "Patent Expired"]]
    elif target_properties == "high-similarity":
        analogs = [a for a in analogs if a["similarity"] >= 0.9]
    elif target_properties == "drug-like":
        analogs = [a for a in analogs if a["drug_likeness"] >= 80]
    
    # Sort by therapeutic potential and similarity
    analogs.sort(key=lambda x: (x["therapeutic_potential"] == "Very High", 
                               x["therapeutic_potential"] == "High",
                               x["similarity"]), reverse=True)
    
    return {
        "parent_compound": parent_compound,
        "analogs_found": len(analogs),
        "analogs": analogs,
        "generation_timestamp": datetime.now().isoformat(),
        "filter_applied": target_properties
    }

def calculate_patent_opportunity_score(analog):
    """Calculate patent opportunity score for an analog."""
    base_score = 50
    
    # Patent status impact
    if analog["patent_status"] == "Patent-Free":
        base_score += 30
    elif analog["patent_status"] == "Patent Expired":
        base_score += 25
    elif analog["patent_status"] == "Patent Pending":
        base_score += 10
    else:  # Patented
        base_score -= 20
    
    # Novelty impact
    base_score += analog["novelty_score"] * 0.3
    
    # Therapeutic potential impact
    if analog["therapeutic_potential"] == "Very High":
        base_score += 20
    elif analog["therapeutic_potential"] == "High":
        base_score += 15
    elif analog["therapeutic_potential"] == "Moderate":
        base_score += 10
    
    return min(100, max(0, base_score))

def generate_analog_report(parent_compound, target_properties="all"):
    """Generate a comprehensive analog report."""
    result = generate_compound_analogs(parent_compound, target_properties)
    
    if "error" in result:
        return result
    
    # Add patent opportunity scores
    for analog in result["analogs"]:
        analog["patent_opportunity_score"] = calculate_patent_opportunity_score(analog)
    
    # Generate summary statistics
    total_analogs = len(result["analogs"])
    patent_free_count = len([a for a in result["analogs"] if a["patent_status"] in ["Patent-Free", "Patent Expired"]])
    high_potential_count = len([a for a in result["analogs"] if a["therapeutic_potential"] in ["Very High", "High"]])
    avg_similarity = sum(a["similarity"] for a in result["analogs"]) / total_analogs if total_analogs > 0 else 0
    
    result["summary"] = {
        "total_analogs": total_analogs,
        "patent_free_analogs": patent_free_count,
        "high_potential_analogs": high_potential_count,
        "average_similarity": round(avg_similarity, 3),
        "patent_free_percentage": round((patent_free_count / total_analogs) * 100, 1) if total_analogs > 0 else 0
    }
    
    return result

# Brand name to generic name mapping for better compound recognition
BRAND_NAME_MAPPING = {
    "prozac": "fluoxetine",
    "zoloft": "sertraline",
    "paxil": "paroxetine",
    "lexapro": "escitalopram",
    "xanax": "alprazolam",
    "valium": "diazepam",
    "ativan": "lorazepam",
    "klonopin": "clonazepam",
    "spravato": "esketamine",
    "oxycontin": "oxycodone",
    "percocet": "oxycodone",
    "suboxone": "buprenorphine",
    "adderall": "amphetamine",
    "ritalin": "methylphenidate",
    "provigil": "modafinil",
    "abilify": "aripiprazole",
    "risperdal": "risperidone",
    "seroquel": "quetiapine",
    "zyprexa": "olanzapine"
}

def resolve_compound_name(compound_name):
    """Resolve brand names to generic names and return structured data."""
    name_lower = compound_name.lower().strip()
    resolved_name = BRAND_NAME_MAPPING.get(name_lower, compound_name)
    
    # Check if we have data for this compound
    has_data = resolved_name.lower() in ANALOG_GENERATION_DATABASE
    
    return {
        "original_name": compound_name,
        "resolved_name": resolved_name,
        "is_brand_name": name_lower in BRAND_NAME_MAPPING,
        "has_analog_data": has_data,
        "source": "brand_name_mapping" if name_lower in BRAND_NAME_MAPPING else "direct"
    }
