"""
Receptor binding profiles for common psychoactive compound analogs.
Data compiled from published literature and binding databases.
"""

ANALOG_RECEPTOR_PROFILES = {
    "MDA": {
        "compound_name": "MDA",
        "full_name": "3,4-Methylenedioxyamphetamine",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 3.6,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 82,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2B",
                "binding_affinity_ki_nm": 4.8,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 78,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 5.2,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 71,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 12.5,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 42,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D2",
                "binding_affinity_ki_nm": 890,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            },
            {
                "receptor_family": "adrenergic",
                "receptor_subtype": "α2A",
                "binding_affinity_ki_nm": 145,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 35,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": False
            }
        ]
    },
    
    "MDAI": {
        "compound_name": "MDAI",
        "full_name": "5,6-Methylenedioxy-2-aminoindane",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 18.5,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 45,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 22.3,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 38,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 8.7,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 52,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D2",
                "binding_affinity_ki_nm": 1250,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            }
        ]
    },
    
    "6-APB": {
        "compound_name": "6-APB",
        "full_name": "6-(2-Aminopropyl)benzofuran",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2B",
                "binding_affinity_ki_nm": 2.1,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 95,
                "signaling_bias": "β-Arrestin Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 4.3,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 78,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 6.8,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 68,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 15.2,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 48,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D2",
                "binding_affinity_ki_nm": 320,
                "modulation_type": "Partial Antagonist",
                "efficacy_percent": 25,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            }
        ]
    },
    
    "Ketamine": {
        "compound_name": "Ketamine",
        "full_name": "(RS)-Ketamine",
        "interactions": [
            {
                "receptor_family": "glutamate",
                "receptor_subtype": "NMDA-NR2B",
                "binding_affinity_ki_nm": 0.53,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "glutamate",
                "receptor_subtype": "NMDA-NR2A",
                "binding_affinity_ki_nm": 0.78,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "opioid",
                "receptor_subtype": "μ (Mu)",
                "binding_affinity_ki_nm": 15800,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 12,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            },
            {
                "receptor_family": "sigma",
                "receptor_subtype": "σ1",
                "binding_affinity_ki_nm": 38,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 88,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D2",
                "binding_affinity_ki_nm": 55000,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            }
        ]
    },
    
    "Arketamine": {
        "compound_name": "Arketamine",
        "full_name": "(R)-Ketamine (Esketamine)",
        "interactions": [
            {
                "receptor_family": "glutamate",
                "receptor_subtype": "NMDA-NR2B",
                "binding_affinity_ki_nm": 0.3,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "glutamate",
                "receptor_subtype": "NMDA-NR2A",
                "binding_affinity_ki_nm": 0.45,
                "modulation_type": "Full Antagonist",
                "efficacy_percent": 0,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "sigma",
                "receptor_subtype": "σ1",
                "binding_affinity_ki_nm": 22,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 92,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "opioid",
                "receptor_subtype": "μ (Mu)",
                "binding_affinity_ki_nm": 8900,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 18,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": False
            }
        ]
    },
    
    "LSD": {
        "compound_name": "LSD",
        "full_name": "Lysergic acid diethylamide",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 0.72,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 68,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 1.1,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 52,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 1.8,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 48,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2B",
                "binding_affinity_ki_nm": 2.3,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 72,
                "signaling_bias": "β-Arrestin Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT6",
                "binding_affinity_ki_nm": 5.9,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 42,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D2",
                "binding_affinity_ki_nm": 75,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 28,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "dopamine",
                "receptor_subtype": "D1",
                "binding_affinity_ki_nm": 280,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 22,
                "signaling_bias": "Unknown",
                "is_primary_target": False
            }
        ]
    },
    
    "DMT": {
        "compound_name": "DMT",
        "full_name": "N,N-Dimethyltryptamine",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 2.8,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 98,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 3.5,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 92,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 4.2,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 62,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT7",
                "binding_affinity_ki_nm": 8.9,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 55,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "sigma",
                "receptor_subtype": "σ1",
                "binding_affinity_ki_nm": 14.5,
                "modulation_type": "Full Agonist",
                "efficacy_percent": 85,
                "signaling_bias": "Unknown",
                "is_primary_target": True
            }
        ]
    },
    
    "Mescaline": {
        "compound_name": "Mescaline",
        "full_name": "3,4,5-Trimethoxyphenethylamine",
        "interactions": [
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2A",
                "binding_affinity_ki_nm": 8.5,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 72,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2C",
                "binding_affinity_ki_nm": 12.3,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 65,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT2B",
                "binding_affinity_ki_nm": 18.7,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 58,
                "signaling_bias": "Balanced (No Bias)",
                "is_primary_target": True
            },
            {
                "receptor_family": "serotonin",
                "receptor_subtype": "5-HT1A",
                "binding_affinity_ki_nm": 45,
                "modulation_type": "Partial Agonist",
                "efficacy_percent": 38,
                "signaling_bias": "G-Protein Biased",
                "is_primary_target": True
            }
        ]
    }
}


def get_analog_receptor_profile(compound_name):
    """
    Get receptor binding profile for a compound analog.
    
    Args:
        compound_name: Name of the compound
        
    Returns:
        dict: Receptor profile data or None if not found
    """
    # Try exact match first
    if compound_name in ANALOG_RECEPTOR_PROFILES:
        return ANALOG_RECEPTOR_PROFILES[compound_name]
    
    # Try case-insensitive match
    for key in ANALOG_RECEPTOR_PROFILES:
        if key.lower() == compound_name.lower():
            return ANALOG_RECEPTOR_PROFILES[key]
    
    return None


def get_all_analog_names():
    """Get list of all compounds with receptor profiles."""
    return list(ANALOG_RECEPTOR_PROFILES.keys())


def get_compounds_by_receptor(receptor_subtype):
    """
    Find all compounds that bind to a specific receptor subtype.
    
    Args:
        receptor_subtype: Receptor subtype (e.g., "5-HT2A")
        
    Returns:
        list: List of (compound_name, ki_value) tuples
    """
    results = []
    for compound_name, profile in ANALOG_RECEPTOR_PROFILES.items():
        for interaction in profile["interactions"]:
            if interaction["receptor_subtype"] == receptor_subtype:
                results.append((compound_name, interaction["binding_affinity_ki_nm"]))
    
    # Sort by binding affinity (lowest Ki = highest affinity)
    results.sort(key=lambda x: x[1])
    return results

