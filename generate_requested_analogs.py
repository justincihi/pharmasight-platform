"""
Generate Analogs for Kava Lactones, MDAI, Mescaline HCl, and Muscimol
Uses RDKit to create novel analogs with IP protection
"""

import sys
sys.path.append('/home/ubuntu/pharmasight-latest/src')

from rdkit_analog_generator import generate_novel_analogs
import json
from datetime import datetime

# Parent compounds with SMILES
PARENT_COMPOUNDS = {
    "Kavain": {
        "smiles": "COc1ccc(\\C=C\\C(=O)N2CCCCC2)cc1",
        "class": "Kava Lactone",
        "therapeutic_area": "Anxiolytic/Sedative"
    },
    "Yangonin": {
        "smiles": "COc1cc(OC)c2c(c1)C=CC(=O)O2",
        "class": "Kava Lactone",
        "therapeutic_area": "Anxiolytic/Sedative"
    },
    "Methysticin": {
        "smiles": "COc1cc2c(cc1OC)C=CC(=O)O2",
        "class": "Kava Lactone",
        "therapeutic_area": "Anxiolytic/Sedative"
    },
    "MDAI": {
        "smiles": "CC(Cc1ccc2ncccc2c1)N",
        "class": "Entactogen",
        "therapeutic_area": "PTSD/Depression"
    },
    "Mescaline HCl": {
        "smiles": "COc1cc(CCN)cc(OC)c1OC",
        "class": "Psychedelic",
        "therapeutic_area": "Depression/Addiction"
    },
    "Muscimol": {
        "smiles": "C1=C(C(=NO1)N)CO",
        "class": "GABA-A Agonist",
        "therapeutic_area": "Anxiolytic/Sleep"
    }
}

def generate_all_requested_analogs():
    """Generate analogs for all requested compounds"""
    
    all_discoveries = {}
    
    print("=" * 80)
    print("GENERATING ANALOGS FOR REQUESTED COMPOUNDS")
    print("=" * 80)
    
    for compound_name, data in PARENT_COMPOUNDS.items():
        print(f"\n🔬 Generating analogs for {compound_name}...")
        print(f"   Parent SMILES: {data['smiles']}")
        print(f"   Class: {data['class']}")
        print(f"   Therapeutic Area: {data['therapeutic_area']}")
        
        # Generate analogs
        analogs = generate_novel_analogs(
            parent_smiles=data['smiles'],
            parent_name=compound_name,
            num_analogs=15
        )
        
        if analogs and len(analogs) > 0:
            print(f"   ✅ Generated {len(analogs)} novel analogs")
            
            # Add therapeutic area to each analog
            for analog in analogs:
                analog['parent_class'] = data['class']
                analog['therapeutic_area'] = data['therapeutic_area']
            
            all_discoveries[compound_name] = {
                "parent_compound": compound_name,
                "parent_smiles": data['smiles'],
                "parent_class": data['class'],
                "therapeutic_area": data['therapeutic_area'],
                "discovery_timestamp": datetime.now().isoformat(),
                "total_analogs": len(analogs),
                "analogs": analogs
            }
            
            # Show top 3 analogs
            print(f"\n   📊 Top 3 Analogs:")
            for i, analog in enumerate(analogs[:3], 1):
                print(f"      {i}. {analog['name']}")
                print(f"         SMILES: {analog['smiles']}")
                print(f"         Drug Likeness: {analog['drug_likeness']}%")
                print(f"         IP Opportunity: {analog['patent_opportunity_score']}/100")
        else:
            print(f"   ❌ Error: {result.get('error', 'Unknown error')}")
    
    return all_discoveries

def save_discoveries_to_files(discoveries):
    """Save discoveries to individual and master files"""
    
    # Save individual discovery files
    for compound_name, data in discoveries.items():
        filename = f"/home/ubuntu/pharmasight-latest/{compound_name.upper().replace(' ', '_')}_ANALOG_DISCOVERIES.json"
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"\n📁 Saved: {filename}")
    
    # Update master discovery log
    master_log_path = '/home/ubuntu/pharmasight-latest/MASTER_ANALOG_DISCOVERIES.json'
    try:
        with open(master_log_path, 'r') as f:
            master_data = json.load(f)
    except FileNotFoundError:
        master_data = {
            "registry_name": "PharmaSight™ Master Analog Discovery Registry",
            "created": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "total_sessions": 0,
            "total_analogs": 0,
            "discovery_sessions": []
        }
    
    # Add new discovery sessions
    for compound_name, data in discoveries.items():
        session = {
            "session_id": f"SESSION-{compound_name.upper().replace(' ', '-')}-{datetime.now().strftime('%Y%m%d')}",
            "parent_compound": compound_name,
            "parent_smiles": data['parent_smiles'],
            "parent_class": data['parent_class'],
            "therapeutic_area": data['therapeutic_area'],
            "discovery_timestamp": data['discovery_timestamp'],
            "discovery_method": "RDKit-based Computational Generation",
            "total_analogs": data['total_analogs'],
            "analogs": data['analogs']
        }
        master_data['discovery_sessions'].append(session)
    
    # Update metadata
    master_data['total_sessions'] = len(master_data['discovery_sessions'])
    master_data['total_analogs'] = sum(s['total_analogs'] for s in master_data['discovery_sessions'])
    master_data['last_updated'] = datetime.now().isoformat()
    
    # Save master log
    with open(master_log_path, 'w') as f:
        json.dump(master_data, f, indent=2)
    
    print(f"\n📁 Updated: {master_log_path}")
    print(f"   Total Sessions: {master_data['total_sessions']}")
    print(f"   Total Analogs: {master_data['total_analogs']}")
    
    return master_data

if __name__ == "__main__":
    print("\n🧪 PharmaSight™ Analog Generation System")
    print("=" * 80)
    
    # Generate all analogs
    discoveries = generate_all_requested_analogs()
    
    # Save to files
    print("\n" + "=" * 80)
    print("SAVING DISCOVERIES TO FILES")
    print("=" * 80)
    master_data = save_discoveries_to_files(discoveries)
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION SUMMARY")
    print("=" * 80)
    print(f"\n✅ Successfully generated analogs for {len(discoveries)} compounds:")
    for compound_name, data in discoveries.items():
        print(f"   • {compound_name}: {data['total_analogs']} novel analogs")
    
    print(f"\n📊 Master Registry Statistics:")
    print(f"   Total Discovery Sessions: {master_data['total_sessions']}")
    print(f"   Total Analog Discoveries: {master_data['total_analogs']}")
    
    # Count patent-free analogs
    patent_free = sum(1 for s in master_data['discovery_sessions'] 
                      for a in s['analogs'] 
                      if a.get('patent_status') == 'Patent-Free (Novel)')
    print(f"   Patent-Free Analogs: {patent_free}")
    
    print("\n✅ All discoveries saved and ready for IP protection!")
    print("=" * 80)

