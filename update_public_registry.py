"""
Update Public Analog Discovery Registry
Generates comprehensive public registry from master discovery log
"""

import json
from datetime import datetime

def generate_public_registry():
    """Generate public registry markdown from master log"""
    
    # Load master discovery log
    with open('/home/ubuntu/pharmasight-latest/MASTER_ANALOG_DISCOVERIES.json', 'r') as f:
        master_data = json.load(f)
    
    # Start building the registry document
    registry_md = f"""# PharmaSight™ Public Analog Discovery Registry

## Purpose
This document serves as a **publicly accessible registry** of pharmaceutical compound analogs discovered and compiled through the PharmaSight™ platform. Publication of this registry establishes **prior art** and **intellectual property rights** for the discoverer.

## Legal Notice
All compounds listed in this registry are claimed as discoveries by the PharmaSight™ research team. This public disclosure establishes prior art and may be used to support patent applications or defend against third-party patent claims.

**Registry Established:** {master_data.get('created', datetime.now().isoformat())[:10]}  
**Last Updated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}  
**Total Discovery Sessions:** {master_data['total_sessions']}  
**Total Analog Discoveries:** {master_data['total_analogs']}

---

"""
    
    # Add each discovery session
    for idx, session in enumerate(master_data['discovery_sessions'], 1):
        parent = session['parent_compound']
        parent_smiles = session.get('parent_smiles', 'N/A')
        discovery_date = session['discovery_timestamp'][:19]
        method = session.get('discovery_method', 'RDKit-based Computational Generation')
        total_analogs = session['total_analogs']
        
        registry_md += f"""## Discovery Session {idx}: {parent} Analogs

**Parent Compound:** {parent}  
**Parent SMILES:** `{parent_smiles}`  
**Discovery Date:** {discovery_date}  
**Discovery Method:** {method}  
**Total Analogs:** {total_analogs}

"""
        
        # Add each analog in the session
        for aidx, analog in enumerate(session['analogs'], 1):
            registry_md += f"""### Analog {idx}.{aidx}: {analog['name']}
- **SMILES:** `{analog['smiles']}`
- **Transformation:** {analog['transformation_applied']}
- **Molecular Weight:** {analog['molecular_weight']} g/mol
- **LogP:** {analog['logp']}
- **TPSA:** {analog['tpsa']} Ų
- **H-Bond Donors:** {analog['h_bond_donors']}
- **H-Bond Acceptors:** {analog['h_bond_acceptors']}
- **Rotatable Bonds:** {analog['rotatable_bonds']}
- **Drug Likeness:** {analog['drug_likeness']}%
- **Therapeutic Potential:** {analog['therapeutic_potential']}
- **Patent Status:** {analog['patent_status']}
- **IP Opportunity Score:** {analog['patent_opportunity_score']}/100
- **Safety Score:** {analog['safety_score']}%
- **Efficacy Score:** {analog['efficacy_score']}%
- **Novelty Score:** {analog['novelty_score']}%
- **Estimated Value:** {analog['estimated_value']}
- **Similarity to Parent:** {round(analog.get('similarity_to_parent', 0) * 100, 1)}%
"""
            
            # Add receptor binding if available
            if analog.get('receptor_binding'):
                registry_md += f"- **Receptor Binding:** "
                bindings = [f"{k}: {v}%" for k, v in analog['receptor_binding'].items()]
                registry_md += ", ".join(bindings) + "\n"
            
            registry_md += f"- **Discovery Timestamp:** {analog.get('discovery_timestamp', datetime.now().isoformat())[:19]}\n\n"
        
        registry_md += "---\n\n"
    
    # Add summary statistics
    registry_md += f"""## Summary Statistics

### Overall Discovery Metrics
- **Total Discovery Sessions:** {master_data['total_sessions']}
- **Total Analog Discoveries:** {master_data['total_analogs']}
- **Parent Compounds Explored:** {', '.join([s['parent_compound'] for s in master_data['discovery_sessions']])}
- **Registry Created:** {master_data.get('created', datetime.now().isoformat())[:10]}
- **Last Updated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}

### Patent Status Distribution
"""
    
    # Calculate patent status distribution
    patent_free = 0
    patent_expired = 0
    patent_pending = 0
    patented = 0
    
    for session in master_data['discovery_sessions']:
        for analog in session['analogs']:
            status = analog['patent_status']
            if status == "Patent-Free":
                patent_free += 1
            elif status == "Patent Expired":
                patent_expired += 1
            elif status == "Patent Pending":
                patent_pending += 1
            else:
                patented += 1
    
    total = master_data['total_analogs']
    registry_md += f"""- **Patent-Free:** {patent_free} ({round(patent_free/total*100, 1)}%)
- **Patent Expired:** {patent_expired} ({round(patent_expired/total*100, 1)}%)
- **Patent Pending:** {patent_pending} ({round(patent_pending/total*100, 1)}%)
- **Patented:** {patented} ({round(patented/total*100, 1)}%)

### High IP Opportunity Analogs
"""
    
    # Count high IP opportunity analogs
    high_ip = sum(1 for s in master_data['discovery_sessions'] 
                  for a in s['analogs'] 
                  if a['patent_opportunity_score'] >= 90)
    
    registry_md += f"""- **IP Opportunity ≥ 90:** {high_ip} analogs ({round(high_ip/total*100, 1)}%)

### Therapeutic Potential Distribution
"""
    
    # Calculate therapeutic potential distribution
    very_high = sum(1 for s in master_data['discovery_sessions'] 
                    for a in s['analogs'] 
                    if a['therapeutic_potential'] == "Very High")
    high = sum(1 for s in master_data['discovery_sessions'] 
               for a in s['analogs'] 
               if a['therapeutic_potential'] == "High")
    moderate = sum(1 for s in master_data['discovery_sessions'] 
                   for a in s['analogs'] 
                   if a['therapeutic_potential'] == "Moderate")
    
    registry_md += f"""- **Very High:** {very_high} analogs ({round(very_high/total*100, 1)}%)
- **High:** {high} analogs ({round(high/total*100, 1)}%)
- **Moderate:** {moderate} analogs ({round(moderate/total*100, 1)}%)

---

## How to Cite This Registry

When referencing compounds from this registry in publications or patent applications, please use:

```
PharmaSight™ Public Analog Discovery Registry. [Compound Name], SMILES: [SMILES String], 
Discovered: [Discovery Date]. Available at: https://github.com/justincihi/pharmasight-platform
```

---

## Contact & IP Inquiries

For licensing inquiries, collaboration opportunities, or IP-related questions regarding compounds in this registry, please contact the PharmaSight™ research team through the GitHub repository.

---

**Document Version:** 2.0  
**Registry Maintained By:** PharmaSight™ Platform  
**GitHub Repository:** https://github.com/justincihi/pharmasight-platform  
**Branch:** analog-discoveries-ip-protected  
**Last Registry Update:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}

---

*This registry is automatically updated as new analog discoveries are made through the PharmaSight™ platform. All compounds listed herein are claimed as discoveries and are published for the purpose of establishing prior art and intellectual property rights.*
"""
    
    # Save the registry
    with open('/home/ubuntu/pharmasight-latest/PUBLIC_ANALOG_DISCOVERY_REGISTRY.md', 'w') as f:
        f.write(registry_md)
    
    return {
        "total_sessions": master_data['total_sessions'],
        "total_analogs": master_data['total_analogs'],
        "patent_free": patent_free,
        "high_ip": high_ip
    }

if __name__ == "__main__":
    print("Updating Public Analog Discovery Registry...")
    print("=" * 70)
    
    result = generate_public_registry()
    
    print(f"\n✅ Public registry updated successfully!")
    print(f"\n📊 Registry Statistics:")
    print(f"   Total Sessions: {result['total_sessions']}")
    print(f"   Total Analogs: {result['total_analogs']}")
    print(f"   Patent-Free Analogs: {result['patent_free']}")
    print(f"   High IP Opportunity (≥90): {result['high_ip']}")
    print(f"\n📁 Updated: PUBLIC_ANALOG_DISCOVERY_REGISTRY.md")

