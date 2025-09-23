"""
Populate PharmaSight™ with sample data for demonstration
Creates sample discoveries, DOI references, and research topics
"""

from discovery_logging_system import get_discovery_logger
from confidence_doi_system import get_confidence_doi_manager
from admin_research_manager import get_research_topic_manager
import datetime

def populate_sample_data():
    """Populate the system with comprehensive sample data."""
    
    print("Populating PharmaSight™ with sample data...")
    
    # Get managers
    discovery_logger = get_discovery_logger()
    confidence_manager = get_confidence_doi_manager()
    research_manager = get_research_topic_manager()
    
    # 1. Add sample research topics
    print("Adding research topics...")
    
    topics = [
        {
            'title': 'Novel Psychedelic Therapeutics for Treatment-Resistant Depression',
            'description': 'Research into psilocybin and LSD analogs with improved safety profiles',
            'keywords': ['psychedelics', 'depression', 'serotonin', '5-HT2A', 'neuroplasticity'],
            'priority': 5,
            'confidence_threshold': 0.85,
            'ip_focus': True,
            'therapeutic_areas': ['Mental Health', 'Neurology'],
            'target_compounds': ['psilocybin', 'LSD', 'DMT']
        },
        {
            'title': 'Safer Opioid Alternatives with Reduced Addiction Potential',
            'description': 'Development of μ-opioid receptor modulators with biased signaling',
            'keywords': ['opioids', 'addiction', 'pain management', 'biased signaling', 'G-protein'],
            'priority': 5,
            'confidence_threshold': 0.90,
            'ip_focus': True,
            'therapeutic_areas': ['Pain Management', 'Addiction Medicine'],
            'target_compounds': ['morphine', 'fentanyl', 'tramadol']
        },
        {
            'title': 'Next-Generation Anxiolytics Beyond Benzodiazepines',
            'description': 'GABA-A receptor modulators with improved selectivity and reduced tolerance',
            'keywords': ['anxiety', 'GABA', 'benzodiazepines', 'allosteric modulation', 'tolerance'],
            'priority': 4,
            'confidence_threshold': 0.80,
            'ip_focus': True,
            'therapeutic_areas': ['Mental Health', 'Neurology'],
            'target_compounds': ['alprazolam', 'diazepam', 'lorazepam']
        },
        {
            'title': 'Cognitive Enhancement Compounds for Neurodegenerative Diseases',
            'description': 'Nootropics targeting acetylcholine, dopamine, and glutamate systems',
            'keywords': ['nootropics', 'cognition', 'Alzheimer', 'acetylcholine', 'NMDA'],
            'priority': 4,
            'confidence_threshold': 0.75,
            'ip_focus': False,
            'therapeutic_areas': ['Neurology', 'Geriatrics'],
            'target_compounds': ['modafinil', 'piracetam', 'donepezil']
        },
        {
            'title': 'Empathogenic Compounds for PTSD and Social Anxiety',
            'description': 'MDMA analogs with enhanced therapeutic index and reduced neurotoxicity',
            'keywords': ['MDMA', 'PTSD', 'empathogen', 'serotonin', 'oxytocin', 'social bonding'],
            'priority': 3,
            'confidence_threshold': 0.80,
            'ip_focus': True,
            'therapeutic_areas': ['Mental Health', 'Trauma Therapy'],
            'target_compounds': ['MDMA', 'MDA', '6-APB']
        }
    ]
    
    for topic in topics:
        research_manager.create_topic(**topic)
    
    # 2. Add sample compound discoveries
    print("Adding compound discoveries...")
    
    discoveries = [
        {
            'compound_name': 'PSI-2024-A1',
            'smiles': 'CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12',
            'category': 'Psychedelic-2024-A1',
            'confidence_score': 0.95,
            'ip_status': 'High Opportunity - No existing patents',
            'source_research_topic_id': 1,
            'molecular_weight': 284.25,
            'logp': 1.74,
            'receptor_activity': {'5-HT2A': 8.5, '5-HT2C': 7.2, '5-HT1A': 6.8},
            'binding_affinity': 6.2,
            'notes': 'Psilocybin analog with improved oral bioavailability'
        },
        {
            'compound_name': 'OPI-2024-B2',
            'smiles': 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',
            'category': 'Opioid-2024-B2',
            'confidence_score': 0.92,
            'ip_status': 'Moderate Opportunity - Related patents exist',
            'source_research_topic_id': 2,
            'molecular_weight': 285.34,
            'logp': 1.23,
            'receptor_activity': {'μ-opioid': 9.1, 'δ-opioid': 6.5},
            'binding_affinity': 1.2,
            'notes': 'Morphine analog with reduced respiratory depression'
        },
        {
            'compound_name': 'ANX-2024-C3',
            'smiles': 'C1=CC=C(C=C1)C2=NC(=NC=N2)N3CCCC3',
            'category': 'Anxiolytic-2024-C3',
            'confidence_score': 0.88,
            'ip_status': 'High Opportunity - Novel mechanism',
            'source_research_topic_id': 3,
            'molecular_weight': 252.31,
            'logp': 2.45,
            'receptor_activity': {'GABA-A': 7.8, 'α1-adrenergic': 5.2},
            'binding_affinity': 15.6,
            'notes': 'Non-benzodiazepine GABA-A modulator'
        },
        {
            'compound_name': 'COG-2024-D4',
            'smiles': 'CC(C)(C)NC(=O)C1=CC=C(C=C1)S(=O)(=O)N',
            'category': 'Nootropic-2024-D4',
            'confidence_score': 0.82,
            'ip_status': 'Low Opportunity - Crowded field',
            'source_research_topic_id': 4,
            'molecular_weight': 284.36,
            'logp': 0.89,
            'receptor_activity': {'AMPA': 6.5, 'NMDA': 5.8},
            'binding_affinity': 125.4,
            'notes': 'Modafinil analog with enhanced cognitive effects'
        },
        {
            'compound_name': 'EMP-2024-E5',
            'smiles': 'CC(CC1=CC2=C(C=C1)OCO2)NC',
            'category': 'Empathogen-2024-E5',
            'confidence_score': 0.90,
            'ip_status': 'High Opportunity - Improved safety profile',
            'source_research_topic_id': 5,
            'molecular_weight': 193.24,
            'logp': 1.85,
            'receptor_activity': {'5-HT2A': 6.5, 'D2': 5.8, 'NET': 7.2},
            'binding_affinity': 158.0,
            'notes': 'MDMA analog with reduced neurotoxicity'
        },
        {
            'compound_name': 'PSI-2024-A2',
            'smiles': 'CCN(CC)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12',
            'category': 'Psychedelic-2024-A2',
            'confidence_score': 0.87,
            'ip_status': 'Moderate Opportunity',
            'source_research_topic_id': 1,
            'molecular_weight': 312.31,
            'logp': 2.15,
            'receptor_activity': {'5-HT2A': 7.8, '5-HT2C': 6.9},
            'binding_affinity': 12.4,
            'notes': 'N-ethyl psilocybin analog'
        },
        {
            'compound_name': 'ANX-2024-C4',
            'smiles': 'C1CN(CCN1)C2=NC3=CC=CC=C3N=C2Cl',
            'category': 'Anxiolytic-2024-C4',
            'confidence_score': 0.75,
            'ip_status': 'Low Opportunity',
            'source_research_topic_id': 3,
            'molecular_weight': 248.71,
            'logp': 1.95,
            'receptor_activity': {'GABA-A': 6.8},
            'binding_affinity': 28.7,
            'notes': 'Benzodiazepine derivative with reduced tolerance'
        }
    ]
    
    for discovery in discoveries:
        discovery_logger.log_discovery(**discovery)
    
    # 3. Add sample DOI references
    print("Adding DOI references...")
    
    doi_references = [
        {
            'doi': '10.1038/s41586-2024-07123-4',
            'title': 'Psilocybin for treatment-resistant depression: a randomised controlled trial',
            'authors': ['Goodwin, G.M.', 'Aaronson, S.T.', 'Alvarez, O.', 'Arden, P.C.'],
            'journal': 'Nature',
            'publication_year': 2024,
            'abstract': 'Psilocybin shows significant efficacy in treating treatment-resistant depression...',
            'keywords': ['psilocybin', 'depression', 'clinical trial', 'psychedelics'],
            'citation_count': 156,
            'impact_factor': 49.962,
            'research_area': 'Psychedelic Research',
            'compound_relevance': 'High',
            'confidence_score': 0.95
        },
        {
            'doi': '10.1016/j.neuron.2024.02.018',
            'title': 'Biased μ-opioid receptor signaling reduces addiction liability',
            'authors': ['Schmid, C.L.', 'Kennedy, N.M.', 'Ross, N.C.', 'Lovell, K.M.'],
            'journal': 'Neuron',
            'publication_year': 2024,
            'abstract': 'G-protein biased μ-opioid receptor agonists show reduced addiction potential...',
            'keywords': ['opioids', 'addiction', 'biased signaling', 'G-protein'],
            'citation_count': 89,
            'impact_factor': 17.173,
            'research_area': 'Opioid Research',
            'compound_relevance': 'High',
            'confidence_score': 0.92
        },
        {
            'doi': '10.1021/acs.jmedchem.2024.00456',
            'title': 'Novel GABA-A receptor modulators for anxiety disorders',
            'authors': ['Miller, P.S.', 'Aricescu, A.R.', 'Smart, T.G.'],
            'journal': 'Journal of Medicinal Chemistry',
            'publication_year': 2024,
            'abstract': 'Structure-based design of selective GABA-A receptor modulators...',
            'keywords': ['GABA-A', 'anxiety', 'allosteric modulation', 'selectivity'],
            'citation_count': 67,
            'impact_factor': 7.446,
            'research_area': 'Anxiolytic Research',
            'compound_relevance': 'High',
            'confidence_score': 0.88
        },
        {
            'doi': '10.1073/pnas.2024.121.e2315678121',
            'title': 'MDMA-assisted psychotherapy for PTSD: mechanisms and safety',
            'authors': ['Mitchell, J.M.', 'Bogenschutz, M.', 'Lilienstein, A.'],
            'journal': 'Proceedings of the National Academy of Sciences',
            'publication_year': 2024,
            'abstract': 'MDMA-assisted psychotherapy shows remarkable efficacy for PTSD treatment...',
            'keywords': ['MDMA', 'PTSD', 'psychotherapy', 'empathogen'],
            'citation_count': 134,
            'impact_factor': 12.779,
            'research_area': 'Empathogen Research',
            'compound_relevance': 'High',
            'confidence_score': 0.90
        },
        {
            'doi': '10.1016/j.tips.2024.03.007',
            'title': 'Cognitive enhancers for neurodegenerative diseases: current status and future directions',
            'authors': ['Ballard, C.', 'Aarsland, D.', 'Cummings, J.'],
            'journal': 'Trends in Pharmacological Sciences',
            'publication_year': 2024,
            'abstract': 'Review of current and emerging cognitive enhancement strategies...',
            'keywords': ['nootropics', 'Alzheimer', 'cognition', 'acetylcholine'],
            'citation_count': 78,
            'impact_factor': 13.065,
            'research_area': 'Nootropic Research',
            'compound_relevance': 'Medium',
            'confidence_score': 0.82
        }
    ]
    
    for doi_ref in doi_references:
        confidence_manager.add_doi_reference(**doi_ref)
    
    print("Sample data population complete!")
    print(f"Added {len(topics)} research topics")
    print(f"Added {len(discoveries)} compound discoveries")
    print(f"Added {len(doi_references)} DOI references")
    
    # Generate summary statistics
    high_confidence_dataset = confidence_manager.generate_confidence_dataset(0.9)
    all_dataset = confidence_manager.generate_confidence_dataset(0.0)
    
    print(f"\nSummary Statistics:")
    print(f"Total discoveries: {all_dataset.total_compounds}")
    print(f"High-confidence discoveries (>90%): {high_confidence_dataset.total_compounds}")
    print(f"Average confidence: {all_dataset.avg_confidence:.3f}")
    print(f"IP opportunities: {all_dataset.ip_opportunities}")
    print(f"Therapeutic areas: {list(all_dataset.therapeutic_areas.keys())}")

if __name__ == "__main__":
    populate_sample_data()
