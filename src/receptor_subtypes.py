#!/usr/bin/env python3
"""
Comprehensive Receptor Subtypes Database
Contains detailed subunit compositions and isoforms for major receptor families
"""

RECEPTOR_SUBTYPES = {
    # ========== GABA-A RECEPTOR SUBTYPES ==========
    "GABA-A_alpha1": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α1β2γ2 (most common)",
        "type": "Ion channel",
        "location": "Cerebral cortex, thalamus, hippocampus",
        "function": "Sedation, anticonvulsant, amnesia",
        "selective_ligands": {
            "agonists": ["Zolpidem", "Zaleplon", "Eszopiclone"],
            "antagonists": ["Beta-CCE", "3-PBC"]
        },
        "clinical_relevance": "Insomnia, epilepsy, sedation",
        "percentage_in_brain": "~60%"
    },
    "GABA-A_alpha2": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α2β3γ2",
        "type": "Ion channel",
        "location": "Limbic system, hippocampus, amygdala",
        "function": "Anxiolysis, muscle relaxation",
        "selective_ligands": {
            "agonists": ["L-838417", "TPA023", "MRK-409"],
            "antagonists": ["Flumazenil"]
        },
        "clinical_relevance": "Anxiety disorders, muscle spasticity",
        "percentage_in_brain": "~15-20%"
    },
    "GABA-A_alpha3": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α3β2/3γ2",
        "type": "Ion channel",
        "location": "Reticular thalamic nucleus, cortex",
        "function": "Sleep regulation, sensorimotor processing",
        "selective_ligands": {
            "agonists": ["TP003", "L-838417"],
            "modulators": ["NS11394"]
        },
        "clinical_relevance": "Absence epilepsy, autism spectrum disorders",
        "percentage_in_brain": "~10-15%"
    },
    "GABA-A_alpha4": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α4βδ (extrasynaptic)",
        "type": "Ion channel",
        "location": "Thalamus, dentate gyrus",
        "function": "Tonic inhibition, neurosteroid sensitivity",
        "selective_ligands": {
            "agonists": ["Gaboxadol", "DS2"],
            "modulators": ["Neurosteroids", "Allopregnanolone"]
        },
        "clinical_relevance": "Insomnia, premenstrual dysphoric disorder",
        "percentage_in_brain": "~5%"
    },
    "GABA-A_alpha5": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α5β3γ2",
        "type": "Ion channel",
        "location": "Hippocampus (25% of receptors)",
        "function": "Learning, memory, cognitive flexibility",
        "selective_ligands": {
            "inverse_agonists": ["α5IA", "RO4938581", "MRK-016"],
            "antagonists": ["L-655708", "PWZ-029"]
        },
        "clinical_relevance": "Cognitive enhancement, Down syndrome",
        "percentage_in_brain": "~5%"
    },
    "GABA-A_alpha6": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α6βδ or α6β2/3γ2",
        "type": "Ion channel",
        "location": "Cerebellar granule cells (exclusive)",
        "function": "Motor coordination, ethanol sensitivity",
        "selective_ligands": {
            "agonists": ["Ro15-4513 (α6-selective)"],
            "modulators": ["Furosemide"]
        },
        "clinical_relevance": "Ataxia, alcohol use disorder",
        "percentage_in_brain": "Cerebellum-specific"
    },
    "GABA-A_delta": {
        "family": "GABA",
        "parent": "GABA-A",
        "subunit_composition": "α4/6βδ",
        "type": "Ion channel",
        "location": "Extrasynaptic sites",
        "function": "Tonic inhibition, neurosteroid super-sensitivity",
        "selective_ligands": {
            "agonists": ["THIP (Gaboxadol)", "Muscimol"],
            "modulators": ["Allopregnanolone", "THDOC"]
        },
        "clinical_relevance": "PMS, postpartum depression, epilepsy",
        "note": "Insensitive to benzodiazepines"
    },
    "GABA-A_rho": {
        "family": "GABA",
        "parent": "GABA-C",
        "subunit_composition": "ρ1-3 homomers",
        "type": "Ion channel",
        "location": "Retina, superior colliculus",
        "function": "Visual processing, retinal signaling",
        "selective_ligands": {
            "agonists": ["CACA", "CAMP"],
            "antagonists": ["TPMPA", "Picrotoxin"]
        },
        "clinical_relevance": "Visual disorders, myopia",
        "note": "Previously called GABA-C receptors"
    },
    
    # ========== NMDA RECEPTOR SUBTYPES ==========
    "NMDA_GluN2A": {
        "family": "Glutamate",
        "parent": "NMDA",
        "subunit_composition": "GluN1/GluN2A",
        "type": "Ion channel",
        "location": "Cortex, hippocampus (adult)",
        "function": "Synaptic plasticity, LTP",
        "selective_ligands": {
            "antagonists": ["TCN-201", "MPX-004"],
            "modulators": ["GNE-8324"]
        },
        "clinical_relevance": "Cognitive enhancement, stroke",
        "kinetics": "Fast deactivation"
    },
    "NMDA_GluN2B": {
        "family": "Glutamate",
        "parent": "NMDA",
        "subunit_composition": "GluN1/GluN2B",
        "type": "Ion channel",
        "location": "Forebrain, extrasynaptic (early development)",
        "function": "Excitotoxicity, chronic pain",
        "selective_ligands": {
            "antagonists": ["Ifenprodil", "Ro 25-6981", "CP-101606", "EVT-101"],
            "modulators": ["Eliprodil"]
        },
        "clinical_relevance": "Depression, neuropathic pain, Alzheimer's",
        "kinetics": "Slow deactivation"
    },
    "NMDA_GluN2C": {
        "family": "Glutamate",
        "parent": "NMDA",
        "subunit_composition": "GluN1/GluN2C",
        "type": "Ion channel",
        "location": "Cerebellum, thalamus",
        "function": "Motor coordination",
        "selective_ligands": {
            "modulators": ["CIQ", "PYD-106"],
            "antagonists": ["QNZ46", "DQP-1105"]
        },
        "clinical_relevance": "Cerebellar ataxia, schizophrenia",
        "kinetics": "Intermediate"
    },
    "NMDA_GluN2D": {
        "family": "Glutamate",
        "parent": "NMDA",
        "subunit_composition": "GluN1/GluN2D",
        "type": "Ion channel",
        "location": "Brainstem, spinal cord, subthalamic nucleus",
        "function": "Motor control, sensory processing",
        "selective_ligands": {
            "modulators": ["CIQ", "PYD-106"],
            "antagonists": ["QNZ46", "DQP-1105"]
        },
        "clinical_relevance": "Parkinson's disease, pain",
        "kinetics": "Very slow deactivation"
    },
    "NMDA_GluN3": {
        "family": "Glutamate",
        "parent": "NMDA",
        "subunit_composition": "GluN1/GluN3A or GluN3A/B",
        "type": "Ion channel",
        "location": "Developing brain, adult olfactory bulb",
        "function": "Neuroprotection, synapse pruning",
        "selective_ligands": {
            "agonists": ["Glycine only (no glutamate required)"]
        },
        "clinical_relevance": "Neurodevelopmental disorders",
        "note": "Glycine-activated, Ca2+-impermeable"
    },
    
    # ========== AMPA RECEPTOR SUBTYPES ==========
    "AMPA_GluA1": {
        "family": "Glutamate",
        "parent": "AMPA",
        "subunit_composition": "GluA1 homomers or heteromers",
        "type": "Ion channel",
        "location": "Hippocampus CA1, cortex",
        "function": "LTP, synaptic plasticity",
        "calcium_permeability": "High (if lacking GluA2)",
        "selective_modulators": ["CX516", "CX614"],
        "clinical_relevance": "Memory disorders, depression"
    },
    "AMPA_GluA2": {
        "family": "Glutamate",
        "parent": "AMPA",
        "subunit_composition": "GluA2-containing heteromers",
        "type": "Ion channel",
        "location": "Widespread",
        "function": "Fast synaptic transmission",
        "calcium_permeability": "Low (Q/R editing)",
        "selective_modulators": ["LY404187", "PEPA"],
        "clinical_relevance": "ALS, epilepsy",
        "note": "Ca2+-impermeable due to RNA editing"
    },
    "AMPA_GluA3": {
        "family": "Glutamate",
        "parent": "AMPA",
        "subunit_composition": "GluA3-containing",
        "type": "Ion channel",
        "location": "Cortex, hippocampus",
        "function": "Synaptic transmission",
        "clinical_relevance": "Autism, intellectual disability"
    },
    "AMPA_GluA4": {
        "family": "Glutamate",
        "parent": "AMPA",
        "subunit_composition": "GluA4-containing",
        "type": "Ion channel",
        "location": "Cerebellum, early development",
        "function": "Fast kinetics, development",
        "clinical_relevance": "Cerebellar disorders"
    },
    
    # ========== NICOTINIC RECEPTOR SUBTYPES ==========
    "nAChR_alpha4beta2": {
        "family": "Nicotinic",
        "parent": "nAChR",
        "subunit_composition": "(α4)2(β2)3 or (α4)3(β2)2",
        "type": "Ion channel",
        "location": "Throughout brain",
        "function": "Addiction, attention, analgesia",
        "selective_ligands": {
            "agonists": ["Varenicline", "Cytisine", "ABT-594"],
            "antagonists": ["DHβE"]
        },
        "clinical_relevance": "Smoking cessation, ADHD, pain"
    },
    "nAChR_alpha7": {
        "family": "Nicotinic",
        "parent": "nAChR",
        "subunit_composition": "(α7)5 homomer",
        "type": "Ion channel",
        "location": "Hippocampus, cortex",
        "function": "Cognition, anti-inflammation",
        "calcium_permeability": "Very high",
        "selective_ligands": {
            "agonists": ["PNU-282987", "GTS-21", "AR-R17779"],
            "antagonists": ["MLA", "α-bungarotoxin"],
            "PAMs": ["PNU-120596", "AVL-3288"]
        },
        "clinical_relevance": "Alzheimer's, schizophrenia"
    },
    "nAChR_alpha3beta4": {
        "family": "Nicotinic",
        "parent": "nAChR",
        "subunit_composition": "(α3)2(β4)3",
        "type": "Ion channel",
        "location": "Autonomic ganglia, medial habenula",
        "function": "Autonomic function, nicotine aversion",
        "selective_ligands": {
            "antagonists": ["AT-1001", "Mecamylamine"]
        },
        "clinical_relevance": "Nicotine dependence, autonomic disorders"
    },
    "nAChR_alpha6beta2": {
        "family": "Nicotinic",
        "parent": "nAChR",
        "subunit_composition": "α6α4β2β3 or α6β2",
        "type": "Ion channel",
        "location": "Dopaminergic neurons, retina",
        "function": "Dopamine release, vision",
        "selective_ligands": {
            "antagonists": ["α-conotoxin MII"]
        },
        "clinical_relevance": "Parkinson's disease, nicotine addiction"
    },
    "nAChR_muscle": {
        "family": "Nicotinic",
        "parent": "nAChR",
        "subunit_composition": "(α1)2β1δε (adult) or (α1)2β1δγ (fetal)",
        "type": "Ion channel",
        "location": "Neuromuscular junction",
        "function": "Muscle contraction",
        "selective_ligands": {
            "antagonists": ["Tubocurarine", "Pancuronium", "Vecuronium"],
            "agonists": ["Succinylcholine"]
        },
        "clinical_relevance": "Myasthenia gravis, anesthesia"
    },
    
    # ========== 5-HT3 RECEPTOR SUBTYPES ==========
    "5-HT3A": {
        "family": "Serotonin",
        "parent": "5-HT3",
        "subunit_composition": "(5-HT3A)5 homomer",
        "type": "Ion channel",
        "location": "Area postrema, vagus nerve",
        "function": "Nausea, vomiting",
        "selective_ligands": {
            "antagonists": ["Ondansetron", "Granisetron"]
        },
        "clinical_relevance": "Antiemetic therapy"
    },
    "5-HT3AB": {
        "family": "Serotonin",
        "parent": "5-HT3",
        "subunit_composition": "5-HT3A/5-HT3B heteromer",
        "type": "Ion channel",
        "location": "Hippocampus, amygdala",
        "function": "Mood, cognition",
        "conductance": "Higher than 5-HT3A alone",
        "clinical_relevance": "Anxiety, IBS"
    },
    
    # ========== P2X RECEPTOR SUBTYPES ==========
    "P2X1": {
        "family": "Purinergic",
        "parent": "P2X",
        "subunit_composition": "(P2X1)3 homomer",
        "type": "Ion channel",
        "location": "Smooth muscle, platelets",
        "function": "Vasoconstriction, platelet aggregation",
        "desensitization": "Very fast",
        "selective_antagonists": ["NF449", "RO0437626"],
        "clinical_relevance": "Thrombosis, bladder disorders"
    },
    "P2X2": {
        "family": "Purinergic",
        "parent": "P2X",
        "subunit_composition": "(P2X2)3 or P2X2/3",
        "type": "Ion channel",
        "location": "CNS, sensory neurons",
        "function": "Neurotransmission, taste",
        "desensitization": "Slow",
        "clinical_relevance": "Hearing, taste disorders"
    },
    "P2X3": {
        "family": "Purinergic",
        "parent": "P2X",
        "subunit_composition": "(P2X3)3 or P2X2/3",
        "type": "Ion channel",
        "location": "Sensory neurons, nociceptors",
        "function": "Pain sensation, bladder reflexes",
        "selective_antagonists": ["A-317491", "AF-353", "Gefapixant"],
        "clinical_relevance": "Chronic cough, pain, bladder disorders"
    },
    "P2X4": {
        "family": "Purinergic",
        "parent": "P2X",
        "subunit_composition": "(P2X4)3",
        "type": "Ion channel",
        "location": "Microglia, endothelium",
        "function": "Neuropathic pain, inflammation",
        "modulation": "Potentiated by ivermectin",
        "selective_antagonists": ["5-BDBD", "PSB-12062"],
        "clinical_relevance": "Neuropathic pain, inflammation"
    },
    "P2X7": {
        "family": "Purinergic",
        "parent": "P2X",
        "subunit_composition": "(P2X7)3",
        "type": "Ion channel",
        "location": "Immune cells, microglia",
        "function": "Inflammation, cell death",
        "unique_feature": "Forms large pore with prolonged activation",
        "selective_antagonists": ["A-740003", "JNJ-47965567", "AZD9056"],
        "clinical_relevance": "Inflammation, depression, pain"
    },
    
    # ========== KAINATE RECEPTOR SUBTYPES ==========
    "KA_GluK1": {
        "family": "Glutamate",
        "parent": "Kainate",
        "subunit_composition": "GluK1 (formerly GluR5)",
        "type": "Ion channel",
        "location": "Hippocampus, cortex",
        "function": "Synaptic modulation",
        "selective_antagonists": ["UBP302", "ACET"],
        "clinical_relevance": "Epilepsy, pain"
    },
    "KA_GluK2": {
        "family": "Glutamate",
        "parent": "Kainate",
        "subunit_composition": "GluK2 (formerly GluR6)",
        "type": "Ion channel",
        "location": "Cerebellum, hippocampus",
        "function": "Synaptic plasticity",
        "clinical_relevance": "Autism, intellectual disability"
    },
    "KA_GluK3": {
        "family": "Glutamate",
        "parent": "Kainate",
        "subunit_composition": "GluK3 (formerly GluR7)",
        "type": "Ion channel",
        "location": "Cortex, retina",
        "function": "Neurotransmission modulation",
        "clinical_relevance": "Bipolar disorder, schizophrenia"
    }
}

def get_subtype_info(subtype_name: str) -> dict:
    """Get detailed information about a specific receptor subtype"""
    return RECEPTOR_SUBTYPES.get(subtype_name, {})

def get_subtypes_by_parent(parent_receptor: str) -> list:
    """Get all subtypes of a parent receptor"""
    return [name for name, data in RECEPTOR_SUBTYPES.items() 
            if data.get('parent', '') == parent_receptor]

def get_subtypes_by_location(brain_region: str) -> list:
    """Get receptor subtypes expressed in a specific brain region"""
    region = brain_region.lower()
    return [name for name, data in RECEPTOR_SUBTYPES.items()
            if region in data.get('location', '').lower()]

def get_benzodiazepine_sensitive_subtypes() -> list:
    """Get GABA-A subtypes sensitive to benzodiazepines"""
    benzo_sensitive = []
    for name, data in RECEPTOR_SUBTYPES.items():
        if data.get('family') == 'GABA' and 'γ2' in data.get('subunit_composition', ''):
            if 'delta' not in name.lower() and 'rho' not in name.lower():
                benzo_sensitive.append(name)
    return benzo_sensitive

def get_calcium_permeable_subtypes() -> list:
    """Get calcium-permeable receptor subtypes"""
    ca_permeable = []
    for name, data in RECEPTOR_SUBTYPES.items():
        if data.get('calcium_permeability') in ['High', 'Very high']:
            ca_permeable.append(name)
        elif 'NMDA' in name:  # All NMDA receptors are Ca-permeable
            ca_permeable.append(name)
    return ca_permeable

def get_extrasynaptic_subtypes() -> list:
    """Get extrasynaptic receptor subtypes"""
    extrasynaptic = []
    for name, data in RECEPTOR_SUBTYPES.items():
        if 'extrasynaptic' in data.get('location', '').lower() or \
           'extrasynaptic' in data.get('subunit_composition', '').lower():
            extrasynaptic.append(name)
    return extrasynaptic