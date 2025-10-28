#!/usr/bin/env python3
"""
Comprehensive Receptor Database
Contains all major receptor types, subtypes, and their binding profiles
"""

RECEPTOR_DATABASE = {
    # ========== SEROTONIN RECEPTORS ==========
    "5-HT1A": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["7E2X", "7E2Y", "5V54"],
        "function": "Anxiolytic, antidepressant",
        "location": "Raphe nuclei, hippocampus, cortex",
        "agonists": ["Buspirone", "Aripiprazole", "8-OH-DPAT"],
        "antagonists": ["WAY-100635", "Pindolol"],
        "therapeutic_relevance": "Anxiety, depression, schizophrenia"
    },
    "5-HT1B": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4IAQ", "4IAR", "5V54"],
        "function": "Vasoconstriction, anti-migraine",
        "location": "Blood vessels, substantia nigra",
        "agonists": ["Sumatriptan", "Eletriptan", "Zolmitriptan"],
        "antagonists": ["SB-216641", "GR-127935"],
        "therapeutic_relevance": "Migraine, cluster headaches"
    },
    "5-HT1D": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "function": "Neurotransmitter release inhibition",
        "location": "Basal ganglia, hippocampus",
        "agonists": ["Sumatriptan", "PNU-109291"],
        "therapeutic_relevance": "Migraine"
    },
    "5-HT2A": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["6A93", "6A94", "6WGT", "7VOD"],
        "function": "Hallucinations, perception, mood",
        "location": "Cortex, claustrum",
        "agonists": ["LSD", "Psilocin", "DMT", "Mescaline", "DOI"],
        "antagonists": ["Ketanserin", "Risperidone", "MDL-100907"],
        "inverse_agonists": ["Clozapine", "Olanzapine"],
        "therapeutic_relevance": "Psychedelics, schizophrenia, depression"
    },
    "5-HT2B": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4IB4", "5TUD", "6DRX"],
        "function": "Cardiovascular regulation",
        "location": "Heart, stomach",
        "agonists": ["BW723C86", "Ro 60-0175"],
        "antagonists": ["SB-204741", "RS-127445"],
        "therapeutic_relevance": "Heart valve disease (avoid activation)"
    },
    "5-HT2C": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["6BQG", "6BQH", "7VOE"],
        "function": "Appetite, anxiety, addiction",
        "location": "Choroid plexus, limbic system",
        "agonists": ["Lorcaserin", "WAY-163909", "Ro 60-0175"],
        "antagonists": ["Agomelatine", "SB-242084"],
        "therapeutic_relevance": "Obesity, anxiety, addiction"
    },
    "5-HT3": {
        "family": "Serotonin",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["4PIR", "6Y1Z", "6Y20"],
        "function": "Nausea, vomiting, gut motility",
        "location": "Area postrema, GI tract",
        "agonists": ["SR-57227", "m-CPBG"],
        "antagonists": ["Ondansetron", "Granisetron", "Palonosetron"],
        "therapeutic_relevance": "Antiemetic, IBS"
    },
    "5-HT4": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "function": "GI motility, memory, mood",
        "location": "GI tract, hippocampus",
        "agonists": ["Prucalopride", "Tegaserod", "BIMU-8"],
        "antagonists": ["GR-113808", "SB-204070"],
        "therapeutic_relevance": "Constipation, cognitive enhancement"
    },
    "5-HT5A": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "function": "Memory, circadian rhythms",
        "location": "Hippocampus, cortex"
    },
    "5-HT6": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["7XTB", "7XTC"],
        "function": "Cognition, mood",
        "location": "Striatum, cortex",
        "antagonists": ["Idalopirdine", "SB-271046", "Ro 04-6790"],
        "therapeutic_relevance": "Alzheimer's, cognitive enhancement"
    },
    "5-HT7": {
        "family": "Serotonin",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["7XIJ", "7XIK"],
        "function": "Sleep, thermoregulation, mood",
        "location": "Thalamus, hypothalamus",
        "agonists": ["AS-19", "LP-211"],
        "antagonists": ["SB-269970", "DR-4004"],
        "therapeutic_relevance": "Depression, sleep disorders"
    },
    
    # ========== DOPAMINE RECEPTORS ==========
    "D1": {
        "family": "Dopamine",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["7CKW", "7CKX", "7CKY"],
        "function": "Motor control, reward, cognition",
        "location": "Striatum, nucleus accumbens, cortex",
        "agonists": ["SKF-38393", "SKF-81297", "Fenoldopam"],
        "antagonists": ["SCH-23390", "SKF-83566"],
        "therapeutic_relevance": "Parkinson's, ADHD, cognitive disorders"
    },
    "D2": {
        "family": "Dopamine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["6CM4", "6LUQ", "7DFP", "S0S"],
        "function": "Motor control, prolactin release, psychosis",
        "location": "Striatum, pituitary, substantia nigra",
        "agonists": ["Bromocriptine", "Ropinirole", "Pramipexole", "Cabergoline"],
        "antagonists": ["Haloperidol", "Risperidone", "Sulpiride"],
        "partial_agonists": ["Aripiprazole", "Brexpiprazole"],
        "therapeutic_relevance": "Schizophrenia, Parkinson's, hyperprolactinemia"
    },
    "D3": {
        "family": "Dopamine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["3PBL", "7CMV"],
        "function": "Reward, emotion, addiction",
        "location": "Nucleus accumbens, islands of Calleja",
        "agonists": ["Pramipexole", "7-OH-DPAT", "PD-128907"],
        "antagonists": ["SB-277011A", "S33084"],
        "therapeutic_relevance": "Addiction, schizophrenia, Parkinson's"
    },
    "D4": {
        "family": "Dopamine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["5WIU", "5WIV"],
        "function": "Cognition, attention",
        "location": "Frontal cortex, hippocampus",
        "agonists": ["PD-168077", "CP-226269"],
        "antagonists": ["L-745870", "Clozapine"],
        "therapeutic_relevance": "ADHD, schizophrenia"
    },
    "D5": {
        "family": "Dopamine",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "function": "Similar to D1, higher affinity",
        "location": "Hippocampus, hypothalamus",
        "therapeutic_relevance": "Hypertension, neuropsychiatric disorders"
    },
    
    # ========== ADRENERGIC RECEPTORS ==========
    "Alpha1A": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "function": "Vasoconstriction, prostate smooth muscle",
        "agonists": ["Phenylephrine", "Methoxamine"],
        "antagonists": ["Tamsulosin", "Silodosin", "Prazosin"],
        "therapeutic_relevance": "BPH, hypertension"
    },
    "Alpha1B": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "function": "Vascular smooth muscle contraction",
        "therapeutic_relevance": "Hypertension"
    },
    "Alpha2A": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["6KUX", "6KUW"],
        "function": "Presynaptic inhibition, sedation, analgesia",
        "location": "Locus coeruleus, brainstem",
        "agonists": ["Clonidine", "Dexmedetomidine", "Guanfacine"],
        "antagonists": ["Yohimbine", "Atipamezole"],
        "therapeutic_relevance": "Hypertension, ADHD, sedation, pain"
    },
    "Beta1": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["7BU6", "7BU7", "7BVQ"],
        "function": "Cardiac stimulation, renin release",
        "location": "Heart, kidney",
        "agonists": ["Dobutamine", "Xamoterol"],
        "antagonists": ["Atenolol", "Metoprolol", "Bisoprolol"],
        "therapeutic_relevance": "Hypertension, heart failure, arrhythmias"
    },
    "Beta2": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["2RH1", "3NY8", "3NY9", "4LDE"],
        "function": "Bronchodilation, glycogenolysis",
        "location": "Lungs, liver, skeletal muscle",
        "agonists": ["Salbutamol", "Formoterol", "Salmeterol"],
        "antagonists": ["ICI-118551", "Butoxamine"],
        "therapeutic_relevance": "Asthma, COPD"
    },
    "Beta3": {
        "family": "Adrenergic",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["7BTS"],
        "function": "Lipolysis, bladder relaxation",
        "location": "Adipose tissue, bladder",
        "agonists": ["Mirabegron", "Solabegron", "Vibegron"],
        "therapeutic_relevance": "Overactive bladder, obesity"
    },
    
    # ========== OPIOID RECEPTORS ==========
    "MOR": {
        "family": "Opioid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4DKL", "5C1M", "6DDE", "7T2G"],
        "function": "Analgesia, euphoria, respiratory depression",
        "location": "Periaqueductal gray, thalamus, spinal cord",
        "agonists": ["Morphine", "Fentanyl", "Oxycodone", "DAMGO"],
        "antagonists": ["Naloxone", "Naltrexone", "Methylnaltrexone"],
        "partial_agonists": ["Buprenorphine"],
        "biased_agonists": ["PZM21", "SR-17018"],
        "therapeutic_relevance": "Pain, addiction, constipation"
    },
    "DOR": {
        "family": "Opioid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4N6H", "4RWD", "6PT2"],
        "function": "Analgesia, mood, neuroprotection",
        "location": "Cortex, striatum, spinal cord",
        "agonists": ["DPDPE", "SNC80", "Deltorphin"],
        "antagonists": ["Naltrindole", "TIPP"],
        "therapeutic_relevance": "Pain, depression, neuroprotection"
    },
    "KOR": {
        "family": "Opioid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4DJH", "6VI4", "7Y0P"],
        "function": "Analgesia, dysphoria, diuresis",
        "location": "Hypothalamus, claustrum, spinal cord",
        "agonists": ["U-69593", "Salvinorin A", "Dynorphin"],
        "antagonists": ["Norbinaltorphimine", "JDTic", "CERC-501"],
        "therapeutic_relevance": "Pain, addiction, depression"
    },
    "NOP": {
        "family": "Opioid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4EA3", "5DHG"],
        "function": "Analgesia, anxiety modulation",
        "location": "Cortex, amygdala, spinal cord",
        "agonists": ["Nociceptin", "Ro 64-6198", "SCH-221510"],
        "antagonists": ["J-113397", "SB-612111"],
        "therapeutic_relevance": "Pain, anxiety, addiction"
    },
    
    # ========== GABA RECEPTORS ==========
    "GABA-A": {
        "family": "GABA",
        "type": "Ion channel",
        "subtype": "Ligand-gated chloride channel",
        "pdb_ids": ["6HUO", "6HUP", "6X3T"],
        "function": "Fast inhibitory neurotransmission",
        "location": "Throughout CNS",
        "agonists": ["GABA", "Muscimol", "Isoguvacine"],
        "positive_modulators": ["Diazepam", "Alprazolam", "Zolpidem", "Propofol"],
        "antagonists": ["Bicuculline", "Gabazine", "Flumazenil"],
        "therapeutic_relevance": "Anxiety, insomnia, epilepsy, anesthesia"
    },
    "GABA-B": {
        "family": "GABA",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4MQE", "4MR7", "6UO8"],
        "function": "Slow inhibitory neurotransmission",
        "location": "Pre and postsynaptic terminals",
        "agonists": ["Baclofen", "GABA", "SKF-97541"],
        "antagonists": ["Phaclofen", "Saclofen", "CGP-55845"],
        "therapeutic_relevance": "Spasticity, addiction, pain"
    },
    
    # ========== GLUTAMATE RECEPTORS ==========
    "NMDA": {
        "family": "Glutamate",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["4PE5", "5UN0", "6WHT"],
        "function": "Synaptic plasticity, memory, excitotoxicity",
        "location": "Hippocampus, cortex",
        "agonists": ["NMDA", "Glutamate", "Glycine", "D-serine"],
        "antagonists": ["Ketamine", "PCP", "MK-801", "Memantine"],
        "therapeutic_relevance": "Depression, Alzheimer's, neuroprotection"
    },
    "AMPA": {
        "family": "Glutamate",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["5WEO", "6DM1", "6QKC"],
        "function": "Fast excitatory neurotransmission",
        "location": "Throughout CNS",
        "agonists": ["AMPA", "Glutamate", "Quisqualate"],
        "positive_modulators": ["Aniracetam", "CX614", "Ampakines"],
        "antagonists": ["CNQX", "NBQX", "Perampanel"],
        "therapeutic_relevance": "Cognitive enhancement, epilepsy"
    },
    "Kainate": {
        "family": "Glutamate",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["3G3F", "4E0X", "5KUF"],
        "function": "Modulation of synaptic transmission",
        "agonists": ["Kainate", "Domoate", "ATPA"],
        "antagonists": ["CNQX", "UBP302"],
        "therapeutic_relevance": "Epilepsy, pain"
    },
    "mGluR1": {
        "family": "Glutamate",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4OR2"],
        "function": "Synaptic plasticity, motor coordination",
        "location": "Cerebellum, hippocampus",
        "agonists": ["Quisqualate", "DHPG"],
        "antagonists": ["JNJ16259685", "BAY 36-7620"],
        "therapeutic_relevance": "Ataxia, neurodegeneration"
    },
    "mGluR5": {
        "family": "Glutamate",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4OO9", "5CGC", "6FFI"],
        "function": "Synaptic plasticity, anxiety, addiction",
        "location": "Hippocampus, cortex, striatum",
        "positive_modulators": ["CDPPB", "VU0409551"],
        "negative_modulators": ["MPEP", "Fenobam", "Mavoglurant"],
        "therapeutic_relevance": "Fragile X, autism, addiction, anxiety"
    },
    
    # ========== CANNABINOID RECEPTORS ==========
    "CB1": {
        "family": "Cannabinoid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["5TGZ", "5XRA", "6KPC"],
        "function": "Appetite, pain, memory, mood",
        "location": "Brain (highest density in basal ganglia)",
        "agonists": ["THC", "Anandamide", "2-AG", "WIN 55,212-2"],
        "antagonists": ["Rimonabant", "Taranabant", "AM251"],
        "therapeutic_relevance": "Pain, appetite, epilepsy, PTSD"
    },
    "CB2": {
        "family": "Cannabinoid",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["5ZTY", "6KPC", "6KPF"],
        "function": "Inflammation, immune modulation, pain",
        "location": "Immune cells, microglia",
        "agonists": ["JWH-133", "HU-308", "AM1241"],
        "antagonists": ["AM630", "SR144528"],
        "therapeutic_relevance": "Inflammation, pain, neuroinflammation"
    },
    
    # ========== CHOLINERGIC RECEPTORS ==========
    "M1": {
        "family": "Muscarinic",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["6OIJ", "6WJG"],
        "function": "Cognition, memory, learning",
        "location": "Cortex, hippocampus",
        "agonists": ["Xanomeline", "Cevimeline", "McN-A-343"],
        "antagonists": ["Pirenzepine", "Telenzepine"],
        "therapeutic_relevance": "Alzheimer's, schizophrenia"
    },
    "M2": {
        "family": "Muscarinic",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["3UON", "4MQS"],
        "function": "Cardiac inhibition, neurotransmitter release",
        "location": "Heart, smooth muscle, CNS",
        "agonists": ["Pilocarpine", "Carbachol"],
        "antagonists": ["Methoctramine", "AF-DX 116"],
        "therapeutic_relevance": "Bradycardia, COPD"
    },
    "M3": {
        "family": "Muscarinic",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4DAJ", "4U14"],
        "function": "Smooth muscle contraction, glandular secretion",
        "location": "Smooth muscle, glands",
        "antagonists": ["Darifenacin", "Solifenacin", "Tiotropium"],
        "therapeutic_relevance": "Overactive bladder, COPD"
    },
    "nAChR_alpha7": {
        "family": "Nicotinic",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["7KOX", "7KOQ"],
        "function": "Cognition, inflammation modulation",
        "location": "Hippocampus, cortex",
        "agonists": ["PNU-282987", "EVP-6124", "GTS-21"],
        "antagonists": ["Alpha-bungarotoxin", "MLA"],
        "positive_modulators": ["PNU-120596", "AVL-3288"],
        "therapeutic_relevance": "Alzheimer's, schizophrenia, inflammation"
    },
    "nAChR_alpha4beta2": {
        "family": "Nicotinic",
        "type": "Ion channel",
        "subtype": "Ligand-gated cation channel",
        "pdb_ids": ["5KXI", "6CNK"],
        "function": "Attention, addiction, analgesia",
        "location": "Throughout brain",
        "agonists": ["Nicotine", "Varenicline", "Cytisine"],
        "antagonists": ["Mecamylamine", "DHβE"],
        "therapeutic_relevance": "Smoking cessation, ADHD, pain"
    },
    
    # ========== HISTAMINE RECEPTORS ==========
    "H1": {
        "family": "Histamine",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["3RZE", "7DFL"],
        "function": "Allergy, inflammation, wakefulness",
        "location": "Smooth muscle, endothelium, CNS",
        "agonists": ["Histamine", "2-Methylhistamine"],
        "antagonists": ["Diphenhydramine", "Cetirizine", "Loratadine"],
        "therapeutic_relevance": "Allergies, insomnia, motion sickness"
    },
    "H2": {
        "family": "Histamine",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "function": "Gastric acid secretion",
        "location": "Gastric parietal cells, heart",
        "agonists": ["Histamine", "Dimaprit"],
        "antagonists": ["Ranitidine", "Famotidine", "Cimetidine"],
        "therapeutic_relevance": "Peptic ulcers, GERD"
    },
    "H3": {
        "family": "Histamine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "function": "Neurotransmitter release modulation",
        "location": "CNS presynaptic terminals",
        "agonists": ["R-α-methylhistamine", "Immepip"],
        "antagonists": ["Pitolisant", "ABT-239"],
        "therapeutic_relevance": "Narcolepsy, cognitive disorders"
    },
    "H4": {
        "family": "Histamine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["7YFN"],
        "function": "Immune modulation, inflammation",
        "location": "Immune cells, bone marrow",
        "agonists": ["4-Methylhistamine", "VUF-8430"],
        "antagonists": ["JNJ-7777120", "Toreforant"],
        "therapeutic_relevance": "Inflammatory diseases, pruritus"
    },
    
    # ========== OTHER IMPORTANT RECEPTORS ==========
    "Adenosine_A1": {
        "family": "Adenosine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["6D9H"],
        "function": "Cardiac protection, sedation",
        "agonists": ["Adenosine", "CPA", "GR79236"],
        "antagonists": ["Caffeine", "Theophylline", "DPCPX"],
        "therapeutic_relevance": "Arrhythmias, sleep, neuroprotection"
    },
    "Adenosine_A2A": {
        "family": "Adenosine",
        "type": "GPCR",
        "subtype": "Gs-coupled",
        "pdb_ids": ["2YDO", "3EML", "5G53"],
        "function": "Vasodilation, immune modulation",
        "agonists": ["Regadenoson", "CGS-21680"],
        "antagonists": ["Istradefylline", "Caffeine", "SCH-58261"],
        "therapeutic_relevance": "Parkinson's, cardiac imaging"
    },
    "TRPV1": {
        "family": "TRP",
        "type": "Ion channel",
        "subtype": "Ligand/heat-gated cation channel",
        "pdb_ids": ["3J5P", "5IRZ", "7L2H"],
        "function": "Pain, heat sensation, inflammation",
        "location": "Sensory neurons",
        "agonists": ["Capsaicin", "Resiniferatoxin", "Anandamide"],
        "antagonists": ["Capsazepine", "AMG-517", "ABT-102"],
        "therapeutic_relevance": "Pain, inflammation, thermoregulation"
    },
    "P2X7": {
        "family": "Purinergic",
        "type": "Ion channel",
        "subtype": "ATP-gated cation channel",
        "pdb_ids": ["5U1U", "5U1V", "6U9V"],
        "function": "Inflammation, cell death, pain",
        "location": "Immune cells, microglia",
        "agonists": ["ATP", "BzATP"],
        "antagonists": ["A-740003", "JNJ-47965567", "AZ-10606120"],
        "therapeutic_relevance": "Inflammation, depression, pain"
    },
    "Orexin_OX1": {
        "family": "Orexin",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4ZJ8", "6TO7"],
        "function": "Wakefulness, appetite, reward",
        "agonists": ["Orexin-A", "Orexin-B"],
        "antagonists": ["SB-334867", "Suvorexant"],
        "therapeutic_relevance": "Insomnia, narcolepsy, addiction"
    },
    "Orexin_OX2": {
        "family": "Orexin",
        "type": "GPCR",
        "subtype": "Gq/11-coupled",
        "pdb_ids": ["4S0V", "5WQC", "7L1U"],
        "function": "Sleep/wake regulation",
        "antagonists": ["Suvorexant", "Lemborexant", "Daridorexant"],
        "therapeutic_relevance": "Insomnia"
    },
    "Melatonin_MT1": {
        "family": "Melatonin",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["6ME2", "6ME3"],
        "function": "Circadian rhythm, sleep",
        "agonists": ["Melatonin", "Ramelteon", "Tasimelteon"],
        "antagonists": ["Luzindole"],
        "therapeutic_relevance": "Insomnia, jet lag"
    },
    "CCR5": {
        "family": "Chemokine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["4MBS", "5UIW", "6AKX"],
        "function": "HIV entry, inflammation",
        "antagonists": ["Maraviroc", "Cenicriviroc", "Vicriviroc"],
        "therapeutic_relevance": "HIV, inflammation"
    },
    "CXCR4": {
        "family": "Chemokine",
        "type": "GPCR",
        "subtype": "Gi/Go-coupled",
        "pdb_ids": ["3ODU", "3OE0", "4RWS"],
        "function": "HIV entry, stem cell homing, cancer",
        "antagonists": ["Plerixafor", "AMD3100", "IT1t"],
        "therapeutic_relevance": "HIV, cancer, stem cell mobilization"
    }
}

# Receptor classification system
RECEPTOR_CLASSIFICATIONS = {
    "by_family": {
        "Monoamine": ["5-HT", "Dopamine", "Adrenergic", "Histamine"],
        "Amino_Acid": ["GABA", "Glutamate", "Glycine"],
        "Peptide": ["Opioid", "Orexin", "Chemokine", "Neuropeptide"],
        "Lipid": ["Cannabinoid", "Prostanoid", "Leukotriene"],
        "Purine": ["Adenosine", "P2X", "P2Y"],
        "Other": ["Melatonin", "Sigma", "TAAR"]
    },
    "by_type": {
        "GPCR": {
            "Gs-coupled": ["Beta1", "Beta2", "Beta3", "D1", "D5", "5-HT4", "5-HT6", "5-HT7", "H2"],
            "Gi/Go-coupled": ["D2", "D3", "D4", "5-HT1A", "5-HT1B", "Alpha2A", "MOR", "DOR", "KOR"],
            "Gq/11-coupled": ["5-HT2A", "5-HT2B", "5-HT2C", "Alpha1A", "M1", "M3", "H1"],
            "G12/13-coupled": ["S1P", "LPA", "Thromboxane"]
        },
        "Ion_Channel": {
            "Ligand-gated": ["GABA-A", "NMDA", "AMPA", "nAChR", "5-HT3", "P2X"],
            "Voltage-gated": ["Nav", "Cav", "Kv", "HCN"],
            "TRP": ["TRPV1", "TRPA1", "TRPM8"]
        },
        "Nuclear": ["ER", "AR", "GR", "MR", "PR", "PPAR", "RXR", "RAR"],
        "Enzyme-linked": ["Insulin", "EGFR", "VEGFR", "BDNF"]
    },
    "by_therapeutic_area": {
        "Psychiatry": ["5-HT2A", "D2", "5-HT1A", "GABA-A", "NMDA"],
        "Pain": ["MOR", "TRPV1", "CB1", "Alpha2A", "Nav1.7"],
        "Cardiovascular": ["Beta1", "Alpha1A", "AT1", "Adenosine_A1"],
        "Respiratory": ["Beta2", "M3", "H1", "CysLT1"],
        "GI": ["5-HT3", "5-HT4", "M3", "H2", "CB1"],
        "Inflammation": ["CB2", "H4", "P2X7", "CCR5", "TNF"],
        "Neurodegenerative": ["M1", "NMDA", "Alpha7", "5-HT6", "Adenosine_A2A"],
        "Addiction": ["MOR", "D2", "D3", "CB1", "nAChR_alpha4beta2"],
        "Sleep": ["GABA-A", "Orexin_OX2", "Melatonin_MT1", "H1", "5-HT2A"]
    }
}

def get_receptor_info(receptor_name: str) -> dict:
    """Get detailed information about a specific receptor"""
    return RECEPTOR_DATABASE.get(receptor_name, {})

def get_receptors_by_family(family: str) -> list:
    """Get all receptors belonging to a specific family"""
    return [name for name, data in RECEPTOR_DATABASE.items() 
            if data.get('family', '').lower() == family.lower()]

def get_receptors_by_therapeutic_area(area: str) -> list:
    """Get receptors relevant to a therapeutic area"""
    return RECEPTOR_CLASSIFICATIONS['by_therapeutic_area'].get(area, [])

def search_receptors(query: str) -> list:
    """Search receptors by name, family, or function"""
    query = query.lower()
    results = []
    
    for receptor_name, data in RECEPTOR_DATABASE.items():
        if (query in receptor_name.lower() or 
            query in data.get('family', '').lower() or
            query in data.get('function', '').lower() or
            query in data.get('therapeutic_relevance', '').lower()):
            results.append(receptor_name)
    
    return results

def get_receptor_ligands(receptor_name: str, ligand_type: str = 'all') -> list:
    """Get ligands for a specific receptor"""
    receptor = RECEPTOR_DATABASE.get(receptor_name, {})
    
    if ligand_type == 'agonists':
        return receptor.get('agonists', [])
    elif ligand_type == 'antagonists':
        return receptor.get('antagonists', [])
    elif ligand_type == 'all':
        ligands = []
        ligands.extend(receptor.get('agonists', []))
        ligands.extend(receptor.get('antagonists', []))
        ligands.extend(receptor.get('partial_agonists', []))
        ligands.extend(receptor.get('inverse_agonists', []))
        ligands.extend(receptor.get('positive_modulators', []))
        ligands.extend(receptor.get('negative_modulators', []))
        return list(set(ligands))
    
    return []

def get_all_receptor_targets() -> dict:
    """Get all receptor targets including subtypes"""
    from src.receptor_subtypes import RECEPTOR_SUBTYPES
    
    # Combine main receptors and subtypes
    all_targets = {}
    all_targets.update(RECEPTOR_DATABASE)
    all_targets.update(RECEPTOR_SUBTYPES)
    
    return all_targets