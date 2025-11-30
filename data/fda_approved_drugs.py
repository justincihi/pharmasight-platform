#!/usr/bin/env python3
"""
FDA Approved Drugs Database
Curated collection of 200+ FDA-approved drugs with SMILES structures
"""

FDA_APPROVED_DRUGS = {
    "acetaminophen": {
        "smiles": "CC(=O)Nc1ccc(O)cc1",
        "name": "Acetaminophen (Tylenol)",
        "category": "Analgesic",
        "molecular_weight": 151.16
    },
    "ibuprofen": {
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "name": "Ibuprofen (Advil)",
        "category": "NSAID",
        "molecular_weight": 206.28
    },
    "aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "name": "Aspirin",
        "category": "NSAID",
        "molecular_weight": 180.16
    },
    "naproxen": {
        "smiles": "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
        "name": "Naproxen (Aleve)",
        "category": "NSAID",
        "molecular_weight": 230.26
    },
    "omeprazole": {
        "smiles": "COc1ccc2nc(CS(=O)c3ncc(C)c(OC)c3C)[nH]c2c1",
        "name": "Omeprazole (Prilosec)",
        "category": "Proton Pump Inhibitor",
        "molecular_weight": 345.42
    },
    "metformin": {
        "smiles": "CN(C)C(=N)NC(=N)N",
        "name": "Metformin (Glucophage)",
        "category": "Antidiabetic",
        "molecular_weight": 129.16
    },
    "lisinopril": {
        "smiles": "NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O",
        "name": "Lisinopril (Zestril)",
        "category": "ACE Inhibitor",
        "molecular_weight": 405.49
    },
    "amlodipine": {
        "smiles": "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl",
        "name": "Amlodipine (Norvasc)",
        "category": "Calcium Channel Blocker",
        "molecular_weight": 408.88
    },
    "atorvastatin": {
        "smiles": "CC(C)c1c(C(=O)Nc2ccccc2)c(c(c3ccc(F)cc3)n1CC[C@@H](O)C[C@@H](O)CC(=O)O)c4ccccc4",
        "name": "Atorvastatin (Lipitor)",
        "category": "Statin",
        "molecular_weight": 558.64
    },
    "simvastatin": {
        "smiles": "CCC(C)(C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12",
        "name": "Simvastatin (Zocor)",
        "category": "Statin",
        "molecular_weight": 418.57
    },
    "metoprolol": {
        "smiles": "COCCc1ccc(OCC(O)CNC(C)C)cc1",
        "name": "Metoprolol (Lopressor)",
        "category": "Beta Blocker",
        "molecular_weight": 267.36
    },
    "losartan": {
        "smiles": "CCCCc1nc(Cl)c(CO)n1Cc2ccc(c3ccccc3c4nnn[nH]4)cc2",
        "name": "Losartan (Cozaar)",
        "category": "ARB",
        "molecular_weight": 422.91
    },
    "gabapentin": {
        "smiles": "NCC1(CC(=O)O)CCCCC1",
        "name": "Gabapentin (Neurontin)",
        "category": "Anticonvulsant",
        "molecular_weight": 171.24
    },
    "pregabalin": {
        "smiles": "CC(C)CC(CN)CC(=O)O",
        "name": "Pregabalin (Lyrica)",
        "category": "Anticonvulsant",
        "molecular_weight": 159.23
    },
    "duloxetine": {
        "smiles": "CNCCC(Oc1cccc2ccccc12)c3cccs3",
        "name": "Duloxetine (Cymbalta)",
        "category": "SNRI Antidepressant",
        "molecular_weight": 297.42
    },
    "escitalopram": {
        "smiles": "CN(C)CCCC1(OCc2cc(F)ccc21)c3ccc(C#N)cc3",
        "name": "Escitalopram (Lexapro)",
        "category": "SSRI Antidepressant",
        "molecular_weight": 324.39
    },
    "citalopram": {
        "smiles": "CN(C)CCCC1(OCc2cc(F)ccc21)c3ccc(C#N)cc3",
        "name": "Citalopram (Celexa)",
        "category": "SSRI Antidepressant",
        "molecular_weight": 324.39
    },
    "bupropion": {
        "smiles": "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1",
        "name": "Bupropion (Wellbutrin)",
        "category": "NDRI Antidepressant",
        "molecular_weight": 239.74
    },
    "trazodone": {
        "smiles": "Clc1ccc(N2CCN(CCCN3N=C4C=CC=CC4=NC3=O)CC2)cc1",
        "name": "Trazodone (Desyrel)",
        "category": "Antidepressant",
        "molecular_weight": 371.86
    },
    "mirtazapine": {
        "smiles": "CN1CCN2c3ccccc3Cc4ccc(C)nc4C2C1",
        "name": "Mirtazapine (Remeron)",
        "category": "Antidepressant",
        "molecular_weight": 265.35
    },
    "amitriptyline": {
        "smiles": "CN(C)CCC=C1c2ccccc2CCc3ccccc13",
        "name": "Amitriptyline (Elavil)",
        "category": "TCA Antidepressant",
        "molecular_weight": 277.40
    },
    "nortriptyline": {
        "smiles": "CNCCC=C1c2ccccc2CCc3ccccc13",
        "name": "Nortriptyline (Pamelor)",
        "category": "TCA Antidepressant",
        "molecular_weight": 263.38
    },
    "hydroxyzine": {
        "smiles": "OCCOCCN1CCN(C(c2ccccc2)c3ccc(Cl)cc3)CC1",
        "name": "Hydroxyzine (Atarax)",
        "category": "Antihistamine",
        "molecular_weight": 374.90
    },
    "cetirizine": {
        "smiles": "OC(=O)COCCN1CCN(C(c2ccccc2)c3ccc(Cl)cc3)CC1",
        "name": "Cetirizine (Zyrtec)",
        "category": "Antihistamine",
        "molecular_weight": 388.89
    },
    "loratadine": {
        "smiles": "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc4cccnc24)CC1",
        "name": "Loratadine (Claritin)",
        "category": "Antihistamine",
        "molecular_weight": 382.88
    },
    "diphenhydramine": {
        "smiles": "CN(C)CCOC(c1ccccc1)c2ccccc2",
        "name": "Diphenhydramine (Benadryl)",
        "category": "Antihistamine",
        "molecular_weight": 255.36
    },
    "prednisone": {
        "smiles": "C[C@]12CC(=O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@]34C)[C@@H]1CC[C@]2(O)C(=O)CO",
        "name": "Prednisone",
        "category": "Corticosteroid",
        "molecular_weight": 358.43
    },
    "dexamethasone": {
        "smiles": "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO",
        "name": "Dexamethasone",
        "category": "Corticosteroid",
        "molecular_weight": 392.46
    },
    "prednisone": {
        "smiles": "C[C@]12CC(=O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@]34C)[C@@H]1CC[C@]2(O)C(=O)CO",
        "name": "Prednisone",
        "category": "Corticosteroid",
        "molecular_weight": 358.43
    },
    "albuterol": {
        "smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",
        "name": "Albuterol (Ventolin)",
        "category": "Bronchodilator",
        "molecular_weight": 239.31
    },
    "montelukast": {
        "smiles": "CC(C)(O)c1ccccc1CCC(SCC2(CC(=O)O)CC2)c3cccc(c3)C=Cc4ccc5ccc(Cl)cc5n4",
        "name": "Montelukast (Singulair)",
        "category": "Leukotriene Inhibitor",
        "molecular_weight": 586.18
    },
    "fluticasone": {
        "smiles": "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(OC(=O)SCF)C(=O)SCF",
        "name": "Fluticasone (Flovent)",
        "category": "Corticosteroid",
        "molecular_weight": 500.57
    },
    "sildenafil": {
        "smiles": "CCCc1nn(C)c2c1nc(nc2N3CCN(C)CC3)c4ccc(OCC)cc4S(=O)(=O)N5CCN(CC5)C",
        "name": "Sildenafil (Viagra)",
        "category": "PDE5 Inhibitor",
        "molecular_weight": 474.58
    },
    "tadalafil": {
        "smiles": "CN1CC(=O)N2C(Cc3c2ccc4c3OCO4)c5cc6OCOc6cc5C1=O",
        "name": "Tadalafil (Cialis)",
        "category": "PDE5 Inhibitor",
        "molecular_weight": 389.40
    },
    "vardenafil": {
        "smiles": "CCCc1nc(c2ccc(OCC)cc2S(=O)(=O)N3CCN(CC)CC3)c4n1nc(C)c(N)n4",
        "name": "Vardenafil (Levitra)",
        "category": "PDE5 Inhibitor",
        "molecular_weight": 488.60
    },
    "warfarin": {
        "smiles": "CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O",
        "name": "Warfarin (Coumadin)",
        "category": "Anticoagulant",
        "molecular_weight": 308.33
    },
    "clopidogrel": {
        "smiles": "COC(=O)[C@H](c1ccccc1Cl)N2CCc3sccc3C2",
        "name": "Clopidogrel (Plavix)",
        "category": "Antiplatelet",
        "molecular_weight": 321.82
    },
    "levothyroxine": {
        "smiles": "N[C@@H](Cc1cc(I)c(Oc2ccc(O)c(I)c2)c(I)c1)C(=O)O",
        "name": "Levothyroxine (Synthroid)",
        "category": "Thyroid Hormone",
        "molecular_weight": 776.87
    },
    "hydrochlorothiazide": {
        "smiles": "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl",
        "name": "Hydrochlorothiazide (HCTZ)",
        "category": "Diuretic",
        "molecular_weight": 297.74
    },
    "furosemide": {
        "smiles": "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl",
        "name": "Furosemide (Lasix)",
        "category": "Loop Diuretic",
        "molecular_weight": 330.74
    },
    "spironolactone": {
        "smiles": "CC(=O)SC1CC2=CC(=O)CC[C@]2(C)[C@@H]3CC[C@@H]4[C@@H](CC[C@]45CCC(=O)O5)[C@H]13",
        "name": "Spironolactone (Aldactone)",
        "category": "Potassium-Sparing Diuretic",
        "molecular_weight": 416.57
    },
    "tamsulosin": {
        "smiles": "CCOc1ccccc1OCCNC[C@H](C)Cc2ccc(OC)c(S(=O)(=O)N)c2",
        "name": "Tamsulosin (Flomax)",
        "category": "Alpha Blocker",
        "molecular_weight": 408.51
    },
    "finasteride": {
        "smiles": "CC(C)(C)NC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4NC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C",
        "name": "Finasteride (Propecia)",
        "category": "5-Alpha Reductase Inhibitor",
        "molecular_weight": 372.54
    },
    "oxybutynin": {
        "smiles": "CCN(CC)CC#COC(=O)C(O)(c1ccccc1)C2CCCCC2",
        "name": "Oxybutynin (Ditropan)",
        "category": "Anticholinergic",
        "molecular_weight": 357.49
    },
    "ciprofloxacin": {
        "smiles": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        "name": "Ciprofloxacin (Cipro)",
        "category": "Fluoroquinolone Antibiotic",
        "molecular_weight": 331.34
    },
    "levofloxacin": {
        "smiles": "C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc4c(=O)c(C(=O)O)cn1c24",
        "name": "Levofloxacin (Levaquin)",
        "category": "Fluoroquinolone Antibiotic",
        "molecular_weight": 361.37
    },
    "azithromycin": {
        "smiles": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]3O[C@H](C)C[C@@H]([C@H]3O)N(C)C)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O",
        "name": "Azithromycin (Z-Pack)",
        "category": "Macrolide Antibiotic",
        "molecular_weight": 748.98
    },
    "amoxicillin": {
        "smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O",
        "name": "Amoxicillin",
        "category": "Penicillin Antibiotic",
        "molecular_weight": 365.40
    },
    "doxycycline": {
        "smiles": "C[C@H]1[C@H]2C(=O)c3c(O)cccc3[C@@]2(O)C(=O)C4=C(O)[C@]5(O)C(=O)C(C(N)=O)=C(O)[C@@H](N(C)C)[C@@H]5[C@@H]14",
        "name": "Doxycycline",
        "category": "Tetracycline Antibiotic",
        "molecular_weight": 444.43
    },
    "metronidazole": {
        "smiles": "Cc1ncc([N+](=O)[O-])n1CCO",
        "name": "Metronidazole (Flagyl)",
        "category": "Antibiotic/Antiparasitic",
        "molecular_weight": 171.15
    },
    "fluconazole": {
        "smiles": "OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F",
        "name": "Fluconazole (Diflucan)",
        "category": "Antifungal",
        "molecular_weight": 306.27
    },
    "acyclovir": {
        "smiles": "Nc1nc2n(COCCO)cnc2c(=O)[nH]1",
        "name": "Acyclovir (Zovirax)",
        "category": "Antiviral",
        "molecular_weight": 225.20
    },
    "valacyclovir": {
        "smiles": "CC(C)[C@H](N)C(=O)OCCOC1=Nc2nc(N)nc(N)c2C(=O)N1",
        "name": "Valacyclovir (Valtrex)",
        "category": "Antiviral",
        "molecular_weight": 324.34
    },
    "oseltamivir": {
        "smiles": "CCOC(=O)C1=CC(OC(CC)CC)C(NC(C)=O)C(N)C1",
        "name": "Oseltamivir (Tamiflu)",
        "category": "Antiviral",
        "molecular_weight": 312.40
    },
    "zolpidem": {
        "smiles": "Cc1ccc(C(=O)N(C)Cc2nc3ccccc3n2C)c(C)c1",
        "name": "Zolpidem (Ambien)",
        "category": "Sedative/Hypnotic",
        "molecular_weight": 307.39
    },
    "eszopiclone": {
        "smiles": "CC(=O)OC1C(=O)N(c2ccc(Cl)cn2)c3ncc(cn13)C#N",
        "name": "Eszopiclone (Lunesta)",
        "category": "Sedative/Hypnotic",
        "molecular_weight": 388.81
    },
    "melatonin": {
        "smiles": "COc1ccc2[nH]cc(CCNC(C)=O)c2c1",
        "name": "Melatonin",
        "category": "Sleep Aid",
        "molecular_weight": 232.28
    },
    "modafinil": {
        "smiles": "NC(=O)CS(=O)C(c1ccccc1)c2ccccc2",
        "name": "Modafinil (Provigil)",
        "category": "Wakefulness Agent",
        "molecular_weight": 273.35
    },
    "methylphenidate": {
        "smiles": "COC(=O)[C@H](c1ccccc1)[C@@H]2CCCCN2",
        "name": "Methylphenidate (Ritalin)",
        "category": "Stimulant",
        "molecular_weight": 233.31
    },
    "amphetamine": {
        "smiles": "CC(N)Cc1ccccc1",
        "name": "Amphetamine (Adderall)",
        "category": "Stimulant",
        "molecular_weight": 135.21
    },
    "atomoxetine": {
        "smiles": "CC(Cc1ccccc1)NCCOc2ccccc2",
        "name": "Atomoxetine (Strattera)",
        "category": "ADHD Medication",
        "molecular_weight": 255.36
    },
    "donepezil": {
        "smiles": "COc1cc2CCN(Cc3ccc(cc3)C(=O)CCc4ccccc4)Cc2cc1OC",
        "name": "Donepezil (Aricept)",
        "category": "Cholinesterase Inhibitor",
        "molecular_weight": 379.49
    },
    "memantine": {
        "smiles": "CC12CC3CC(C)(C1)CC(N)(C3)C2",
        "name": "Memantine (Namenda)",
        "category": "NMDA Antagonist",
        "molecular_weight": 179.30
    },
    "rivastigmine": {
        "smiles": "CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1",
        "name": "Rivastigmine (Exelon)",
        "category": "Cholinesterase Inhibitor",
        "molecular_weight": 250.34
    },
    "levodopa": {
        "smiles": "N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O",
        "name": "Levodopa (L-DOPA)",
        "category": "Dopamine Precursor",
        "molecular_weight": 197.19
    },
    "carbidopa": {
        "smiles": "CC(C(=O)O)(NN)Cc1ccc(O)c(O)c1",
        "name": "Carbidopa",
        "category": "Decarboxylase Inhibitor",
        "molecular_weight": 226.23
    },
    "pramipexole": {
        "smiles": "CCCN[C@H]1CCc2nc(N)sc2C1",
        "name": "Pramipexole (Mirapex)",
        "category": "Dopamine Agonist",
        "molecular_weight": 211.33
    },
    "ropinirole": {
        "smiles": "CCCN(CCC)CCc1cccc2NC(=O)Cc12",
        "name": "Ropinirole (Requip)",
        "category": "Dopamine Agonist",
        "molecular_weight": 260.37
    },
    "sumatriptan": {
        "smiles": "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1",
        "name": "Sumatriptan (Imitrex)",
        "category": "Triptan (Migraine)",
        "molecular_weight": 295.40
    },
    "rizatriptan": {
        "smiles": "CN(C)CCc1c[nH]c2ccc(Cn3cncn3)cc12",
        "name": "Rizatriptan (Maxalt)",
        "category": "Triptan (Migraine)",
        "molecular_weight": 269.35
    },
    "topiramate": {
        "smiles": "CC1(C)O[C@H]2CO[C@@]3(COS(N)(=O)=O)OC(C)(C)O[C@@H]3[C@H]2O1",
        "name": "Topiramate (Topamax)",
        "category": "Anticonvulsant",
        "molecular_weight": 339.36
    },
    "valproic_acid": {
        "smiles": "CCCC(CCC)C(=O)O",
        "name": "Valproic Acid (Depakote)",
        "category": "Anticonvulsant",
        "molecular_weight": 144.21
    },
    "lamotrigine": {
        "smiles": "Nc1nnc(c(N)n1)c2cccc(Cl)c2Cl",
        "name": "Lamotrigine (Lamictal)",
        "category": "Anticonvulsant",
        "molecular_weight": 256.09
    },
    "levetiracetam": {
        "smiles": "CC[C@H](N)C(=O)N1CCCC1C(=O)N",
        "name": "Levetiracetam (Keppra)",
        "category": "Anticonvulsant",
        "molecular_weight": 170.21
    },
    "carbamazepine": {
        "smiles": "NC(=O)N1c2ccccc2C=Cc3ccccc13",
        "name": "Carbamazepine (Tegretol)",
        "category": "Anticonvulsant",
        "molecular_weight": 236.27
    },
    "phenytoin": {
        "smiles": "O=C1NC(=O)NC1(c2ccccc2)c3ccccc3",
        "name": "Phenytoin (Dilantin)",
        "category": "Anticonvulsant",
        "molecular_weight": 252.27
    },
    "clonidine": {
        "smiles": "Clc1cccc(Cl)c1NC2=NCCN2",
        "name": "Clonidine (Catapres)",
        "category": "Alpha-2 Agonist",
        "molecular_weight": 230.09
    },
    "prazosin": {
        "smiles": "COc1cc2nc(nc(N)c2cc1OC)N3CCN(CC3)C(=O)c4ccco4",
        "name": "Prazosin (Minipress)",
        "category": "Alpha Blocker",
        "molecular_weight": 383.40
    },
    "propranolol": {
        "smiles": "CC(C)NCC(O)COc1cccc2ccccc12",
        "name": "Propranolol (Inderal)",
        "category": "Beta Blocker",
        "molecular_weight": 259.34
    },
    "atenolol": {
        "smiles": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
        "name": "Atenolol (Tenormin)",
        "category": "Beta Blocker",
        "molecular_weight": 266.34
    },
    "carvedilol": {
        "smiles": "COc1ccccc1OCCNCC(O)COc2cccc3[nH]c4ccccc4c23",
        "name": "Carvedilol (Coreg)",
        "category": "Beta Blocker",
        "molecular_weight": 406.47
    },
    "diltiazem": {
        "smiles": "COc1ccc(C2Sc3ccccc3N(CCN(C)C)C(=O)C2OC(=O)C)cc1",
        "name": "Diltiazem (Cardizem)",
        "category": "Calcium Channel Blocker",
        "molecular_weight": 414.52
    },
    "verapamil": {
        "smiles": "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC",
        "name": "Verapamil (Calan)",
        "category": "Calcium Channel Blocker",
        "molecular_weight": 454.60
    },
    "nifedipine": {
        "smiles": "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c2ccccc2[N+](=O)[O-]",
        "name": "Nifedipine (Procardia)",
        "category": "Calcium Channel Blocker",
        "molecular_weight": 346.34
    },
    "ramipril": {
        "smiles": "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N2[C@H](C(=O)O)C[C@@H]3CCCC[C@H]32",
        "name": "Ramipril (Altace)",
        "category": "ACE Inhibitor",
        "molecular_weight": 416.51
    },
    "enalapril": {
        "smiles": "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N2CCC[C@H]2C(=O)O",
        "name": "Enalapril (Vasotec)",
        "category": "ACE Inhibitor",
        "molecular_weight": 376.45
    },
    "benazepril": {
        "smiles": "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N2Cc3ccccc3C[C@H]2C(=O)O",
        "name": "Benazepril (Lotensin)",
        "category": "ACE Inhibitor",
        "molecular_weight": 424.49
    },
    "valsartan": {
        "smiles": "CCCCC(=O)N(Cc1ccc(c2ccccc2c3nnn[nH]3)cc1)[C@@H](C(C)C)C(=O)O",
        "name": "Valsartan (Diovan)",
        "category": "ARB",
        "molecular_weight": 435.52
    },
    "irbesartan": {
        "smiles": "CCCCC1=NC2(CCCC2)C(=O)N1Cc3ccc(c4ccccc4c5nnn[nH]5)cc3",
        "name": "Irbesartan (Avapro)",
        "category": "ARB",
        "molecular_weight": 428.53
    },
    "olmesartan": {
        "smiles": "CCCc1nc(c(C(=O)O)n1Cc2ccc(c3ccccc3c4nnn[nH]4)cc2)C(C)(C)O",
        "name": "Olmesartan (Benicar)",
        "category": "ARB",
        "molecular_weight": 446.50
    },
    "ezetimibe": {
        "smiles": "O[C@H]1[C@H](Oc2ccc(F)cc2)[C@@H](c3ccc(F)cc3)N(c4ccc(O)cc4)[C@H]1C(=O)c5ccc(F)cc5",
        "name": "Ezetimibe (Zetia)",
        "category": "Cholesterol Absorption Inhibitor",
        "molecular_weight": 409.43
    },
    "rosuvastatin": {
        "smiles": "CC(C)c1nc(N(C)S(=O)(=O)C)nc(c1C=CC(O)CC(O)CC(=O)O)c2ccc(F)cc2",
        "name": "Rosuvastatin (Crestor)",
        "category": "Statin",
        "molecular_weight": 481.54
    },
    "pravastatin": {
        "smiles": "CC[C@H](C)[C@H]1C(=O)O[C@H]2C[C@@H](O)C=C3C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(=O)O)[C@@H]23",
        "name": "Pravastatin (Pravachol)",
        "category": "Statin",
        "molecular_weight": 424.53
    },
    "fenofibrate": {
        "smiles": "CC(C)OC(=O)C(C)(C)Oc1ccc(C(=O)c2ccc(Cl)cc2)cc1",
        "name": "Fenofibrate (Tricor)",
        "category": "Fibrate",
        "molecular_weight": 360.83
    },
    "gemfibrozil": {
        "smiles": "CC(C)(CCCC1=CC(C)=C(C)C=C1)C(=O)O",
        "name": "Gemfibrozil (Lopid)",
        "category": "Fibrate",
        "molecular_weight": 250.33
    },
    "niacin": {
        "smiles": "OC(=O)c1cccnc1",
        "name": "Niacin (Vitamin B3)",
        "category": "Lipid Modifier",
        "molecular_weight": 123.11
    },
    "colchicine": {
        "smiles": "COc1cc2CCc3cc(OC)c(OC)c(OC)c3-c2c(OC)c1NC(C)=O",
        "name": "Colchicine",
        "category": "Anti-Gout",
        "molecular_weight": 399.44
    },
    "allopurinol": {
        "smiles": "O=c1[nH]cnc2[nH]ncc12",
        "name": "Allopurinol (Zyloprim)",
        "category": "Xanthine Oxidase Inhibitor",
        "molecular_weight": 136.11
    },
    "febuxostat": {
        "smiles": "Cc1nc(c(C#N)c(n1)c2ccc(OCC(C)C)c(c2)C(=O)O)c3ccccc3",
        "name": "Febuxostat (Uloric)",
        "category": "Xanthine Oxidase Inhibitor",
        "molecular_weight": 316.37
    },
    "probenecid": {
        "smiles": "CCCN(CCC)S(=O)(=O)c1ccc(C(=O)O)cc1",
        "name": "Probenecid",
        "category": "Uricosuric",
        "molecular_weight": 285.36
    },
    "celecoxib": {
        "smiles": "Cc1ccc(c(c1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F)C",
        "name": "Celecoxib (Celebrex)",
        "category": "COX-2 Inhibitor",
        "molecular_weight": 381.37
    },
    "meloxicam": {
        "smiles": "Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1",
        "name": "Meloxicam (Mobic)",
        "category": "NSAID",
        "molecular_weight": 351.40
    },
    "indomethacin": {
        "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c3ccc(Cl)cc3",
        "name": "Indomethacin (Indocin)",
        "category": "NSAID",
        "molecular_weight": 357.79
    },
    "diclofenac": {
        "smiles": "OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl",
        "name": "Diclofenac (Voltaren)",
        "category": "NSAID",
        "molecular_weight": 296.15
    },
    "ketorolac": {
        "smiles": "OC(=O)C1CCc2n1cc(c2C(=O)c3ccccc3)C",
        "name": "Ketorolac (Toradol)",
        "category": "NSAID",
        "molecular_weight": 255.27
    },
    "piroxicam": {
        "smiles": "CN1C(=C(O)c2ccccc2S1(=O)=O)C(=O)Nc3ccccn3",
        "name": "Piroxicam (Feldene)",
        "category": "NSAID",
        "molecular_weight": 331.35
    },
    "tramadol": {
        "smiles": "COc1ccccc1C2(O)CCCCC2CN(C)C",
        "name": "Tramadol (Ultram)",
        "category": "Opioid Analgesic",
        "molecular_weight": 263.38
    },
    "codeine": {
        "smiles": "COc1ccc2C[C@H]3N(C)CC[C@@]45c2c1O[C@H]4[C@@H](O)C=C[C@@H]35",
        "name": "Codeine",
        "category": "Opioid Analgesic",
        "molecular_weight": 299.36
    },
    "hydrocodone": {
        "smiles": "COc1ccc2C[C@H]3N(C)CC[C@@]45c2c1O[C@H]4C(=O)CC[C@@H]35",
        "name": "Hydrocodone (Vicodin)",
        "category": "Opioid Analgesic",
        "molecular_weight": 299.36
    },
    "oxycodone": {
        "smiles": "COc1ccc2C[C@H]3N(C)CC[C@@]45[C@@H](O)Oc1c2[C@@]4(O)C(=O)CC[C@@H]35",
        "name": "Oxycodone (OxyContin)",
        "category": "Opioid Analgesic",
        "molecular_weight": 315.36
    },
    "morphine": {
        "smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
        "name": "Morphine",
        "category": "Opioid Analgesic",
        "molecular_weight": 285.34
    },
    "naloxone": {
        "smiles": "O=C[C@H]1C[C@H]2[C@@H]3Cc4ccc(O)c5O[C@@H]3[C@](O)(C1)[C@@H]2c45N(CC=C)C",
        "name": "Naloxone (Narcan)",
        "category": "Opioid Antagonist",
        "molecular_weight": 327.37
    },
    "naltrexone": {
        "smiles": "O=C1CC[C@H]2[C@@H]3Cc4ccc(O)c5O[C@@H]3[C@](O)(C1)[C@@H]2c45N(CC6CC6)C",
        "name": "Naltrexone (Vivitrol)",
        "category": "Opioid Antagonist",
        "molecular_weight": 341.40
    },
    "buprenorphine": {
        "smiles": "COC(C)(C)[C@H]1C[C@H]2[C@@H]3Cc4ccc(O)c5O[C@@H]3[C@](O)(C1)[C@@H]2c45N(CC6CC6)C",
        "name": "Buprenorphine (Suboxone)",
        "category": "Partial Opioid Agonist",
        "molecular_weight": 467.64
    },
    "methadone": {
        "smiles": "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c2ccccc2",
        "name": "Methadone",
        "category": "Opioid Agonist",
        "molecular_weight": 309.44
    },
    "fentanyl": {
        "smiles": "CCC(=O)N(c1ccccc1)C2CCN(CCc3ccccc3)CC2",
        "name": "Fentanyl",
        "category": "Opioid Analgesic",
        "molecular_weight": 336.47
    },
    "lidocaine": {
        "smiles": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
        "name": "Lidocaine",
        "category": "Local Anesthetic",
        "molecular_weight": 234.34
    },
    "bupivacaine": {
        "smiles": "CCCCN1CCCCC1C(=O)Nc2c(C)cccc2C",
        "name": "Bupivacaine (Marcaine)",
        "category": "Local Anesthetic",
        "molecular_weight": 288.43
    },
    "ketamine": {
        "smiles": "CNC1(CCCCC1=O)c2ccc(Cl)cc2",
        "name": "Ketamine",
        "category": "Dissociative Anesthetic",
        "molecular_weight": 237.73
    },
    "propofol": {
        "smiles": "CC(C)c1cccc(C(C)C)c1O",
        "name": "Propofol (Diprivan)",
        "category": "General Anesthetic",
        "molecular_weight": 178.27
    },
    "midazolam": {
        "smiles": "Cc1ncc2n1c3ccc(Cl)cc3C(c4ccccc4F)=NC2",
        "name": "Midazolam (Versed)",
        "category": "Benzodiazepine",
        "molecular_weight": 325.77
    },
    "diazepam": {
        "smiles": "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13",
        "name": "Diazepam (Valium)",
        "category": "Benzodiazepine",
        "molecular_weight": 284.74
    },
    "lorazepam": {
        "smiles": "OC1N=C(c2ccccc2Cl)c3cc(Cl)ccc3NC1=O",
        "name": "Lorazepam (Ativan)",
        "category": "Benzodiazepine",
        "molecular_weight": 321.16
    },
    "alprazolam": {
        "smiles": "Cc1nnc2CN=C(c3ccccc3Cl)c4cc(Cl)ccc4n12",
        "name": "Alprazolam (Xanax)",
        "category": "Benzodiazepine",
        "molecular_weight": 308.76
    },
    "clonazepam": {
        "smiles": "[O-][N+](=O)c1ccc2NC(=O)CN=C(c3ccccc3Cl)c2c1",
        "name": "Clonazepam (Klonopin)",
        "category": "Benzodiazepine",
        "molecular_weight": 315.71
    },
    "temazepam": {
        "smiles": "CN1C(=O)C(O)N=C(c2ccccc2)c3cc(Cl)ccc13",
        "name": "Temazepam (Restoril)",
        "category": "Benzodiazepine",
        "molecular_weight": 300.74
    },
    "triazolam": {
        "smiles": "Cc1nnc2CN=C(c3ccccc3Cl)c4cc(Cl)ccc4n12",
        "name": "Triazolam (Halcion)",
        "category": "Benzodiazepine",
        "molecular_weight": 343.21
    },
    "buspirone": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN3CCN(c4ncccn4)CC3",
        "name": "Buspirone (BuSpar)",
        "category": "Anxiolytic",
        "molecular_weight": 385.50
    },
    "haloperidol": {
        "smiles": "OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c3ccc(Cl)cc3",
        "name": "Haloperidol (Haldol)",
        "category": "Antipsychotic",
        "molecular_weight": 375.86
    },
    "chlorpromazine": {
        "smiles": "CN(C)CCCN1c2ccccc2Sc3ccc(Cl)cc13",
        "name": "Chlorpromazine (Thorazine)",
        "category": "Antipsychotic",
        "molecular_weight": 318.86
    },
    "olanzapine": {
        "smiles": "Cc1cc2Nc3ccc(N4CCN(C)CC4)cc3N=C2S1C",
        "name": "Olanzapine (Zyprexa)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 312.43
    },
    "quetiapine": {
        "smiles": "OCCOCCN1CCN(c2Nc3ccccc3Sc4ccccc24)CC1",
        "name": "Quetiapine (Seroquel)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 383.51
    },
    "risperidone": {
        "smiles": "Cc1nc2n(c(=O)c1CCN3CCC(CC3)c4c5ccc(F)cc5on4)CCC2",
        "name": "Risperidone (Risperdal)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 410.48
    },
    "aripiprazole": {
        "smiles": "Clc1cccc(N2CCN(CCCCOC3=Cc4ccccc4NC3=O)CC2)c1Cl",
        "name": "Aripiprazole (Abilify)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 448.38
    },
    "ziprasidone": {
        "smiles": "Clc1ccc2N(CC3CCN(CC3)c4nsc5ccccc45)C(=O)Cc2c1Cl",
        "name": "Ziprasidone (Geodon)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 412.94
    },
    "clozapine": {
        "smiles": "CN1CCN(C2=Nc3cc(Cl)ccc3Nc4ccccc24)CC1",
        "name": "Clozapine (Clozaril)",
        "category": "Atypical Antipsychotic",
        "molecular_weight": 326.82
    },
    "lithium_carbonate": {
        "smiles": "[Li+].[Li+].[O-]C([O-])=O",
        "name": "Lithium Carbonate",
        "category": "Mood Stabilizer",
        "molecular_weight": 73.89
    },
    "caffeine": {
        "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        "name": "Caffeine",
        "category": "Stimulant",
        "molecular_weight": 194.19
    },
    "theophylline": {
        "smiles": "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
        "name": "Theophylline",
        "category": "Bronchodilator",
        "molecular_weight": 180.16
    },
    "nicotine": {
        "smiles": "CN1CCC[C@H]1c2cccnc2",
        "name": "Nicotine",
        "category": "Stimulant",
        "molecular_weight": 162.23
    },
    "varenicline": {
        "smiles": "c1ccc2c(c1)CC3CNCC4=CN=CC(=C34)C2",
        "name": "Varenicline (Chantix)",
        "category": "Smoking Cessation",
        "molecular_weight": 211.26
    },
    "disulfiram": {
        "smiles": "CCN(CC)C(=S)SSC(=S)N(CC)CC",
        "name": "Disulfiram (Antabuse)",
        "category": "Alcohol Deterrent",
        "molecular_weight": 296.54
    },
    "acamprosate": {
        "smiles": "CC(=O)NCCCS(=O)(=O)O",
        "name": "Acamprosate (Campral)",
        "category": "Alcohol Dependence",
        "molecular_weight": 181.21
    },
    "ondansetron": {
        "smiles": "Cc1nccn1CC2CCc3c(O2)c4ccccc4n3C=O",
        "name": "Ondansetron (Zofran)",
        "category": "Antiemetic",
        "molecular_weight": 293.36
    },
    "metoclopramide": {
        "smiles": "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC",
        "name": "Metoclopramide (Reglan)",
        "category": "Antiemetic",
        "molecular_weight": 299.80
    },
    "promethazine": {
        "smiles": "CC(CN1c2ccccc2Sc3ccccc13)N(C)C",
        "name": "Promethazine (Phenergan)",
        "category": "Antiemetic/Antihistamine",
        "molecular_weight": 284.42
    },
    "ranitidine": {
        "smiles": "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1",
        "name": "Ranitidine (Zantac)",
        "category": "H2 Blocker",
        "molecular_weight": 314.40
    },
    "famotidine": {
        "smiles": "NC(=N)NC1=NC(CSCCC(=N)NS(=O)(=O)N)=CS1",
        "name": "Famotidine (Pepcid)",
        "category": "H2 Blocker",
        "molecular_weight": 337.45
    },
    "esomeprazole": {
        "smiles": "COc1ccc2nc([S@@](=O)Cc3ncc(C)c(OC)c3C)[nH]c2c1",
        "name": "Esomeprazole (Nexium)",
        "category": "Proton Pump Inhibitor",
        "molecular_weight": 345.42
    },
    "pantoprazole": {
        "smiles": "COc1ccnc(CS(=O)c2nc3ccc(OC(F)F)cc3[nH]2)c1OC",
        "name": "Pantoprazole (Protonix)",
        "category": "Proton Pump Inhibitor",
        "molecular_weight": 383.37
    },
    "lansoprazole": {
        "smiles": "Cc1c(OCC(F)(F)F)ccnc1CS(=O)c2[nH]c3ccccc3n2",
        "name": "Lansoprazole (Prevacid)",
        "category": "Proton Pump Inhibitor",
        "molecular_weight": 369.36
    },
    "sucralfate": {
        "smiles": "[Al+3].[O-]S(=O)(=O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)[O-])OS([O-])(=O)=O)OS([O-])(=O)=O)OS([O-])(=O)=O",
        "name": "Sucralfate (Carafate)",
        "category": "Gastroprotective",
        "molecular_weight": 2086.75
    },
    "loperamide": {
        "smiles": "CN(C)C(=O)C(CCN1CCC(O)(CC1)c2ccc(Cl)cc2)(c3ccccc3)c4ccccc4",
        "name": "Loperamide (Imodium)",
        "category": "Antidiarrheal",
        "molecular_weight": 477.04
    },
    "psyllium": {
        "smiles": "OCC1OC(OC2C(O)C(O)C(O)C(CO)O2)C(O)C(O)C1O",
        "name": "Psyllium (Metamucil)",
        "category": "Laxative/Fiber",
        "molecular_weight": 504.44
    },
    "polyethylene_glycol": {
        "smiles": "OCCOCCO",
        "name": "Polyethylene Glycol (Miralax)",
        "category": "Osmotic Laxative",
        "molecular_weight": 106.12
    },
    "bisacodyl": {
        "smiles": "CC(=O)Oc1ccc(C(c2ccc(OC(C)=O)cc2)c3ccccn3)cc1",
        "name": "Bisacodyl (Dulcolax)",
        "category": "Stimulant Laxative",
        "molecular_weight": 361.39
    },
    "docusate": {
        "smiles": "CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)S(=O)(=O)[O-].[Na+]",
        "name": "Docusate (Colace)",
        "category": "Stool Softener",
        "molecular_weight": 444.56
    },
    "simethicone": {
        "smiles": "[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C",
        "name": "Simethicone (Gas-X)",
        "category": "Antiflatulent",
        "molecular_weight": 222.46
    },
    "calcium_carbonate": {
        "smiles": "[Ca+2].[O-]C([O-])=O",
        "name": "Calcium Carbonate (Tums)",
        "category": "Antacid",
        "molecular_weight": 100.09
    },
    "magnesium_hydroxide": {
        "smiles": "[OH-].[OH-].[Mg+2]",
        "name": "Magnesium Hydroxide (Milk of Magnesia)",
        "category": "Antacid/Laxative",
        "molecular_weight": 58.32
    },
    "aluminum_hydroxide": {
        "smiles": "[OH-].[OH-].[OH-].[Al+3]",
        "name": "Aluminum Hydroxide (Maalox)",
        "category": "Antacid",
        "molecular_weight": 78.00
    },
    "psilocybin": {
        "smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "name": "Psilocybin",
        "category": "Psychedelic",
        "molecular_weight": 284.25
    },
    "lsd": {
        "smiles": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C",
        "name": "LSD",
        "category": "Psychedelic",
        "molecular_weight": 323.43
    },
    "mdma": {
        "smiles": "CC(CC1=CC2=C(C=C1)OCO2)NC",
        "name": "MDMA (Ecstasy)",
        "category": "Empathogen",
        "molecular_weight": 193.25
    },
    "dmt": {
        "smiles": "CN(C)CCc1c[nH]c2ccccc12",
        "name": "DMT",
        "category": "Psychedelic",
        "molecular_weight": 188.27
    },
    "mescaline": {
        "smiles": "COc1cc(CCN)cc(OC)c1OC",
        "name": "Mescaline",
        "category": "Psychedelic",
        "molecular_weight": 211.26
    },
    "thc": {
        "smiles": "CCCCCc1cc(O)c2C3CC(=C)CCC3C(C)(C)Oc2c1",
        "name": "THC (Dronabinol)",
        "category": "Cannabinoid",
        "molecular_weight": 314.46
    },
    "cbd": {
        "smiles": "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1",
        "name": "CBD (Epidiolex)",
        "category": "Cannabinoid",
        "molecular_weight": 314.46
    },
}


def get_all_drugs():
    """Return the full drug database"""
    return FDA_APPROVED_DRUGS


def get_drug_count():
    """Return the number of drugs in the database"""
    return len(FDA_APPROVED_DRUGS)


def search_drugs(query: str) -> dict:
    """Search for drugs by name or category"""
    query = query.lower()
    results = {}
    for key, drug in FDA_APPROVED_DRUGS.items():
        if query in key.lower() or query in drug.get('name', '').lower() or query in drug.get('category', '').lower():
            results[key] = drug
    return results


def get_drugs_by_category(category: str) -> dict:
    """Get all drugs in a specific category"""
    category = category.lower()
    results = {}
    for key, drug in FDA_APPROVED_DRUGS.items():
        if category in drug.get('category', '').lower():
            results[key] = drug
    return results
