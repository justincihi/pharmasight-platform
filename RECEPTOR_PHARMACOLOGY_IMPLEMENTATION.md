# Receptor Pharmacology Module - Implementation Summary

**Date:** October 27, 2025  
**Module:** `src/receptor_pharmacology.py`  
**Status:** ✅ Complete and Tested

---

## Overview

Comprehensive receptor targeting and pharmacology module with **70 receptor subtypes** across 14 neurotransmitter families, supporting detailed pharmacological profiling including modulation types and signaling bias.

---

## Receptor Families Implemented

### 1. Serotonin Receptors (11 subtypes)
- **5-HT1A** - Anxiolytic, antidepressant (Gi/o-coupled)
- **5-HT1B** - Mood regulation, aggression (Gi/o-coupled)
- **5-HT1D** - Migraine treatment (Gi/o-coupled)
- **5-HT2A** - Psychedelic effects, mood, cognition (Gq/11-coupled) ⭐
- **5-HT2B** - Cardiac valve regulation (Gq/11-coupled)
- **5-HT2C** - Appetite, mood, cognition (Gq/11-coupled) ⭐
- **5-HT3** - Nausea, vomiting (Ion channel)
- **5-HT4** - GI motility, cognition (Gs-coupled)
- **5-HT5A** - Orphan receptor (Gi/o-coupled)
- **5-HT6** - Cognition, memory (Gs-coupled)
- **5-HT7** - Circadian rhythm, mood (Gs-coupled)

### 2. Dopamine Receptors (5 subtypes)
- **D1** - Motor control, reward, cognition (Gs-coupled)
- **D2** - Motor control, prolactin regulation (Gi/o-coupled)
- **D3** - Reward, addiction (Gi/o-coupled)
- **D4** - Cognition, attention (Gi/o-coupled)
- **D5** - Cognition, renal function (Gs-coupled)

### 3. GABA Receptors (7 subtypes)
- **GABAA-α1** - Sedation, amnesia (Cl- channel) ⭐
- **GABAA-α2** - Anxiolytic effects (Cl- channel) ⭐
- **GABAA-α3** - Anxiolytic, muscle relaxation (Cl- channel) ⭐
- **GABAA-α4** - Tonic inhibition (Cl- channel)
- **GABAA-α5** - Cognition, memory (Cl- channel)
- **GABAA-α6** - Motor coordination (Cl- channel)
- **GABAB** - Muscle relaxation, anxiolytic (Gi/o-coupled)

### 4. Glutamate Receptors (13 subtypes)
**NMDA Receptors:**
- **NMDA-GluN1** - Essential subunit (Ca2+/Na+ channel)
- **NMDA-GluN2A** - Synaptic plasticity (Ca2+/Na+ channel)
- **NMDA-GluN2B** - Pain, depression, ketamine target (Ca2+/Na+ channel) ⭐
- **NMDA-GluN2C** - Motor coordination (Ca2+/Na+ channel)
- **NMDA-GluN2D** - Sensory processing (Ca2+/Na+ channel)

**AMPA Receptors:**
- **AMPA-GluA1** - Fast excitatory transmission (Na+/K+ channel)
- **AMPA-GluA2** - Fast excitatory transmission (Na+/K+ channel)
- **AMPA-GluA3** - Synaptic plasticity (Na+/K+ channel)
- **AMPA-GluA4** - Synaptic plasticity (Na+/K+ channel)

**Kainate Receptors:**
- **Kainate-GluK1** - Modulation of transmission (Na+/K+ channel)
- **Kainate-GluK2** - Modulation of transmission (Na+/K+ channel)

**Metabotropic Receptors:**
- **mGluR1** - Synaptic plasticity (Gq/11-coupled)
- **mGluR5** - Addiction, pain, anxiety (Gq/11-coupled)

### 5. Opioid Receptors (4 subtypes)
- **μ (Mu)** - Analgesia, euphoria, respiratory depression (Gi/o-coupled) ⭐
- **δ (Delta)** - Analgesia, antidepressant (Gi/o-coupled) ⭐
- **κ (Kappa)** - Analgesia, dysphoria (Gi/o-coupled) ⭐
- **NOP** - Pain modulation, anxiety (Gi/o-coupled)

### 6. Cannabinoid Receptors (2 subtypes)
- **CB1** - Psychoactive effects, appetite, pain (Gi/o-coupled)
- **CB2** - Immune modulation, inflammation (Gi/o-coupled)

### 7. Adrenergic Receptors (9 subtypes)
**Alpha Receptors:**
- **α1A** - Vasoconstriction, urinary retention (Gq/11-coupled)
- **α1B** - Vasoconstriction (Gq/11-coupled)
- **α1D** - Vasoconstriction (Gq/11-coupled)
- **α2A** - Sedation, analgesia, hypotension (Gi/o-coupled)
- **α2B** - Vasoconstriction (Gi/o-coupled)
- **α2C** - Mood regulation (Gi/o-coupled)

**Beta Receptors:**
- **β1** - Increased heart rate and contractility (Gs-coupled)
- **β2** - Bronchodilation, vasodilation (Gs-coupled)
- **β3** - Lipolysis, bladder relaxation (Gs-coupled)

### 8. Muscarinic Receptors (5 subtypes)
- **M1** - Cognition, memory (Gq/11-coupled)
- **M2** - Decreased heart rate (Gi/o-coupled)
- **M3** - Glandular secretion, smooth muscle (Gq/11-coupled)
- **M4** - Motor control, cognition (Gi/o-coupled)
- **M5** - Dopamine release modulation (Gq/11-coupled)

### 9. Nicotinic Receptors (3 subtypes)
- **α4β2** - Cognition, addiction (Na+/K+ channel)
- **α7** - Cognition, neuroprotection (Na+/K+ channel)
- **α3β4** - Autonomic transmission (Na+/K+ channel)

### 10. Orexin Receptors (2 subtypes)
- **OX1** - Wakefulness, arousal (Gq/11-coupled)
- **OX2** - Sleep-wake regulation (Gq/11-coupled)

### 11. mTOR Pathway (2 complexes) ⭐
- **mTORC1** - Protein synthesis, neuroplasticity, synaptic growth
- **mTORC2** - Cell survival, cytoskeletal organization

### 12. Histamine Receptors (4 subtypes)
- **H1** - Wakefulness, allergic response (Gq/11-coupled)
- **H2** - Gastric acid secretion (Gs-coupled)
- **H3** - Neurotransmitter release modulation (Gi/o-coupled)
- **H4** - Immune modulation (Gi/o-coupled)

### 13. Sigma Receptors (2 subtypes)
- **σ1** - Neuroprotection, ion channel modulation
- **σ2** - Cell proliferation, apoptosis

### 14. Trace Amine Receptors (1 subtype)
- **TAAR1** - Modulation of monoaminergic systems (Gs-coupled)

---

## Pharmacological Modulation Types

### Agonism/Antagonism:
1. **Full Agonist** - 100% receptor activation
2. **Partial Agonist** - Submaximal receptor activation (e.g., 50-80%)
3. **Full Antagonist** - Blocks receptor activation
4. **Partial Antagonist** - Partial blockade
5. **Inverse Agonist** - Reduces basal receptor activity

### Allosteric Modulation:
6. **Positive Allosteric Modulator (PAM)** - Enhances agonist effects
7. **Negative Allosteric Modulator (NAM)** - Reduces agonist effects
8. **Silent Allosteric Modulator (SAM)** - Binds without effect

---

## Signaling Bias Options

### G-Protein vs β-Arrestin Pathway:
1. **G-Protein Biased** - Preferential G-protein signaling
2. **β-Arrestin Biased** - Preferential β-arrestin signaling
3. **Balanced** - No signaling bias
4. **Unknown** - Bias not determined

**Clinical Relevance:**
- G-protein biased opioids may have reduced respiratory depression
- β-arrestin biased compounds may have different side effect profiles
- Important for optimizing therapeutic index

---

## Key Features

### ReceptorDatabase Class:
- ✅ 70 receptor subtypes with complete metadata
- ✅ Receptor family organization
- ✅ Subtype search functionality
- ✅ Location and function information
- ✅ Signaling pathway details

### ReceptorProfile Class:
- ✅ Compound-specific receptor binding profiles
- ✅ Ki (binding affinity) tracking in nM
- ✅ Efficacy percentage (0-100%)
- ✅ Modulation type specification
- ✅ Signaling bias annotation
- ✅ Primary vs off-target classification
- ✅ JSON export capability

### Enums:
- ✅ ModulationType - 8 modulation types
- ✅ SignalingBias - 4 bias categories

---

## Usage Examples

### 1. Get All Receptor Families:
```python
from src.receptor_pharmacology import ReceptorDatabase

db = ReceptorDatabase()
receptors = db.get_all_receptors()

for family_key, family_data in receptors.items():
    print(f"{family_data['name']}: {len(family_data['subtypes'])} subtypes")
```

### 2. Search for Specific Receptors:
```python
# Search for NMDA receptors
results = db.search_receptors("NMDA")

# Search by function
results = db.search_receptors("depression")
```

### 3. Create Receptor Profile:
```python
from src.receptor_pharmacology import ReceptorProfile, ModulationType, SignalingBias

profile = ReceptorProfile("Psilocybin")

profile.add_interaction(
    receptor_family="serotonin",
    receptor_subtype="5-HT2A",
    modulation_type=ModulationType.PARTIAL_AGONIST,
    binding_affinity=6.0,  # Ki in nM
    efficacy=80.0,  # 80% efficacy
    signaling_bias=SignalingBias.G_PROTEIN_BIASED,
    notes="Primary psychedelic target"
)

# Get primary targets (Ki < 100 nM)
primary_targets = profile.get_primary_targets()

# Export to JSON
json_data = profile.export_to_json()
```

### 4. Get Specific Receptor Info:
```python
# Get 5-HT2A receptor details
receptor_info = db.get_receptor_subtype("serotonin", "5-HT2A")
print(receptor_info['function'])
print(receptor_info['signaling'])
```

---

## Integration Points

### For Compound Analysis:
- Display receptor binding profiles
- Show primary targets vs off-targets
- Indicate modulation type and efficacy
- Highlight signaling bias

### For Analog Generation:
- Optimize for specific receptor selectivity
- Target desired modulation types
- Avoid problematic off-targets (e.g., 5-HT2B)
- Enhance signaling bias

### For PKPD Simulations:
- Model receptor occupancy over time
- Simulate concentration-effect relationships
- Account for partial agonism
- Consider allosteric modulation

### For Docking Studies:
- Validate binding predictions with known Ki values
- Correlate docking scores with experimental data
- Predict modulation type from binding pose
- Estimate efficacy from structural features

---

## Clinical Applications

### Psychedelic Therapy:
- 5-HT2A partial agonism (psychedelic effects)
- 5-HT2C modulation (mood effects)
- 5-HT1A agonism (anxiolytic effects)
- Avoid 5-HT2B agonism (cardiac risk)

### Depression Treatment:
- NMDA-GluN2B antagonism (rapid antidepressant)
- 5-HT1A agonism (anxiolytic, antidepressant)
- mTORC1 activation (neuroplasticity)
- Avoid excessive β-arrestin bias

### Pain Management:
- μ-opioid agonism with G-protein bias (analgesia, reduced side effects)
- NMDA antagonism (neuropathic pain)
- CB1 modulation (pain relief)
- α2A agonism (adjuvant analgesia)

### Anxiety Disorders:
- GABAA-α2/α3 PAM (anxiolytic without sedation)
- 5-HT1A agonism (anxiolytic)
- Avoid GABAA-α1 (sedation)

---

## Testing Results

✅ **Module Tested:** All functions working  
✅ **Total Receptors:** 70 subtypes across 14 families  
✅ **Modulation Types:** 8 types implemented  
✅ **Signaling Bias:** 4 categories  
✅ **Profile Creation:** Working  
✅ **JSON Export:** Working  
✅ **Search Function:** Working  

---

## Future Enhancements

### Potential Additions:
1. **Quantitative SAR** - Structure-activity relationships
2. **Receptor Occupancy Modeling** - Time-dependent occupancy
3. **Polypharmacology Optimization** - Multi-target optimization
4. **Adverse Effect Prediction** - Based on off-target profile
5. **Drug-Drug Interactions** - Receptor competition modeling
6. **Functional Selectivity** - Pathway-specific efficacy
7. **Species Differences** - Human vs rodent receptor profiles
8. **Receptor Polymorphisms** - Genetic variants

---

## Documentation

### Each Receptor Includes:
- Full name
- Receptor family (GPCR, ion channel, etc.)
- G-protein coupling or ion selectivity
- Anatomical location
- Physiological function
- Available signaling pathways
- Clinical relevance

### Each Profile Includes:
- Compound name
- Receptor family and subtype
- Modulation type (agonist, antagonist, PAM, etc.)
- Binding affinity (Ki in nM)
- Efficacy percentage
- Signaling bias (G-protein vs β-arrestin)
- Notes and clinical context

---

## Status

✅ **Implementation Complete**  
✅ **Testing Complete**  
✅ **Ready for Integration**  
✅ **Ready to Commit to GitHub**

---

*Module: `src/receptor_pharmacology.py`*  
*Lines of Code: ~1,200*  
*Date: October 27, 2025*

