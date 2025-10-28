"""
Comprehensive Receptor Targeting and Pharmacology Module
Includes all major neurotransmitter receptor types, subtypes, and modulation mechanisms
"""

from typing import Dict, List, Optional, Tuple
from enum import Enum
import json
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ModulationType(Enum):
    """Types of receptor modulation"""
    FULL_AGONIST = "Full Agonist"
    PARTIAL_AGONIST = "Partial Agonist"
    FULL_ANTAGONIST = "Full Antagonist"
    PARTIAL_ANTAGONIST = "Partial Antagonist"
    INVERSE_AGONIST = "Inverse Agonist"
    POSITIVE_ALLOSTERIC_MODULATOR = "Positive Allosteric Modulator (PAM)"
    NEGATIVE_ALLOSTERIC_MODULATOR = "Negative Allosteric Modulator (NAM)"
    SILENT_ALLOSTERIC_MODULATOR = "Silent Allosteric Modulator (SAM)"


class SignalingBias(Enum):
    """G-protein vs β-arrestin signaling bias"""
    G_PROTEIN_BIASED = "G-Protein Biased"
    BETA_ARRESTIN_BIASED = "β-Arrestin Biased"
    BALANCED = "Balanced (No Bias)"
    UNKNOWN = "Unknown"


class ReceptorDatabase:
    """Comprehensive database of neurotransmitter receptors and subtypes"""
    
    RECEPTORS = {
        # Serotonin (5-HT) Receptors
        "serotonin": {
            "name": "Serotonin Receptors",
            "subtypes": {
                "5-HT1A": {
                    "full_name": "5-HT1A Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Raphe nuclei, hippocampus, cortex",
                    "function": "Anxiolytic, antidepressant effects",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "5-HT1B": {
                    "full_name": "5-HT1B Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Basal ganglia, striatum",
                    "function": "Mood regulation, aggression",
                    "signaling": ["Gi/o"],
                },
                "5-HT1D": {
                    "full_name": "5-HT1D Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, blood vessels",
                    "function": "Migraine treatment target",
                    "signaling": ["Gi/o"],
                },
                "5-HT2A": {
                    "full_name": "5-HT2A Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Cortex, platelets",
                    "function": "Psychedelic effects, mood, cognition",
                    "signaling": ["Gq/11", "β-arrestin"],
                },
                "5-HT2B": {
                    "full_name": "5-HT2B Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Heart, GI tract",
                    "function": "Cardiac valve regulation (avoid agonism)",
                    "signaling": ["Gq/11"],
                },
                "5-HT2C": {
                    "full_name": "5-HT2C Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Choroid plexus, cortex",
                    "function": "Appetite, mood, cognition",
                    "signaling": ["Gq/11", "β-arrestin"],
                },
                "5-HT3": {
                    "full_name": "5-HT3 Receptor",
                    "family": "Ligand-gated ion channel",
                    "location": "GI tract, brain stem",
                    "function": "Nausea, vomiting",
                    "signaling": ["Ion channel"],
                },
                "5-HT4": {
                    "full_name": "5-HT4 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "GI tract, hippocampus",
                    "function": "GI motility, cognition",
                    "signaling": ["Gs"],
                },
                "5-HT5A": {
                    "full_name": "5-HT5A Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Cortex, hippocampus",
                    "function": "Unknown (orphan receptor)",
                    "signaling": ["Gi/o"],
                },
                "5-HT6": {
                    "full_name": "5-HT6 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Striatum, cortex, hippocampus",
                    "function": "Cognition, memory",
                    "signaling": ["Gs"],
                },
                "5-HT7": {
                    "full_name": "5-HT7 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Hypothalamus, thalamus",
                    "function": "Circadian rhythm, mood",
                    "signaling": ["Gs"],
                },
            }
        },
        
        # Dopamine Receptors
        "dopamine": {
            "name": "Dopamine Receptors",
            "subtypes": {
                "D1": {
                    "full_name": "Dopamine D1 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Striatum, cortex",
                    "function": "Motor control, reward, cognition",
                    "signaling": ["Gs", "β-arrestin"],
                },
                "D2": {
                    "full_name": "Dopamine D2 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Striatum, pituitary",
                    "function": "Motor control, prolactin regulation",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "D3": {
                    "full_name": "Dopamine D3 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Nucleus accumbens, limbic areas",
                    "function": "Reward, addiction",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "D4": {
                    "full_name": "Dopamine D4 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Frontal cortex, hippocampus",
                    "function": "Cognition, attention",
                    "signaling": ["Gi/o"],
                },
                "D5": {
                    "full_name": "Dopamine D5 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Hippocampus, hypothalamus",
                    "function": "Cognition, renal function",
                    "signaling": ["Gs"],
                },
            }
        },
        
        # GABA Receptors
        "gaba": {
            "name": "GABA Receptors",
            "subtypes": {
                "GABAA-α1": {
                    "full_name": "GABAA Receptor α1 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Widespread CNS",
                    "function": "Sedation, amnesia",
                    "signaling": ["Cl- channel"],
                },
                "GABAA-α2": {
                    "full_name": "GABAA Receptor α2 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus, cortex",
                    "function": "Anxiolytic effects",
                    "signaling": ["Cl- channel"],
                },
                "GABAA-α3": {
                    "full_name": "GABAA Receptor α3 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Cortex, amygdala",
                    "function": "Anxiolytic, muscle relaxation",
                    "signaling": ["Cl- channel"],
                },
                "GABAA-α4": {
                    "full_name": "GABAA Receptor α4 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Thalamus, hippocampus",
                    "function": "Tonic inhibition",
                    "signaling": ["Cl- channel"],
                },
                "GABAA-α5": {
                    "full_name": "GABAA Receptor α5 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus",
                    "function": "Cognition, memory",
                    "signaling": ["Cl- channel"],
                },
                "GABAA-α6": {
                    "full_name": "GABAA Receptor α6 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Cerebellum",
                    "function": "Motor coordination",
                    "signaling": ["Cl- channel"],
                },
                "GABAB": {
                    "full_name": "GABAB Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Widespread CNS",
                    "function": "Muscle relaxation, anxiolytic",
                    "signaling": ["Gi/o"],
                },
            }
        },
        
        # Glutamate Receptors
        "glutamate": {
            "name": "Glutamate Receptors",
            "subtypes": {
                "NMDA-GluN1": {
                    "full_name": "NMDA Receptor GluN1 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Widespread CNS",
                    "function": "Synaptic plasticity, learning, memory",
                    "signaling": ["Ca2+/Na+ channel"],
                },
                "NMDA-GluN2A": {
                    "full_name": "NMDA Receptor GluN2A Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Forebrain",
                    "function": "Synaptic plasticity",
                    "signaling": ["Ca2+/Na+ channel"],
                },
                "NMDA-GluN2B": {
                    "full_name": "NMDA Receptor GluN2B Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Forebrain, spinal cord",
                    "function": "Pain, depression (ketamine target)",
                    "signaling": ["Ca2+/Na+ channel"],
                },
                "NMDA-GluN2C": {
                    "full_name": "NMDA Receptor GluN2C Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Cerebellum, thalamus",
                    "function": "Motor coordination",
                    "signaling": ["Ca2+/Na+ channel"],
                },
                "NMDA-GluN2D": {
                    "full_name": "NMDA Receptor GluN2D Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Midbrain, thalamus",
                    "function": "Sensory processing",
                    "signaling": ["Ca2+/Na+ channel"],
                },
                "AMPA-GluA1": {
                    "full_name": "AMPA Receptor GluA1 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus, cortex",
                    "function": "Fast excitatory transmission",
                    "signaling": ["Na+/K+ channel"],
                },
                "AMPA-GluA2": {
                    "full_name": "AMPA Receptor GluA2 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Widespread CNS",
                    "function": "Fast excitatory transmission",
                    "signaling": ["Na+/K+ channel"],
                },
                "AMPA-GluA3": {
                    "full_name": "AMPA Receptor GluA3 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Cortex, hippocampus",
                    "function": "Synaptic plasticity",
                    "signaling": ["Na+/K+ channel"],
                },
                "AMPA-GluA4": {
                    "full_name": "AMPA Receptor GluA4 Subunit",
                    "family": "Ligand-gated ion channel",
                    "location": "Cerebellum, hippocampus",
                    "function": "Synaptic plasticity",
                    "signaling": ["Na+/K+ channel"],
                },
                "Kainate-GluK1": {
                    "full_name": "Kainate Receptor GluK1",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus, cortex",
                    "function": "Modulation of synaptic transmission",
                    "signaling": ["Na+/K+ channel"],
                },
                "Kainate-GluK2": {
                    "full_name": "Kainate Receptor GluK2",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus, cortex",
                    "function": "Modulation of synaptic transmission",
                    "signaling": ["Na+/K+ channel"],
                },
                "mGluR1": {
                    "full_name": "Metabotropic Glutamate Receptor 1",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Cerebellum, hippocampus",
                    "function": "Synaptic plasticity",
                    "signaling": ["Gq/11"],
                },
                "mGluR5": {
                    "full_name": "Metabotropic Glutamate Receptor 5",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Hippocampus, striatum",
                    "function": "Addiction, pain, anxiety",
                    "signaling": ["Gq/11"],
                },
            }
        },
        
        # Opioid Receptors
        "opioid": {
            "name": "Opioid Receptors",
            "subtypes": {
                "μ (Mu)": {
                    "full_name": "Mu Opioid Receptor (MOR)",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, spinal cord, GI tract",
                    "function": "Analgesia, euphoria, respiratory depression",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "δ (Delta)": {
                    "full_name": "Delta Opioid Receptor (DOR)",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, peripheral tissues",
                    "function": "Analgesia, antidepressant effects",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "κ (Kappa)": {
                    "full_name": "Kappa Opioid Receptor (KOR)",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, spinal cord",
                    "function": "Analgesia, dysphoria, diuresis",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "NOP": {
                    "full_name": "Nociceptin/Orphanin FQ Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, spinal cord",
                    "function": "Pain modulation, anxiety",
                    "signaling": ["Gi/o"],
                },
            }
        },
        
        # Cannabinoid Receptors
        "cannabinoid": {
            "name": "Cannabinoid Receptors",
            "subtypes": {
                "CB1": {
                    "full_name": "Cannabinoid Receptor 1",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, peripheral tissues",
                    "function": "Psychoactive effects, appetite, pain",
                    "signaling": ["Gi/o", "β-arrestin"],
                },
                "CB2": {
                    "full_name": "Cannabinoid Receptor 2",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Immune cells, peripheral tissues",
                    "function": "Immune modulation, inflammation",
                    "signaling": ["Gi/o"],
                },
            }
        },
        
        # Adrenergic Receptors
        "adrenergic": {
            "name": "Adrenergic Receptors",
            "subtypes": {
                "α1A": {
                    "full_name": "Alpha-1A Adrenergic Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Smooth muscle, prostate",
                    "function": "Vasoconstriction, urinary retention",
                    "signaling": ["Gq/11"],
                },
                "α1B": {
                    "full_name": "Alpha-1B Adrenergic Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Smooth muscle, brain",
                    "function": "Vasoconstriction",
                    "signaling": ["Gq/11"],
                },
                "α1D": {
                    "full_name": "Alpha-1D Adrenergic Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Smooth muscle, brain",
                    "function": "Vasoconstriction",
                    "signaling": ["Gq/11"],
                },
                "α2A": {
                    "full_name": "Alpha-2A Adrenergic Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, peripheral tissues",
                    "function": "Sedation, analgesia, hypotension",
                    "signaling": ["Gi/o"],
                },
                "α2B": {
                    "full_name": "Alpha-2B Adrenergic Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Vascular smooth muscle",
                    "function": "Vasoconstriction",
                    "signaling": ["Gi/o"],
                },
                "α2C": {
                    "full_name": "Alpha-2C Adrenergic Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain, peripheral tissues",
                    "function": "Mood regulation",
                    "signaling": ["Gi/o"],
                },
                "β1": {
                    "full_name": "Beta-1 Adrenergic Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Heart",
                    "function": "Increased heart rate and contractility",
                    "signaling": ["Gs", "β-arrestin"],
                },
                "β2": {
                    "full_name": "Beta-2 Adrenergic Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Lungs, smooth muscle",
                    "function": "Bronchodilation, vasodilation",
                    "signaling": ["Gs", "β-arrestin"],
                },
                "β3": {
                    "full_name": "Beta-3 Adrenergic Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Adipose tissue, bladder",
                    "function": "Lipolysis, bladder relaxation",
                    "signaling": ["Gs"],
                },
            }
        },
        
        # Muscarinic Receptors
        "muscarinic": {
            "name": "Muscarinic Acetylcholine Receptors",
            "subtypes": {
                "M1": {
                    "full_name": "Muscarinic M1 Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Cortex, hippocampus, striatum",
                    "function": "Cognition, memory",
                    "signaling": ["Gq/11"],
                },
                "M2": {
                    "full_name": "Muscarinic M2 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Heart, smooth muscle",
                    "function": "Decreased heart rate",
                    "signaling": ["Gi/o"],
                },
                "M3": {
                    "full_name": "Muscarinic M3 Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Smooth muscle, glands",
                    "function": "Glandular secretion, smooth muscle contraction",
                    "signaling": ["Gq/11"],
                },
                "M4": {
                    "full_name": "Muscarinic M4 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Striatum, cortex",
                    "function": "Motor control, cognition",
                    "signaling": ["Gi/o"],
                },
                "M5": {
                    "full_name": "Muscarinic M5 Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Substantia nigra, VTA",
                    "function": "Dopamine release modulation",
                    "signaling": ["Gq/11"],
                },
            }
        },
        
        # Nicotinic Receptors
        "nicotinic": {
            "name": "Nicotinic Acetylcholine Receptors",
            "subtypes": {
                "α4β2": {
                    "full_name": "Nicotinic α4β2 Receptor",
                    "family": "Ligand-gated ion channel",
                    "location": "Brain (widespread)",
                    "function": "Cognition, addiction",
                    "signaling": ["Na+/K+ channel"],
                },
                "α7": {
                    "full_name": "Nicotinic α7 Receptor",
                    "family": "Ligand-gated ion channel",
                    "location": "Hippocampus, cortex",
                    "function": "Cognition, neuroprotection",
                    "signaling": ["Na+/K+ channel"],
                },
                "α3β4": {
                    "full_name": "Nicotinic α3β4 Receptor",
                    "family": "Ligand-gated ion channel",
                    "location": "Autonomic ganglia",
                    "function": "Autonomic transmission",
                    "signaling": ["Na+/K+ channel"],
                },
            }
        },
        
        # Orexin Receptors
        "orexin": {
            "name": "Orexin (Hypocretin) Receptors",
            "subtypes": {
                "OX1": {
                    "full_name": "Orexin Receptor 1",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Hypothalamus, locus coeruleus",
                    "function": "Wakefulness, arousal",
                    "signaling": ["Gq/11"],
                },
                "OX2": {
                    "full_name": "Orexin Receptor 2",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Hypothalamus, widespread",
                    "function": "Sleep-wake regulation",
                    "signaling": ["Gq/11"],
                },
            }
        },
        
        # mTOR (Neuroplasticity)
        "mtor": {
            "name": "mTOR Pathway",
            "subtypes": {
                "mTORC1": {
                    "full_name": "mTOR Complex 1",
                    "family": "Protein kinase complex",
                    "location": "Intracellular (cytoplasm)",
                    "function": "Protein synthesis, neuroplasticity, synaptic growth",
                    "signaling": ["Serine/threonine kinase"],
                },
                "mTORC2": {
                    "full_name": "mTOR Complex 2",
                    "family": "Protein kinase complex",
                    "location": "Intracellular (cytoplasm)",
                    "function": "Cell survival, cytoskeletal organization",
                    "signaling": ["Serine/threonine kinase"],
                },
            }
        },
        
        # Histamine Receptors
        "histamine": {
            "name": "Histamine Receptors",
            "subtypes": {
                "H1": {
                    "full_name": "Histamine H1 Receptor",
                    "family": "GPCR (Gq/11-coupled)",
                    "location": "Brain, smooth muscle",
                    "function": "Wakefulness, allergic response",
                    "signaling": ["Gq/11"],
                },
                "H2": {
                    "full_name": "Histamine H2 Receptor",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Stomach, heart",
                    "function": "Gastric acid secretion",
                    "signaling": ["Gs"],
                },
                "H3": {
                    "full_name": "Histamine H3 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Brain (presynaptic)",
                    "function": "Neurotransmitter release modulation",
                    "signaling": ["Gi/o"],
                },
                "H4": {
                    "full_name": "Histamine H4 Receptor",
                    "family": "GPCR (Gi/o-coupled)",
                    "location": "Immune cells",
                    "function": "Immune modulation",
                    "signaling": ["Gi/o"],
                },
            }
        },
        
        # Sigma Receptors
        "sigma": {
            "name": "Sigma Receptors",
            "subtypes": {
                "σ1": {
                    "full_name": "Sigma-1 Receptor",
                    "family": "Chaperone protein",
                    "location": "ER, mitochondria, plasma membrane",
                    "function": "Neuroprotection, modulation of ion channels",
                    "signaling": ["Chaperone"],
                },
                "σ2": {
                    "full_name": "Sigma-2 Receptor",
                    "family": "Transmembrane protein",
                    "location": "ER, mitochondria",
                    "function": "Cell proliferation, apoptosis",
                    "signaling": ["Unknown"],
                },
            }
        },
        
        # Trace Amine Receptors
        "taar": {
            "name": "Trace Amine-Associated Receptors",
            "subtypes": {
                "TAAR1": {
                    "full_name": "Trace Amine-Associated Receptor 1",
                    "family": "GPCR (Gs-coupled)",
                    "location": "Brain, peripheral tissues",
                    "function": "Modulation of monoaminergic systems",
                    "signaling": ["Gs"],
                },
            }
        },
    }
    
    @classmethod
    def get_all_receptors(cls) -> Dict:
        """Get all receptor families and subtypes"""
        return cls.RECEPTORS
    
    @classmethod
    def get_receptor_family(cls, family_name: str) -> Optional[Dict]:
        """Get specific receptor family"""
        return cls.RECEPTORS.get(family_name.lower())
    
    @classmethod
    def get_receptor_subtype(cls, family_name: str, subtype_name: str) -> Optional[Dict]:
        """Get specific receptor subtype"""
        family = cls.get_receptor_family(family_name)
        if family:
            return family.get('subtypes', {}).get(subtype_name)
        return None
    
    @classmethod
    def list_all_subtypes(cls) -> List[Tuple[str, str, str]]:
        """List all receptor subtypes as (family, subtype, full_name)"""
        subtypes = []
        for family_key, family_data in cls.RECEPTORS.items():
            for subtype_key, subtype_data in family_data.get('subtypes', {}).items():
                subtypes.append((
                    family_key,
                    subtype_key,
                    subtype_data.get('full_name', subtype_key)
                ))
        return subtypes
    
    @classmethod
    def search_receptors(cls, query: str) -> List[Dict]:
        """Search receptors by name or function"""
        results = []
        query_lower = query.lower()
        
        for family_key, family_data in cls.RECEPTORS.items():
            for subtype_key, subtype_data in family_data.get('subtypes', {}).items():
                # Search in name, function, location
                searchable_text = f"{subtype_data.get('full_name', '')} {subtype_data.get('function', '')} {subtype_data.get('location', '')}".lower()
                
                if query_lower in searchable_text:
                    results.append({
                        'family': family_key,
                        'subtype': subtype_key,
                        'data': subtype_data
                    })
        
        return results


class ReceptorProfile:
    """Receptor binding profile for a compound"""
    
    def __init__(self, compound_name: str):
        self.compound_name = compound_name
        self.receptor_interactions = []
    
    def add_interaction(self,
                       receptor_family: str,
                       receptor_subtype: str,
                       modulation_type: ModulationType,
                       binding_affinity: float,  # Ki in nM
                       efficacy: Optional[float] = None,  # % efficacy (0-100)
                       signaling_bias: SignalingBias = SignalingBias.UNKNOWN,
                       notes: str = "") -> None:
        """Add a receptor interaction to the profile"""
        
        interaction = {
            'receptor_family': receptor_family,
            'receptor_subtype': receptor_subtype,
            'modulation_type': modulation_type.value,
            'binding_affinity_ki_nm': binding_affinity,
            'efficacy_percent': efficacy,
            'signaling_bias': signaling_bias.value,
            'notes': notes,
        }
        
        self.receptor_interactions.append(interaction)
    
    def get_profile(self) -> Dict:
        """Get complete receptor profile"""
        return {
            'compound_name': self.compound_name,
            'total_interactions': len(self.receptor_interactions),
            'interactions': self.receptor_interactions,
        }
    
    def get_primary_targets(self, ki_threshold: float = 100.0) -> List[Dict]:
        """Get primary targets (Ki < threshold nM)"""
        return [
            interaction for interaction in self.receptor_interactions
            if interaction['binding_affinity_ki_nm'] < ki_threshold
        ]
    
    def get_off_targets(self, ki_threshold: float = 100.0) -> List[Dict]:
        """Get off-target interactions (Ki >= threshold nM)"""
        return [
            interaction for interaction in self.receptor_interactions
            if interaction['binding_affinity_ki_nm'] >= ki_threshold
        ]
    
    def export_to_json(self) -> str:
        """Export profile to JSON"""
        return json.dumps(self.get_profile(), indent=2)


# Singleton instance
_receptor_database = None

def get_receptor_database() -> ReceptorDatabase:
    """Get or create receptor database instance"""
    global _receptor_database
    if _receptor_database is None:
        _receptor_database = ReceptorDatabase()
    return _receptor_database

