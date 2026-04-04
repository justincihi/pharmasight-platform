/**
 * Comprehensive Receptor Subtype Database
 * Includes species variants, subtypes, and pharmacodynamic modulation modes
 * Data sourced from RCSB PDB, PDSP Database, and neuropharmacology literature
 */

export interface ReceptorSubtype {
  id: string;
  name: string;
  description: string;
  pdbIds: string[];
  species: 'human' | 'rat' | 'mouse' | 'both';
  modulationModes: ModulationMode[];
  ligandTypes: LigandType[];
  notes?: string;
}

export interface ModulationMode {
  id: string;
  name: string;
  description: string;
  examples?: string[];
}

export type LigandType = 'agonist' | 'antagonist' | 'partial-agonist' | 'pam' | 'nam' | 'sam';

export const MODULATION_MODES: Record<string, ModulationMode> = {
  'g-protein': {
    id: 'g-protein',
    name: 'G-Protein Coupled',
    description: 'G-protein signaling pathway activation',
    examples: ['Gs', 'Gi/o', 'Gq/11', 'G12/13'],
  },
  'beta-arrestin': {
    id: 'beta-arrestin',
    name: 'β-Arrestin Biased',
    description: 'β-arrestin signaling pathway activation',
    examples: ['β-arrestin1', 'β-arrestin2'],
  },
  'ion-channel': {
    id: 'ion-channel',
    name: 'Ion Channel Gating',
    description: 'Direct ion channel opening/closing',
    examples: ['NMDA', 'AMPA', 'GABA-A', 'Glycine'],
  },
  'kinase': {
    id: 'kinase',
    name: 'Kinase Activation',
    description: 'Receptor tyrosine kinase signaling',
    examples: ['EGFR', 'PDGFR', 'NGF-TrkA'],
  },
};

export const LIGAND_TYPES: Record<LigandType, string> = {
  'agonist': 'Full Agonist (100% efficacy)',
  'partial-agonist': 'Partial Agonist (20-80% efficacy)',
  'antagonist': 'Antagonist (0% efficacy)',
  'pam': 'Positive Allosteric Modulator',
  'nam': 'Negative Allosteric Modulator',
  'sam': 'Silent Allosteric Modulator',
};

/**
 * NMDA Receptor Subtypes
 * Ionotropic glutamate receptors with distinct pharmacology based on GluN2 subunit composition
 */
export const NMDA_SUBTYPES: ReceptorSubtype[] = [
  {
    id: 'nmda-glun2a',
    name: 'NMDA GluN2A',
    description: 'NMDA receptor with GluN2A subunit (GluN1-GluN2A)',
    pdbIds: ['8XLK'],
    species: 'rat',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Predominantly in mature neurons, associated with learning and memory. Target for cognitive enhancement.',
  },
  {
    id: 'nmda-glun2b',
    name: 'NMDA GluN2B',
    description: 'NMDA receptor with GluN2B subunit (GluN1-GluN2B)',
    pdbIds: ['9JNN'],
    species: 'rat',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Enriched in developing neurons and extrasynaptic locations. Target for neuroprotection and pain.',
  },
  {
    id: 'nmda-glun2a-glun2b',
    name: 'NMDA GluN2A/GluN2B (Tri-heteromeric)',
    description: 'NMDA receptor with both GluN2A and GluN2B subunits',
    pdbIds: ['8XLK'],
    species: 'rat',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Native tri-heteromeric composition found in cortex and hippocampus. Important for synaptic plasticity.',
  },
  {
    id: 'nmda-glun2c',
    name: 'NMDA GluN2C',
    description: 'NMDA receptor with GluN2C subunit',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Primarily in cerebellum and brainstem. Lower conductance than GluN2A/GluN2B.',
  },
  {
    id: 'nmda-glun2d',
    name: 'NMDA GluN2D',
    description: 'NMDA receptor with GluN2D subunit',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Enriched in thalamus and brainstem. Unique pharmacology with reduced sensitivity to Mg2+ block.',
  },
];

/**
 * GABA-A Receptor Subtypes
 * Heteropentameric ion channels with distinct subunit compositions
 */
export const GABA_A_SUBTYPES: ReceptorSubtype[] = [
  {
    id: 'gabaa-alpha1',
    name: 'GABA-A α1-containing',
    description: 'GABA-A receptor with α1 subunit (α1β2γ2)',
    pdbIds: ['6D6T'],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Most abundant GABA-A subtype. Mediates sedation, amnesia, and anticonvulsant effects.',
  },
  {
    id: 'gabaa-alpha2',
    name: 'GABA-A α2-containing',
    description: 'GABA-A receptor with α2 subunit (α2β2γ2)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in anxiety regulation. Selective α2 PAMs show anxiolytic effects without sedation.',
  },
  {
    id: 'gabaa-alpha3',
    name: 'GABA-A α3-containing',
    description: 'GABA-A receptor with α3 subunit (α3β2γ2)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in anxiety and stress responses. Target for anxiolytic development.',
  },
  {
    id: 'gabaa-alpha5',
    name: 'GABA-A α5-containing',
    description: 'GABA-A receptor with α5 subunit (α5β2γ2)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['ion-channel']],
    ligandTypes: ['agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in learning and memory. α5 NAMs show pro-cognitive effects.',
  },
];

/**
 * Serotonin Receptor Subtypes
 * G-protein coupled receptors with diverse pharmacology and tissue distribution
 */
export const SEROTONIN_SUBTYPES: ReceptorSubtype[] = [
  {
    id: '5ht2a',
    name: '5-HT2A',
    description: 'Serotonin 5-HT2A receptor',
    pdbIds: ['7E2X'],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Primary target for psychedelic compounds (LSD, psilocybin). Involved in hallucinogenic effects.',
  },
  {
    id: '5ht2c',
    name: '5-HT2C',
    description: 'Serotonin 5-HT2C receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in appetite regulation and mood. Target for antidepressants and weight loss drugs.',
  },
  {
    id: '5ht1a',
    name: '5-HT1A',
    description: 'Serotonin 5-HT1A receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist'],
    notes: 'Autoreceptor on serotonin neurons. Target for SSRIs and anxiolytics.',
  },
  {
    id: '5ht1b',
    name: '5-HT1B',
    description: 'Serotonin 5-HT1B receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist'],
    notes: 'Presynaptic heteroreceptor. Involved in migraine and mood regulation.',
  },
  {
    id: '5ht1d',
    name: '5-HT1D',
    description: 'Serotonin 5-HT1D receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist'],
    notes: 'Similar to 5-HT1B. Target for migraine treatment (triptans).',
  },
  {
    id: '5ht7',
    name: '5-HT7',
    description: 'Serotonin 5-HT7 receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam'],
    notes: 'Involved in circadian rhythm and mood. Emerging target for depression and sleep disorders.',
  },
];

/**
 * Dopamine Receptor Subtypes
 * G-protein coupled receptors with distinct signaling and tissue distribution
 */
export const DOPAMINE_SUBTYPES: ReceptorSubtype[] = [
  {
    id: 'dopamine-d1',
    name: 'Dopamine D1',
    description: 'Dopamine D1 receptor (D1-like family)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Stimulatory receptor. Target for Parkinsons and addiction.',
  },
  {
    id: 'dopamine-d2',
    name: 'Dopamine D2',
    description: 'Dopamine D2 receptor (D2-like family)',
    pdbIds: ['6A93'],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Primary antipsychotic target. Involved in reward and motor control.',
  },
  {
    id: 'dopamine-d3',
    name: 'Dopamine D3',
    description: 'Dopamine D3 receptor (D2-like family)',
    pdbIds: ['3PBL'],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in reward and motivation. Target for addiction treatment.',
  },
  {
    id: 'dopamine-d4',
    name: 'Dopamine D4',
    description: 'Dopamine D4 receptor (D2-like family)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in attention and impulse control. Emerging ADHD target.',
  },
  {
    id: 'dopamine-d5',
    name: 'Dopamine D5',
    description: 'Dopamine D5 receptor (D1-like family)',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Similar to D1 but with distinct tissue distribution.',
  },
];

/**
 * Opioid Receptor Subtypes
 * G-protein coupled receptors involved in pain and reward
 */
export const OPIOID_SUBTYPES: ReceptorSubtype[] = [
  {
    id: 'opioid-mu',
    name: 'Opioid μ (Mu)',
    description: 'Mu opioid receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam'],
    notes: 'Primary analgesic target. Involved in addiction and reward.',
  },
  {
    id: 'opioid-delta',
    name: 'Opioid δ (Delta)',
    description: 'Delta opioid receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam'],
    notes: 'Involved in pain, mood, and immune function.',
  },
  {
    id: 'opioid-kappa',
    name: 'Opioid κ (Kappa)',
    description: 'Kappa opioid receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam'],
    notes: 'Involved in pain, mood, and stress response. Emerging antidepressant target.',
  },
  {
    id: 'opioid-nociceptin',
    name: 'Opioid Nociceptin (ORL-1)',
    description: 'Nociceptin/orphanin FQ receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam'],
    notes: 'Involved in pain and stress. Emerging target for anxiety and depression.',
  },
];

/**
 * Muscarinic Receptor Subtypes
 * G-protein coupled acetylcholine receptors
 */
export const MUSCARINIC_SUBTYPES: ReceptorSubtype[] = [
  {
    id: 'muscarinic-m1',
    name: 'Muscarinic M1',
    description: 'M1 muscarinic acetylcholine receptor',
    pdbIds: ['5CXV'],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein'], MODULATION_MODES['beta-arrestin']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in cognition and memory. Target for Alzheimers and schizophrenia.',
  },
  {
    id: 'muscarinic-m2',
    name: 'Muscarinic M2',
    description: 'M2 muscarinic acetylcholine receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Autoreceptor on cholinergic neurons. Involved in motor control.',
  },
  {
    id: 'muscarinic-m3',
    name: 'Muscarinic M3',
    description: 'M3 muscarinic acetylcholine receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in secretion and smooth muscle contraction.',
  },
  {
    id: 'muscarinic-m4',
    name: 'Muscarinic M4',
    description: 'M4 muscarinic acetylcholine receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in motor control and cognition. Target for schizophrenia.',
  },
  {
    id: 'muscarinic-m5',
    name: 'Muscarinic M5',
    description: 'M5 muscarinic acetylcholine receptor',
    pdbIds: [],
    species: 'both',
    modulationModes: [MODULATION_MODES['g-protein']],
    ligandTypes: ['agonist', 'partial-agonist', 'antagonist', 'pam', 'nam'],
    notes: 'Involved in reward and motivation. Emerging target for addiction.',
  },
];

/**
 * Get all receptor subtypes for a given family
 */
export function getReceptorSubtypes(family: string): ReceptorSubtype[] {
  switch (family.toLowerCase()) {
    case 'nmda':
    case 'glutamate-nmda':
      return NMDA_SUBTYPES;
    case 'gaba-a':
    case 'gabaa':
      return GABA_A_SUBTYPES;
    case 'serotonin':
    case '5ht':
      return SEROTONIN_SUBTYPES;
    case 'dopamine':
      return DOPAMINE_SUBTYPES;
    case 'opioid':
      return OPIOID_SUBTYPES;
    case 'muscarinic':
      return MUSCARINIC_SUBTYPES;
    default:
      return [];
  }
}

/**
 * Get all available receptor families
 */
export function getAllReceptorFamilies() {
  return [
    { id: 'nmda', name: 'NMDA (Glutamate)', subtypes: NMDA_SUBTYPES },
    { id: 'gaba-a', name: 'GABA-A', subtypes: GABA_A_SUBTYPES },
    { id: 'serotonin', name: 'Serotonin (5-HT)', subtypes: SEROTONIN_SUBTYPES },
    { id: 'dopamine', name: 'Dopamine', subtypes: DOPAMINE_SUBTYPES },
    { id: 'opioid', name: 'Opioid', subtypes: OPIOID_SUBTYPES },
    { id: 'muscarinic', name: 'Muscarinic', subtypes: MUSCARINIC_SUBTYPES },
  ];
}
