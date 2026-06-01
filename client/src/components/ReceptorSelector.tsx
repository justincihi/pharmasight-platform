import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Upload, ChevronRight, Zap } from 'lucide-react';

interface ReceptorSelectorProps {
  onReceptorSelect?: (receptor: SelectedReceptor) => void;
  onPdbUpload?: (file: File) => void;
  allowPdbUpload?: boolean;
}

export interface SelectedReceptor {
  family: string;
  subtype: string;
  species: 'human' | 'mouse' | 'rat' | 'both';
  modulationMode: string;
  ligandType: string;
  pdbId?: string;
}

const RECEPTOR_FAMILIES = [
  // ── Glutamate ──────────────────────────────────────────────────────────────
  {
    id: 'nmda',
    name: 'NMDA Receptor',
    icon: '🧠',
    description: 'Ionotropic glutamate receptor — ketamine, memantine, PCP target',
    subtypes: [
      { id: 'glun1', name: 'GluN1 (Obligate subunit)', pdbId: '5UN1', species: 'both' },
      { id: 'glun2a', name: 'GluN2A (Mature neurons)', pdbId: '5TPB', species: 'rat' },
      { id: 'glun2b', name: 'GluN2B (Extrasynaptic / depression)', pdbId: '5EWM', species: 'rat' },
      { id: 'glun2a-glun2b', name: 'GluN2A/GluN2B (Tri-heteromeric)', pdbId: '8XLK', species: 'rat' },
      { id: 'glun2c', name: 'GluN2C (Cerebellum)', species: 'both' },
      { id: 'glun2d', name: 'GluN2D (Thalamus)', species: 'both' },
    ],
  },
  {
    id: 'ampa',
    name: 'AMPA Receptor',
    icon: '⚡',
    description: 'Fast excitatory ion channel — ampakines, cognition',
    subtypes: [
      { id: 'glua1', name: 'GluA1 (Hippocampus / plasticity)', pdbId: '5L1F', species: 'both' },
      { id: 'glua2', name: 'GluA2 (Ca²⁺ impermeability)', pdbId: '3KG2', species: 'both' },
      { id: 'glua3', name: 'GluA3', species: 'both' },
      { id: 'glua4', name: 'GluA4 (Cerebellum)', species: 'both' },
    ],
  },
  {
    id: 'mglu',
    name: 'Metabotropic Glutamate (mGluR)',
    icon: '🔗',
    description: 'GPCR glutamate modulators — anxiety, schizophrenia',
    subtypes: [
      { id: 'mglu2', name: 'mGluR2 (Presynaptic / anxiety)', pdbId: '7MTQ', species: 'both' },
      { id: 'mglu3', name: 'mGluR3 (Glia / schizophrenia)', species: 'both' },
      { id: 'mglu5', name: 'mGluR5 (Synaptic plasticity)', pdbId: '6FFH', species: 'both' },
      { id: 'mglu1', name: 'mGluR1 (Cerebellum)', species: 'both' },
    ],
  },
  // ── GABA ───────────────────────────────────────────────────────────────────
  {
    id: 'gaba-a',
    name: 'GABA-A Receptor',
    icon: '⚖️',
    description: 'Ionotropic GABA receptor — benzodiazepines, barbiturates, neurosteroids',
    subtypes: [
      { id: 'alpha1', name: 'α1 (Sedation / hypnosis)', pdbId: '6X3Z', species: 'both' },
      { id: 'alpha2', name: 'α2 (Anxiolytic)', species: 'both' },
      { id: 'alpha3', name: 'α3 (Anxiolytic / anticonvulsant)', species: 'both' },
      { id: 'alpha4', name: 'α4 (Tonic inhibition)', species: 'both' },
      { id: 'alpha5', name: 'α5 (Memory / cognition)', species: 'both' },
      { id: 'alpha6', name: 'α6 (Cerebellum)', species: 'both' },
      { id: 'delta', name: 'δ (Extrasynaptic / tonic)', species: 'both' },
    ],
  },
  {
    id: 'gaba-b',
    name: 'GABA-B Receptor',
    icon: '🔒',
    description: 'Gi-coupled GPCR — baclofen, GHB, spasticity',
    subtypes: [
      { id: 'gb1a', name: 'GB1a (Presynaptic)', pdbId: '6UO8', species: 'both' },
      { id: 'gb1b', name: 'GB1b (Postsynaptic)', species: 'both' },
      { id: 'gb2', name: 'GB2 (Obligate partner)', species: 'both' },
    ],
  },
  // ── Serotonin ──────────────────────────────────────────────────────────────
  {
    id: 'serotonin',
    name: 'Serotonin Receptor',
    icon: '😊',
    description: 'GPCR & ion channel — mood, psychedelics, appetite',
    subtypes: [
      { id: '5ht2a', name: '5-HT2A (Psychedelics / antipsychotics)', pdbId: '6A93', species: 'both' },
      { id: '5ht2b', name: '5-HT2B (Cardiac safety screen)', pdbId: '4IB4', species: 'both' },
      { id: '5ht2c', name: '5-HT2C (Appetite / impulsivity)', pdbId: '6BQG', species: 'both' },
      { id: '5ht1a', name: '5-HT1A (Anxiety / autoreceptor)', pdbId: '7E2X', species: 'both' },
      { id: '5ht1b', name: '5-HT1B (Migraine / presynaptic)', pdbId: '4IAQ', species: 'both' },
      { id: '5ht3', name: '5-HT3 (Ion channel / nausea)', pdbId: '6HIN', species: 'both' },
      { id: '5ht6', name: '5-HT6 (Cognitive enhancement)', species: 'both' },
      { id: '5ht7', name: '5-HT7 (Circadian / mood)', species: 'both' },
    ],
  },
  // ── Dopamine ───────────────────────────────────────────────────────────────
  {
    id: 'dopamine',
    name: 'Dopamine Receptor',
    icon: '🎯',
    description: 'Gs/Gi-coupled GPCR — reward, motor control, antipsychotics',
    subtypes: [
      { id: 'd1', name: 'D1 (Gs / reward / motor)', pdbId: '7CKZ', species: 'both' },
      { id: 'd2', name: 'D2 (Gi / antipsychotic target)', pdbId: '6CM4', species: 'both' },
      { id: 'd3', name: 'D3 (Limbic / addiction)', pdbId: '3PBL', species: 'both' },
      { id: 'd4', name: 'D4 (PFC / ADHD)', species: 'both' },
      { id: 'd5', name: 'D5 (Hippocampus / Gs)', species: 'both' },
    ],
  },
  // ── Opioid ─────────────────────────────────────────────────────────────────
  {
    id: 'opioid',
    name: 'Opioid Receptor',
    icon: '💊',
    description: 'Gi-coupled GPCR — analgesia, addiction, mood',
    subtypes: [
      { id: 'mu', name: 'μ (MOR) — Analgesia / euphoria', pdbId: '8EF5', species: 'both' },
      { id: 'delta', name: 'δ (DOR) — Mood / neuroprotection', pdbId: '4N6H', species: 'both' },
      { id: 'kappa', name: 'κ (KOR) — Dysphoria / dissociation', pdbId: '6B73', species: 'both' },
      { id: 'nop', name: 'NOP/ORL-1 — Stress / pain modulation', pdbId: '4EA3', species: 'both' },
    ],
  },
  // ── Adrenergic ─────────────────────────────────────────────────────────────
  {
    id: 'adrenergic',
    name: 'Adrenergic Receptor',
    icon: '⚡',
    description: 'Gq/Gi/Gs-coupled GPCR — PTSD, ADHD, cardiovascular',
    subtypes: [
      { id: 'a1a', name: 'α1A (Vasoconstriction / PTSD)', species: 'both' },
      { id: 'a2a', name: 'α2A (Presynaptic / ADHD / sedation)', pdbId: '6KUX', species: 'both' },
      { id: 'a2b', name: 'α2B (Vascular)', species: 'both' },
      { id: 'a2c', name: 'α2C (Striatum / cognition)', species: 'both' },
      { id: 'b1', name: 'β1 (Cardiac rate)', pdbId: '2VT4', species: 'both' },
      { id: 'b2', name: 'β2 (Bronchodilation)', pdbId: '2RH1', species: 'both' },
    ],
  },
  // ── Cannabinoid ────────────────────────────────────────────────────────────
  {
    id: 'cannabinoid',
    name: 'Cannabinoid Receptor',
    icon: '🌿',
    description: 'Gi-coupled GPCR — THC, pain, nausea, neuroinflammation',
    subtypes: [
      { id: 'cb1', name: 'CB1 (CNS / psychoactive)', pdbId: '5TGZ', species: 'both' },
      { id: 'cb2', name: 'CB2 (Immune / anti-inflammatory)', pdbId: '5ZTY', species: 'both' },
    ],
  },
  // ── Muscarinic ─────────────────────────────────────────────────────────────
  {
    id: 'muscarinic',
    name: 'Muscarinic Receptor',
    icon: '🎭',
    description: 'Gq/Gi-coupled acetylcholine GPCR — cognition, Alzheimer, schizophrenia',
    subtypes: [
      { id: 'm1', name: 'M1 (Cognition / memory)', pdbId: '5CXV', species: 'both' },
      { id: 'm2', name: 'M2 (Cardiac / presynaptic)', pdbId: '3UON', species: 'both' },
      { id: 'm3', name: 'M3 (Smooth muscle / glands)', pdbId: '4U15', species: 'both' },
      { id: 'm4', name: 'M4 (Striatum / schizophrenia)', species: 'both' },
      { id: 'm5', name: 'M5 (Reward / VTA)', species: 'both' },
    ],
  },
  // ── Nicotinic ──────────────────────────────────────────────────────────────
  {
    id: 'nicotinic',
    name: 'Nicotinic Receptor',
    icon: '🔌',
    description: 'Ligand-gated ion channel — cognition, smoking cessation',
    subtypes: [
      { id: 'a7', name: 'α7 (Cognitive enhancement)', pdbId: '7KOO', species: 'both' },
      { id: 'a4b2', name: 'α4β2 (Smoking cessation)', pdbId: '5KXI', species: 'both' },
      { id: 'a3b4', name: 'α3β4 (Autonomic)', species: 'both' },
    ],
  },
  // ── Orexin ─────────────────────────────────────────────────────────────────
  {
    id: 'orexin',
    name: 'Orexin Receptor',
    icon: '😴',
    description: 'Gq-coupled GPCR — wakefulness, appetite, addiction',
    subtypes: [
      { id: 'ox1r', name: 'OX1R (Addiction / arousal)', pdbId: '6TOD', species: 'both' },
      { id: 'ox2r', name: 'OX2R (Sleep-wake / insomnia)', pdbId: '5WQC', species: 'both' },
    ],
  },
  // ── Sigma ──────────────────────────────────────────────────────────────────
  {
    id: 'sigma',
    name: 'Sigma Receptor',
    icon: '🔮',
    description: 'ER chaperone/receptor — neuroplasticity, neuropathic pain',
    subtypes: [
      { id: 'sigma1', name: 'σ1 (ER-mitochondria / neuroplasticity)', pdbId: '5HK1', species: 'both' },
      { id: 'sigma2', name: 'σ2 / TMEM97 (Cancer / Alzheimer)', pdbId: '7M95', species: 'both' },
    ],
  },
  // ── TAAR ───────────────────────────────────────────────────────────────────
  {
    id: 'taar',
    name: 'Trace Amine-Associated Receptor',
    icon: '🧬',
    description: 'Gs-coupled GPCR — novel antipsychotic target, monoamine modulation',
    subtypes: [
      { id: 'taar1', name: 'TAAR1 (Schizophrenia / depression)', species: 'both' },
      { id: 'taar5', name: 'TAAR5 (Social behavior)', species: 'both' },
    ],
  },
  // ── Histamine ──────────────────────────────────────────────────────────────
  {
    id: 'histamine',
    name: 'Histamine Receptor',
    icon: '🤧',
    description: 'Gq/Gs/Gi-coupled GPCR — allergy, sedation, cognition',
    subtypes: [
      { id: 'h1', name: 'H1 (Allergy / sedation)', pdbId: '3RZE', species: 'both' },
      { id: 'h2', name: 'H2 (Gastric acid)', species: 'both' },
      { id: 'h3', name: 'H3 (Presynaptic / cognition)', species: 'both' },
      { id: 'h4', name: 'H4 (Immune / inflammation)', species: 'both' },
    ],
  },
  // ── Transporters ───────────────────────────────────────────────────────────
  {
    id: 'transporter',
    name: 'Monoamine Transporter',
    icon: '🚌',
    description: 'Sodium-dependent reuptake transporters — SSRI, SNRI, stimulant targets',
    subtypes: [
      { id: 'sert', name: 'SERT (Serotonin — SSRI target)', pdbId: '6DZZ', species: 'both' },
      { id: 'dat', name: 'DAT (Dopamine — stimulant target)', pdbId: '4XP4', species: 'both' },
      { id: 'net', name: 'NET (Norepinephrine — SNRI target)', species: 'both' },
    ],
  },
  // ── Enzymes ────────────────────────────────────────────────────────────────
  {
    id: 'mao',
    name: 'Monoamine Oxidase (MAO)',
    icon: '🔬',
    description: 'Mitochondrial enzyme — MAOI antidepressants, Parkinson',
    subtypes: [
      { id: 'mao-a', name: 'MAO-A (Serotonin / NE metabolism)', pdbId: '2Z5X', species: 'both' },
      { id: 'mao-b', name: 'MAO-B (Dopamine / PEA metabolism)', pdbId: '2V5Z', species: 'both' },
    ],
  },
  // ── Kinases / Neuroplasticity ───────────────────────────────────────────────
  {
    id: 'kinase',
    name: 'Kinase / Neuroplasticity Target',
    icon: '🔑',
    description: 'mTOR, TrkB — downstream of ketamine / BDNF signaling',
    subtypes: [
      { id: 'mtor', name: 'mTOR (Synaptogenesis / neuroplasticity)', pdbId: '4JSV', species: 'both' },
      { id: 'trkb', name: 'TrkB / NTRK2 (BDNF receptor)', pdbId: '4AT3', species: 'both' },
    ],
  },
  // ── Cardiac Safety ─────────────────────────────────────────────────────────
  {
    id: 'cardiac-safety',
    name: 'Cardiac Safety Screen',
    icon: '❤️',
    description: 'hERG channel — QT prolongation risk screening',
    subtypes: [
      { id: 'herg', name: 'hERG / Kv11.1 (QT prolongation)', pdbId: '5VA1', species: 'both' },
    ],
  },
  // ── CYP450 ─────────────────────────────────────────────────────────────────
  {
    id: 'cyp450',
    name: 'CYP450 Metabolic Screen',
    icon: '⚗️',
    description: 'Drug-metabolizing enzymes — DDI and metabolic stability screening',
    subtypes: [
      { id: 'cyp3a4', name: 'CYP3A4 (~50% of drugs)', pdbId: '1TQN', species: 'both' },
      { id: 'cyp2d6', name: 'CYP2D6 (Opioids / antidepressants)', pdbId: '2F9Q', species: 'both' },
      { id: 'cyp2c9', name: 'CYP2C9 (Warfarin / NSAIDs)', pdbId: '1OG5', species: 'both' },
      { id: 'cyp1a2', name: 'CYP1A2 (Caffeine / clozapine)', species: 'both' },
    ],
  },
  // ── Melatonin ──────────────────────────────────────────────────────────────
  {
    id: 'melatonin',
    name: 'Melatonin Receptor',
    icon: '🌙',
    description: 'Gi-coupled GPCR — circadian rhythm, sleep, depression',
    subtypes: [
      { id: 'mt1', name: 'MT1 (Circadian / sleep onset)', species: 'both' },
      { id: 'mt2', name: 'MT2 (Phase-shifting)', species: 'both' },
    ],
  },
  // ── Adenosine ──────────────────────────────────────────────────────────────
  {
    id: 'adenosine',
    name: 'Adenosine Receptor',
    icon: '☕',
    description: 'Gs/Gi-coupled GPCR — caffeine target, neuroinflammation, Parkinson',
    subtypes: [
      { id: 'a1', name: 'A1 (Sedation / neuroprotection)', species: 'both' },
      { id: 'a2a', name: 'A2A (Parkinson / neuroinflammation)', pdbId: '3EML', species: 'both' },
      { id: 'a2b', name: 'A2B (Inflammation)', species: 'both' },
    ],
  },
  // ── Glycine ────────────────────────────────────────────────────────────────
  {
    id: 'glycine',
    name: 'Glycine Receptor',
    icon: '🧊',
    description: 'Ligand-gated Cl⁻ channel — inhibitory neurotransmission, pain',
    subtypes: [
      { id: 'glra1', name: 'GlyRα1 (Spinal / pain)', pdbId: '3JAD', species: 'both' },
      { id: 'glra2', name: 'GlyRα2 (Developing brain)', species: 'both' },
    ],
  },
];

const MODULATION_MODES = [
  { id: 'ion-channel', name: 'Ion Channel Gating', description: 'Direct ion channel opening/closing' },
  { id: 'g-protein', name: 'G-Protein Coupled', description: 'Gs, Gi/o, Gq/11, G12/13 pathways' },
  { id: 'beta-arrestin', name: 'β-Arrestin Biased', description: 'β-arrestin signaling pathway' },
  { id: 'kinase', name: 'Kinase Activation', description: 'Receptor tyrosine kinase signaling' },
];

const LIGAND_TYPES = [
  { id: 'agonist', name: 'Full Agonist', description: '100% efficacy' },
  { id: 'partial-agonist', name: 'Partial Agonist', description: '20-80% efficacy' },
  { id: 'antagonist', name: 'Antagonist', description: '0% efficacy' },
  { id: 'pam', name: 'Positive Allosteric Modulator', description: 'Enhances endogenous activity' },
  { id: 'nam', name: 'Negative Allosteric Modulator', description: 'Reduces endogenous activity' },
];

const SPECIES_OPTIONS = [
  { id: 'human', name: 'Human', emoji: '👤' },
  { id: 'mouse', name: 'Mouse', emoji: '🐭' },
  { id: 'rat', name: 'Rat', emoji: '🐀' },
];

export function ReceptorSelector({ onReceptorSelect, onPdbUpload, allowPdbUpload = true }: ReceptorSelectorProps) {
  const [selectedFamily, setSelectedFamily] = useState<string | null>(null);
  const [selectedSubtype, setSelectedSubtype] = useState<string | null>(null);
  const [selectedSpecies, setSelectedSpecies] = useState<'human' | 'mouse' | 'rat'>('human');
  const [selectedModulation, setSelectedModulation] = useState<string>('ion-channel');
  const [selectedLigand, setSelectedLigand] = useState<string>('agonist');
  const [uploadedFile, setUploadedFile] = useState<File | null>(null);

  const family = RECEPTOR_FAMILIES.find((f) => f.id === selectedFamily);
  const subtype = family?.subtypes.find((s) => s.id === selectedSubtype);

  const handleSelectReceptor = () => {
    if (selectedFamily && selectedSubtype) {
      onReceptorSelect?.({
        family: selectedFamily,
        subtype: selectedSubtype,
        species: selectedSpecies,
        modulationMode: selectedModulation,
        ligandType: selectedLigand,
        pdbId: subtype?.pdbId,
      });
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      setUploadedFile(file);
      onPdbUpload?.(file);
    }
  };

  return (
    <div className="w-full space-y-6">
      <Tabs defaultValue="prebuilt" className="w-full">
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="prebuilt">Pre-built Receptors</TabsTrigger>
          {allowPdbUpload && <TabsTrigger value="upload">Upload PDB</TabsTrigger>}
        </TabsList>

        {/* Pre-built Receptors Tab */}
        <TabsContent value="prebuilt" className="space-y-6">
          {/* Step 1: Select Family */}
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <span className="text-2xl">1️⃣</span>
                Select Receptor Family
              </CardTitle>
              <CardDescription>Choose the primary receptor target for your docking</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
                {RECEPTOR_FAMILIES.map((fam) => (
                  <button
                    key={fam.id}
                    onClick={() => {
                      setSelectedFamily(fam.id);
                      setSelectedSubtype(null);
                    }}
                    className={`p-4 rounded-lg border-2 transition-all text-left ${
                      selectedFamily === fam.id
                        ? 'border-blue-500 bg-blue-50'
                        : 'border-gray-200 hover:border-gray-300'
                    }`}
                  >
                    <div className="flex items-start gap-3">
                      <span className="text-3xl">{fam.icon}</span>
                      <div>
                        <h3 className="font-semibold">{fam.name}</h3>
                        <p className="text-sm text-gray-600">{fam.description}</p>
                      </div>
                    </div>
                  </button>
                ))}
              </div>
            </CardContent>
          </Card>

          {/* Step 2: Select Subtype */}
          {selectedFamily && family && (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <span className="text-2xl">2️⃣</span>
                  Select Subtype
                </CardTitle>
                <CardDescription>Choose the specific receptor subtype</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-2">
                  {family.subtypes.map((subtype) => (
                    <button
                      key={subtype.id}
                      onClick={() => setSelectedSubtype(subtype.id)}
                      className={`w-full p-4 rounded-lg border-2 transition-all text-left flex items-center justify-between ${
                        selectedSubtype === subtype.id
                          ? 'border-blue-500 bg-blue-50'
                          : 'border-gray-200 hover:border-gray-300'
                      }`}
                    >
                      <div className="flex-1">
                        <h4 className="font-semibold">{subtype.name}</h4>
                        {subtype.pdbId && (
                          <Badge variant="secondary" className="mt-1">
                            PDB: {subtype.pdbId}
                          </Badge>
                        )}
                      </div>
                      {selectedSubtype === subtype.id && <ChevronRight className="w-5 h-5 text-blue-500" />}
                    </button>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Step 3: Select Species */}
          {selectedSubtype && (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <span className="text-2xl">3️⃣</span>
                  Select Species
                </CardTitle>
                <CardDescription>Choose between human, mouse, or rat structures</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-3 gap-3">
                  {SPECIES_OPTIONS.map((species) => (
                    <button
                      key={species.id}
                      onClick={() => setSelectedSpecies(species.id as any)}
                      className={`p-4 rounded-lg border-2 transition-all text-center ${
                        selectedSpecies === species.id
                          ? 'border-blue-500 bg-blue-50'
                          : 'border-gray-200 hover:border-gray-300'
                      }`}
                    >
                      <div className="text-3xl mb-2">{species.emoji}</div>
                      <div className="font-semibold">{species.name}</div>
                    </button>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Step 4: Select Modulation Mode */}
          {selectedSubtype && (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <span className="text-2xl">4️⃣</span>
                  Signaling Pathway
                </CardTitle>
                <CardDescription>Choose the modulation mode for your ligand</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                  {MODULATION_MODES.map((mode) => (
                    <button
                      key={mode.id}
                      onClick={() => setSelectedModulation(mode.id)}
                      className={`p-4 rounded-lg border-2 transition-all text-left ${
                        selectedModulation === mode.id
                          ? 'border-blue-500 bg-blue-50'
                          : 'border-gray-200 hover:border-gray-300'
                      }`}
                    >
                      <h4 className="font-semibold">{mode.name}</h4>
                      <p className="text-sm text-gray-600">{mode.description}</p>
                    </button>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Step 5: Select Ligand Type */}
          {selectedSubtype && (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <span className="text-2xl">5️⃣</span>
                  Ligand Type
                </CardTitle>
                <CardDescription>Specify the pharmacological action of your compound</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                  {LIGAND_TYPES.map((ligand) => (
                    <button
                      key={ligand.id}
                      onClick={() => setSelectedLigand(ligand.id)}
                      className={`p-4 rounded-lg border-2 transition-all text-left ${
                        selectedLigand === ligand.id
                          ? 'border-blue-500 bg-blue-50'
                          : 'border-gray-200 hover:border-gray-300'
                      }`}
                    >
                      <h4 className="font-semibold">{ligand.name}</h4>
                      <p className="text-sm text-gray-600">{ligand.description}</p>
                    </button>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Summary & Submit */}
          {selectedSubtype && (
            <Card className="bg-gradient-to-r from-blue-50 to-indigo-50 border-blue-200">
              <CardHeader>
                <CardTitle>Configuration Summary</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid grid-cols-2 md:grid-cols-3 gap-4">
                  <div>
                    <p className="text-sm text-gray-600">Family</p>
                    <p className="font-semibold">{family?.name}</p>
                  </div>
                  <div>
                    <p className="text-sm text-gray-600">Subtype</p>
                    <p className="font-semibold">{subtype?.name}</p>
                  </div>
                  <div>
                    <p className="text-sm text-gray-600">Species</p>
                    <p className="font-semibold capitalize">{selectedSpecies}</p>
                  </div>
                  <div>
                    <p className="text-sm text-gray-600">Signaling</p>
                    <p className="font-semibold">
                      {MODULATION_MODES.find((m) => m.id === selectedModulation)?.name}
                    </p>
                  </div>
                  <div>
                    <p className="text-sm text-gray-600">Ligand Type</p>
                    <p className="font-semibold">
                      {LIGAND_TYPES.find((l) => l.id === selectedLigand)?.name}
                    </p>
                  </div>
                  {subtype?.pdbId && (
                    <div>
                      <p className="text-sm text-gray-600">PDB ID</p>
                      <p className="font-semibold">{subtype.pdbId}</p>
                    </div>
                  )}
                </div>

                <Button
                  onClick={handleSelectReceptor}
                  className="w-full bg-blue-600 hover:bg-blue-700 text-white"
                  size="lg"
                >
                  <Zap className="w-4 h-4 mr-2" />
                  Start Docking with This Receptor
                </Button>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        {/* Upload PDB Tab */}
        {allowPdbUpload && (
          <TabsContent value="upload" className="space-y-6">
            <Card>
              <CardHeader>
                <CardTitle>Upload Custom PDB File</CardTitle>
                <CardDescription>Upload your own receptor structure file</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="border-2 border-dashed border-gray-300 rounded-lg p-8 text-center">
                  <Upload className="w-12 h-12 mx-auto text-gray-400 mb-4" />
                  <h3 className="font-semibold mb-2">Upload PDB File</h3>
                  <p className="text-sm text-gray-600 mb-4">Drag and drop or click to select</p>

                  <input
                    type="file"
                    accept=".pdb,.pdbqt"
                    onChange={handleFileUpload}
                    className="hidden"
                    id="pdb-upload"
                  />
                  <label htmlFor="pdb-upload">
                    <Button variant="outline" className="cursor-pointer">
                      Select File
                    </Button>
                  </label>

                  {uploadedFile && (
                    <div className="mt-4 p-3 bg-green-50 border border-green-200 rounded">
                      <p className="text-sm text-green-800">
                        ✓ File selected: <strong>{uploadedFile.name}</strong>
                      </p>
                      <p className="text-xs text-green-700 mt-1">
                        Size: {(uploadedFile.size / 1024).toFixed(2)} KB
                      </p>
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>
          </TabsContent>
        )}
      </Tabs>
    </div>
  );
}
