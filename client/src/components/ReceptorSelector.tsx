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
  {
    id: 'nmda',
    name: 'NMDA Receptor',
    icon: '🧠',
    description: 'Ionotropic glutamate receptor - Ion channel gating',
    subtypes: [
      { id: 'glun2a', name: 'GluN2A (Mature neurons)', pdbId: '8XLK', species: 'rat' },
      { id: 'glun2b', name: 'GluN2B (Developing neurons)', pdbId: '9JNN', species: 'rat' },
      { id: 'glun2a-glun2b', name: 'GluN2A/GluN2B (Tri-heteromeric)', pdbId: '8XLK', species: 'rat' },
      { id: 'glun2c', name: 'GluN2C (Cerebellum)', species: 'both' },
      { id: 'glun2d', name: 'GluN2D (Thalamus)', species: 'both' },
    ],
  },
  {
    id: 'gaba-a',
    name: 'GABA-A Receptor',
    icon: '⚖️',
    description: 'Ionotropic GABA receptor - Inhibitory ion channel',
    subtypes: [
      { id: 'alpha1', name: 'α1-containing (Sedation)', pdbId: '6D6T', species: 'both' },
      { id: 'alpha2', name: 'α2-containing (Anxiety)', species: 'both' },
      { id: 'alpha3', name: 'α3-containing (Stress)', species: 'both' },
      { id: 'alpha5', name: 'α5-containing (Memory)', species: 'both' },
    ],
  },
  {
    id: 'serotonin',
    name: 'Serotonin Receptor',
    icon: '😊',
    description: 'G-protein coupled receptor - Mood & perception',
    subtypes: [
      { id: '5ht2a', name: '5-HT2A (Psychedelics)', pdbId: '7E2X', species: 'both' },
      { id: '5ht2c', name: '5-HT2C (Appetite)', species: 'both' },
      { id: '5ht1a', name: '5-HT1A (Anxiety)', species: 'both' },
      { id: '5ht1b', name: '5-HT1B (Migraine)', species: 'both' },
      { id: '5ht7', name: '5-HT7 (Sleep)', species: 'both' },
    ],
  },
  {
    id: 'dopamine',
    name: 'Dopamine Receptor',
    icon: '🎯',
    description: 'G-protein coupled receptor - Reward & motivation',
    subtypes: [
      { id: 'd1', name: 'D1 (Stimulatory)', species: 'both' },
      { id: 'd2', name: 'D2 (Antipsychotic)', pdbId: '6A93', species: 'both' },
      { id: 'd3', name: 'D3 (Addiction)', pdbId: '3PBL', species: 'both' },
      { id: 'd4', name: 'D4 (ADHD)', species: 'both' },
      { id: 'd5', name: 'D5 (Cognition)', species: 'both' },
    ],
  },
  {
    id: 'opioid',
    name: 'Opioid Receptor',
    icon: '💊',
    description: 'G-protein coupled receptor - Pain & reward',
    subtypes: [
      { id: 'mu', name: 'μ (Mu) - Analgesia', species: 'both' },
      { id: 'delta', name: 'δ (Delta) - Mood', species: 'both' },
      { id: 'kappa', name: 'κ (Kappa) - Stress', species: 'both' },
      { id: 'nociceptin', name: 'Nociceptin (ORL-1)', species: 'both' },
    ],
  },
  {
    id: 'muscarinic',
    name: 'Muscarinic Receptor',
    icon: '🎭',
    description: 'G-protein coupled acetylcholine receptor',
    subtypes: [
      { id: 'm1', name: 'M1 (Cognition)', pdbId: '5CXV', species: 'both' },
      { id: 'm2', name: 'M2 (Motor)', species: 'both' },
      { id: 'm3', name: 'M3 (Secretion)', species: 'both' },
      { id: 'm4', name: 'M4 (Schizophrenia)', species: 'both' },
      { id: 'm5', name: 'M5 (Reward)', species: 'both' },
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
