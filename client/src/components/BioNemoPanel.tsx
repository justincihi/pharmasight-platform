import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Textarea } from '@/components/ui/textarea';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Loader2, Dna, Zap, AlertCircle, CheckCircle2, Activity } from 'lucide-react';
import { toast } from 'sonner';

interface BioNemoPanelProps {
  ligandSmiles?: string;
}

// Common receptor sequences for quick testing
const PRESET_SEQUENCES = {
  'GABA-A alpha1': 'MRPSGFNSTVSRQPLDKFNWTIDMQSLNFHSQKITVNRLLGDAYLFQNLSQTLQNLSQTLQNLSQTLQNLSQ',
  '5-HT2A': 'MDILCEENTSLSSTTNSLMQLNDDTRLYSNDFNSGEANTSDAFNWTVDSENRTNLSCEGCLSPSYQSVPQELNRY',
  'NMDA GluN1': 'MTEQRQNKDGIQKENIQDNQIAQNKDGIQKENIQDNQIAQNKDGIQKENIQDNQIAQNKDGIQKENIQDNQIA',
  'Dopamine D2': 'MDPLNLSWYDDDLERQNWSRPFNGSDGKADRPHYNYYATLLTLLIAVIVFGNVLVCMAVSREKALQTTTNYLIT',
  'Kappa Opioid': 'MDSPIQIFRGEPGPTCSPGALDDGPGSWLNLSHVDGNQSDPCGLPLGPRDGLGGARPGPESAGQLLEAARAGCP',
};

export function BioNemoPanel({ ligandSmiles }: BioNemoPanelProps) {
  const [sequence, setSequence] = useState('');
  const [task, setTask] = useState<'embedding' | 'structure' | 'function' | 'binding_sites'>('binding_sites');
  const [analysisResult, setAnalysisResult] = useState<any>(null);
  const [bindingResult, setBindingResult] = useState<any>(null);
  const [structureResult, setStructureResult] = useState<any>(null);
  const [activeMode, setActiveMode] = useState<'analyze' | 'binding' | 'structure'>('analyze');

  const analyzeMutation = trpc.bionemo.analyzeProtein.useMutation({
    onSuccess: (data) => {
      setAnalysisResult(data);
      toast.success('Protein analysis complete');
    },
    onError: (e) => toast.error(`Analysis failed: ${e.message}`),
  });

  const bindingMutation = trpc.bionemo.predictBinding.useMutation({
    onSuccess: (data) => {
      setBindingResult(data);
      toast.success('Binding prediction complete');
    },
    onError: (e) => toast.error(`Binding prediction failed: ${e.message}`),
  });

  const structureMutation = trpc.bionemo.predictStructure.useMutation({
    onSuccess: (data) => {
      setStructureResult(data);
      toast.success('Structure prediction complete');
    },
    onError: (e) => toast.error(`Structure prediction failed: ${e.message}`),
  });

  const isLoading = analyzeMutation.isPending || bindingMutation.isPending || structureMutation.isPending;

  const handlePreset = (key: string) => {
    setSequence(PRESET_SEQUENCES[key as keyof typeof PRESET_SEQUENCES]);
  };

  const handleRun = () => {
    if (!sequence.trim()) {
      toast.error('Please enter a protein sequence');
      return;
    }
    if (activeMode === 'analyze') {
      analyzeMutation.mutate({ sequence: sequence.trim(), task });
    } else if (activeMode === 'binding') {
      if (!ligandSmiles) {
        toast.error('No ligand SMILES available. Select a compound first.');
        return;
      }
      bindingMutation.mutate({ proteinSequence: sequence.trim(), ligandSmiles });
    } else {
      structureMutation.mutate({ sequence: sequence.trim() });
    }
  };

  return (
    <div className="space-y-4">
      {/* Header */}
      <div className="flex items-start gap-4">
        <Dna className="h-8 w-8 text-purple-500 mt-1 shrink-0" />
        <div>
          <h3 className="font-semibold text-lg">BioNemo Protein Analysis</h3>
          <p className="text-sm text-muted-foreground">
            NVIDIA BioNemo-powered protein analysis: embeddings, structure prediction, binding site detection, and protein-ligand affinity.
          </p>
        </div>
      </div>

      {/* Mode selector */}
      <div className="flex gap-2 flex-wrap">
        {[
          { id: 'analyze', label: 'Protein Analysis', icon: Activity },
          { id: 'binding', label: 'Binding Affinity', icon: Zap },
          { id: 'structure', label: 'Structure Prediction', icon: Dna },
        ].map(({ id, label, icon: Icon }) => (
          <Button
            key={id}
            variant={activeMode === id ? 'default' : 'outline'}
            size="sm"
            onClick={() => setActiveMode(id as any)}
          >
            <Icon className="h-4 w-4 mr-1" />
            {label}
          </Button>
        ))}
      </div>

      {/* Sequence Input */}
      <Card>
        <CardHeader className="pb-3">
          <CardTitle className="text-sm font-medium">Protein Sequence (FASTA or raw)</CardTitle>
          <CardDescription className="text-xs">
            Paste a protein sequence or choose a preset receptor
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-3">
          {/* Preset buttons */}
          <div className="flex flex-wrap gap-1">
            {Object.keys(PRESET_SEQUENCES).map((key) => (
              <Button
                key={key}
                variant="outline"
                size="sm"
                className="text-xs h-7"
                onClick={() => handlePreset(key)}
              >
                {key}
              </Button>
            ))}
          </div>

          <Textarea
            placeholder="Enter protein sequence (e.g. MTEQRQNKDGI...)"
            value={sequence}
            onChange={(e) => setSequence(e.target.value)}
            rows={4}
            className="font-mono text-xs"
          />

          {activeMode === 'analyze' && (
            <Select value={task} onValueChange={(v) => setTask(v as any)}>
              <SelectTrigger className="w-full">
                <SelectValue placeholder="Select analysis task" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="binding_sites">Binding Site Detection</SelectItem>
                <SelectItem value="function">Function Prediction</SelectItem>
                <SelectItem value="embedding">Protein Embedding</SelectItem>
                <SelectItem value="structure">Secondary Structure</SelectItem>
              </SelectContent>
            </Select>
          )}

          {activeMode === 'binding' && ligandSmiles && (
            <div className="p-2 bg-blue-50 dark:bg-blue-950/20 rounded text-xs">
              <span className="font-medium">Ligand SMILES:</span>{' '}
              <code className="text-blue-600">{ligandSmiles.slice(0, 60)}{ligandSmiles.length > 60 ? '...' : ''}</code>
            </div>
          )}

          <Button
            onClick={handleRun}
            disabled={isLoading || !sequence.trim()}
            className="w-full bg-gradient-to-r from-purple-600 to-blue-600 hover:from-purple-700 hover:to-blue-700"
          >
            {isLoading ? (
              <><Loader2 className="h-4 w-4 mr-2 animate-spin" /> Running BioNemo Analysis...</>
            ) : (
              <><Dna className="h-4 w-4 mr-2" /> Run BioNemo Analysis</>
            )}
          </Button>
        </CardContent>
      </Card>

      {/* Analysis Results */}
      {analysisResult && activeMode === 'analyze' && (
        <Card>
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle className="text-sm font-medium flex items-center gap-2">
                <CheckCircle2 className="h-4 w-4 text-green-500" />
                Protein Analysis Results
              </CardTitle>
              {analysisResult.isDemo && (
                <Badge variant="outline" className="text-xs text-amber-600 border-amber-300">
                  Demo Mode
                </Badge>
              )}
            </div>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Predicted Function */}
            {analysisResult.predicted_function && (
              <div>
                <p className="text-xs font-medium text-muted-foreground mb-1">Predicted Function</p>
                <Badge className="capitalize">{analysisResult.predicted_function.replace(/_/g, ' ')}</Badge>
              </div>
            )}

            {/* Binding Sites */}
            {analysisResult.binding_sites?.length > 0 && (
              <div>
                <p className="text-xs font-medium text-muted-foreground mb-2">Predicted Binding Sites</p>
                <div className="space-y-2">
                  {analysisResult.binding_sites.map((site: any, i: number) => (
                    <div key={i} className="flex items-center justify-between p-2 bg-muted/50 rounded text-xs">
                      <div>
                        <span className="font-mono font-medium">{site.residue}</span>
                        <span className="text-muted-foreground ml-2">Position {site.position}</span>
                      </div>
                      <div className="flex items-center gap-2">
                        <div className="w-20 bg-gray-200 rounded-full h-1.5">
                          <div
                            className="bg-purple-500 h-1.5 rounded-full"
                            style={{ width: `${site.confidence * 100}%` }}
                          />
                        </div>
                        <span className="text-muted-foreground">{(site.confidence * 100).toFixed(0)}%</span>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Secondary Structure */}
            {analysisResult.secondary_structure && (
              <div>
                <p className="text-xs font-medium text-muted-foreground mb-1">Secondary Structure</p>
                <code className="text-xs font-mono bg-muted p-2 rounded block overflow-x-auto whitespace-nowrap">
                  {analysisResult.secondary_structure}
                </code>
                <div className="flex gap-3 mt-1 text-xs text-muted-foreground">
                  <span><span className="text-blue-500 font-bold">H</span> = Helix</span>
                  <span><span className="text-green-500 font-bold">E</span> = Sheet</span>
                  <span><span className="text-gray-500 font-bold">C</span> = Coil</span>
                </div>
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Binding Results */}
      {bindingResult && activeMode === 'binding' && (
        <Card>
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle className="text-sm font-medium flex items-center gap-2">
                <Zap className="h-4 w-4 text-yellow-500" />
                Binding Affinity Prediction
              </CardTitle>
              {bindingResult.isDemo && (
                <Badge variant="outline" className="text-xs text-amber-600 border-amber-300">
                  Demo Mode
                </Badge>
              )}
            </div>
          </CardHeader>
          <CardContent className="space-y-3">
            <div className="grid grid-cols-2 gap-3">
              <div className="p-3 bg-purple-50 dark:bg-purple-950/20 rounded-lg text-center">
                <p className="text-xs text-muted-foreground">Predicted Affinity</p>
                <p className="text-xl font-bold text-purple-600">
                  {bindingResult.predicted_affinity_kcal?.toFixed(2)} kcal/mol
                </p>
              </div>
              <div className="p-3 bg-blue-50 dark:bg-blue-950/20 rounded-lg text-center">
                <p className="text-xs text-muted-foreground">Confidence</p>
                <p className="text-xl font-bold text-blue-600">
                  {(bindingResult.confidence * 100).toFixed(0)}%
                </p>
              </div>
            </div>

            <div className="flex items-center gap-2">
              <span className="text-xs text-muted-foreground">Binding Mode:</span>
              <Badge variant="secondary" className="capitalize">{bindingResult.binding_mode}</Badge>
            </div>

            {bindingResult.key_interactions?.length > 0 && (
              <div>
                <p className="text-xs font-medium text-muted-foreground mb-2">Key Interactions</p>
                <div className="space-y-1">
                  {bindingResult.key_interactions.map((int: any, i: number) => (
                    <div key={i} className="flex justify-between text-xs p-2 bg-muted/50 rounded">
                      <span className="font-mono font-medium">{int.residue}</span>
                      <span className="text-muted-foreground capitalize">{int.type.replace(/_/g, ' ')}</span>
                      <span className="text-muted-foreground">{int.distance?.toFixed(1)} Å</span>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Structure Results */}
      {structureResult && activeMode === 'structure' && (
        <Card>
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle className="text-sm font-medium flex items-center gap-2">
                <Dna className="h-4 w-4 text-green-500" />
                Structure Prediction
              </CardTitle>
              {structureResult.isDemo && (
                <Badge variant="outline" className="text-xs text-amber-600 border-amber-300">
                  Demo Mode
                </Badge>
              )}
            </div>
          </CardHeader>
          <CardContent className="space-y-3">
            <div className="grid grid-cols-2 gap-3">
              <div className="p-3 bg-green-50 dark:bg-green-950/20 rounded-lg text-center">
                <p className="text-xs text-muted-foreground">Sequence Length</p>
                <p className="text-xl font-bold text-green-600">{structureResult.length} aa</p>
              </div>
              <div className="p-3 bg-blue-50 dark:bg-blue-950/20 rounded-lg text-center">
                <p className="text-xs text-muted-foreground">pLDDT Score</p>
                <p className="text-xl font-bold text-blue-600">{structureResult.plddt_score?.toFixed(1)}</p>
              </div>
            </div>

            {structureResult.predicted_secondary_structure && (
              <div>
                <p className="text-xs font-medium text-muted-foreground mb-1">Secondary Structure</p>
                <code className="text-xs font-mono bg-muted p-2 rounded block overflow-x-auto whitespace-nowrap">
                  {structureResult.predicted_secondary_structure}
                </code>
                <div className="flex gap-3 mt-1 text-xs text-muted-foreground">
                  <span><span className="text-blue-500 font-bold">H</span> = α-Helix</span>
                  <span><span className="text-green-500 font-bold">E</span> = β-Sheet</span>
                  <span><span className="text-gray-500 font-bold">C</span> = Coil/Loop</span>
                </div>
              </div>
            )}
          </CardContent>
        </Card>
      )}
    </div>
  );
}
