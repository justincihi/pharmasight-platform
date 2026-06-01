import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import DashboardLayout from '@/components/DashboardLayout';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Loader2, Stethoscope, ExternalLink, TrendingUp, Beaker, Target } from 'lucide-react';
import { toast } from 'sonner';

const PHASE_COLOR: Record<number, string> = {
  4: 'text-green-700 bg-green-50 border-green-200',
  3: 'text-blue-700 bg-blue-50 border-blue-200',
  2: 'text-yellow-700 bg-yellow-50 border-yellow-200',
  1: 'text-orange-700 bg-orange-50 border-orange-200',
  0: 'text-gray-600 bg-gray-50 border-gray-200',
};

const PHASE_LABEL: Record<number, string> = {
  4: 'Approved',
  3: 'Phase III',
  2: 'Phase II',
  1: 'Phase I',
  0: 'Preclinical',
};

export default function ClinicalComparability() {
  const [smiles, setSmiles] = useState('');
  const [compoundName, setCompoundName] = useState('');
  const [therapeuticArea, setTherapeuticArea] = useState('');
  const [result, setResult] = useState<any>(null);

  // Also allow selecting from master list
  const analogsQuery = trpc.analog.list.useQuery({ limit: 100, offset: 0 });
  const [selectedAnalogId, setSelectedAnalogId] = useState<string>('');

  const compareMutation = trpc.clinical.compare.useMutation({
    onSuccess: (data) => {
      setResult(data);
      toast.success('Clinical comparability analysis complete');
    },
    onError: (e) => toast.error(`Analysis failed: ${e.message}`),
  });

  const handleSelectAnalog = (id: string) => {
    setSelectedAnalogId(id);
    const analog = analogsQuery.data?.find((a: any) => String(a.id) === id);
    if (analog) {
      setSmiles(analog.smiles);
      setCompoundName(analog.compoundName);
    }
  };

  const handleCompare = () => {
    if (!smiles.trim()) {
      toast.error('Enter a SMILES string');
      return;
    }
    compareMutation.mutate({
      smiles: smiles.trim(),
      compoundName: compoundName || undefined,
      therapeuticArea: therapeuticArea || undefined,
    });
  };

  return (
    <DashboardLayout>
      <div className="p-6 space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3">
          <Stethoscope className="h-7 w-7 text-teal-500" />
          <div>
            <h1 className="text-2xl font-bold">Clinical Comparability Predictor</h1>
            <p className="text-sm text-muted-foreground">
              Compare your analogs to approved drugs using ChEMBL structural similarity and AI-predicted therapeutic class
            </p>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Input */}
          <div className="space-y-4">
            <Card>
              <CardHeader className="pb-3">
                <CardTitle className="text-sm font-medium">Compound Input</CardTitle>
                <CardDescription className="text-xs">Select from master list or enter SMILES manually</CardDescription>
              </CardHeader>
              <CardContent className="space-y-3">
                {/* Select from master list */}
                <div>
                  <p className="text-xs text-muted-foreground mb-1">From Master List</p>
                  <Select value={selectedAnalogId} onValueChange={handleSelectAnalog}>
                    <SelectTrigger>
                      <SelectValue placeholder="Select an analog..." />
                    </SelectTrigger>
                    <SelectContent>
                      {(analogsQuery.data ?? []).map((a: any) => (
                        <SelectItem key={a.id} value={String(a.id)}>
                          {a.compoundName}
                        </SelectItem>
                      ))}
                    </SelectContent>
                  </Select>
                </div>

                <div className="flex items-center gap-2">
                  <div className="flex-1 h-px bg-border" />
                  <span className="text-xs text-muted-foreground">or</span>
                  <div className="flex-1 h-px bg-border" />
                </div>

                <div>
                  <p className="text-xs text-muted-foreground mb-1">SMILES String</p>
                  <Input
                    placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    className="font-mono text-xs"
                  />
                </div>

                <div className="grid grid-cols-2 gap-2">
                  <div>
                    <p className="text-xs text-muted-foreground mb-1">Compound Name (optional)</p>
                    <Input
                      placeholder="e.g. PharmaSight-001"
                      value={compoundName}
                      onChange={(e) => setCompoundName(e.target.value)}
                      className="text-sm"
                    />
                  </div>
                  <div>
                    <p className="text-xs text-muted-foreground mb-1">Therapeutic Area (optional)</p>
                    <Input
                      placeholder="e.g. CNS, Oncology"
                      value={therapeuticArea}
                      onChange={(e) => setTherapeuticArea(e.target.value)}
                      className="text-sm"
                    />
                  </div>
                </div>

                <Button
                  onClick={handleCompare}
                  disabled={compareMutation.isPending || !smiles.trim()}
                  className="w-full bg-teal-600 hover:bg-teal-700 text-white"
                >
                  {compareMutation.isPending ? (
                    <><Loader2 className="h-4 w-4 mr-2 animate-spin" /> Comparing...</>
                  ) : (
                    <><Stethoscope className="h-4 w-4 mr-2" /> Compare to Approved Drugs</>
                  )}
                </Button>
              </CardContent>
            </Card>

            {/* Info card */}
            <Card className="bg-teal-50/50 dark:bg-teal-950/20 border-teal-200 dark:border-teal-800">
              <CardContent className="pt-4 space-y-2">
                <p className="text-xs font-medium text-teal-700 dark:text-teal-300">How it works</p>
                <ul className="text-xs text-teal-600 dark:text-teal-400 space-y-1">
                  <li>• Searches ChEMBL for structurally similar approved drugs (Tanimoto ≥ 60%)</li>
                  <li>• AI predicts therapeutic class, mechanism of action, and likely indication</li>
                  <li>• Shows development stage and clinical comparators</li>
                  <li>• Links to ChEMBL and literature for each comparator</li>
                </ul>
              </CardContent>
            </Card>
          </div>

          {/* Results */}
          <div className="space-y-4">
            {!result && !compareMutation.isPending && (
              <Card className="h-64 flex items-center justify-center">
                <CardContent className="text-center">
                  <Stethoscope className="h-10 w-10 text-muted-foreground mx-auto mb-2 opacity-30" />
                  <p className="text-muted-foreground text-sm">Enter a compound to find clinical comparators</p>
                </CardContent>
              </Card>
            )}

            {compareMutation.isPending && (
              <Card className="h-64 flex items-center justify-center">
                <CardContent className="text-center">
                  <Loader2 className="h-8 w-8 animate-spin text-teal-500 mx-auto mb-2" />
                  <p className="text-muted-foreground text-sm">Searching ChEMBL and predicting therapeutic class...</p>
                </CardContent>
              </Card>
            )}

            {result && (
              <>
                {/* AI Therapeutic Prediction */}
                {result.therapeuticPrediction && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <Target className="h-4 w-4 text-teal-500" />
                        AI Therapeutic Profile — {result.compoundName}
                      </CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-3">
                      <div className="grid grid-cols-2 gap-3">
                        <div className="p-2 bg-teal-50 dark:bg-teal-950/20 rounded">
                          <p className="text-xs text-muted-foreground">Therapeutic Class</p>
                          <p className="text-sm font-medium">{result.therapeuticPrediction.therapeuticClass}</p>
                        </div>
                        <div className="p-2 bg-teal-50 dark:bg-teal-950/20 rounded">
                          <p className="text-xs text-muted-foreground">Development Stage</p>
                          <p className="text-sm font-medium">{result.therapeuticPrediction.developmentStage}</p>
                        </div>
                      </div>
                      <div>
                        <p className="text-xs text-muted-foreground">Mechanism of Action</p>
                        <p className="text-sm">{result.therapeuticPrediction.mechanismOfAction}</p>
                      </div>
                      <div>
                        <p className="text-xs text-muted-foreground">Predicted Indication</p>
                        <p className="text-sm">{result.therapeuticPrediction.indication}</p>
                      </div>
                      {result.therapeuticPrediction.similarApprovedDrugs?.length > 0 && (
                        <div>
                          <p className="text-xs text-muted-foreground mb-1">AI-Predicted Comparators</p>
                          <div className="flex flex-wrap gap-1">
                            {result.therapeuticPrediction.similarApprovedDrugs.map((d: any, i: number) => (
                              <Badge key={i} variant="outline" className="text-xs">
                                {d.name} ({d.similarity})
                              </Badge>
                            ))}
                          </div>
                        </div>
                      )}
                    </CardContent>
                  </Card>
                )}

                {/* ChEMBL Similar Drugs */}
                {result.similarDrugs?.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <Beaker className="h-4 w-4 text-purple-500" />
                        ChEMBL Structural Comparators ({result.similarDrugs.length})
                      </CardTitle>
                      <CardDescription className="text-xs">
                        Approved drugs with ≥ 60% Tanimoto structural similarity
                      </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-2">
                      {result.similarDrugs.map((drug: any, i: number) => (
                        <div key={i} className="flex items-start justify-between gap-3 p-3 border rounded-lg hover:bg-muted/50">
                          <div className="flex-1 min-w-0">
                            <div className="flex items-center gap-2 flex-wrap">
                              <span className="font-medium text-sm">{drug.name || drug.chemblId}</span>
                              <Badge variant="outline" className={`text-xs h-4 px-1 ${PHASE_COLOR[drug.maxPhase] ?? PHASE_COLOR[0]}`}>
                                {PHASE_LABEL[drug.maxPhase] ?? `Phase ${drug.maxPhase}`}
                              </Badge>
                            </div>
                            <div className="flex items-center gap-3 mt-1 text-xs text-muted-foreground">
                              <span>Similarity: <strong>{(drug.similarity * 100).toFixed(0)}%</strong></span>
                              {drug.mw && <span>MW: {Number(drug.mw).toFixed(0)}</span>}
                              {drug.logp && <span>LogP: {Number(drug.logp).toFixed(1)}</span>}
                              {drug.molecularFormula && <span className="font-mono">{drug.molecularFormula}</span>}
                            </div>
                            {drug.therapeuticClass && drug.therapeuticClass !== 'Unknown' && (
                              <p className="text-xs text-muted-foreground mt-0.5">{drug.therapeuticClass}</p>
                            )}
                          </div>
                          <div className="flex items-center gap-2 shrink-0">
                            {/* Similarity bar */}
                            <div className="w-16 bg-gray-200 dark:bg-gray-700 rounded-full h-1.5">
                              <div
                                className="bg-teal-500 h-1.5 rounded-full"
                                style={{ width: `${Math.min(100, drug.similarity * 100)}%` }}
                              />
                            </div>
                            <a
                              href={`https://www.ebi.ac.uk/chembl/compound_report_card/${drug.chemblId}/`}
                              target="_blank"
                              rel="noopener noreferrer"
                              className="text-muted-foreground hover:text-foreground"
                            >
                              <ExternalLink className="h-3.5 w-3.5" />
                            </a>
                          </div>
                        </div>
                      ))}
                    </CardContent>
                  </Card>
                )}

                {result.similarDrugs?.length === 0 && !result.therapeuticPrediction && (
                  <Card>
                    <CardContent className="pt-4 text-center">
                      <TrendingUp className="h-8 w-8 text-muted-foreground mx-auto mb-2 opacity-30" />
                      <p className="text-sm text-muted-foreground">No structurally similar approved drugs found in ChEMBL at ≥ 60% similarity. This compound may be highly novel.</p>
                    </CardContent>
                  </Card>
                )}
              </>
            )}
          </div>
        </div>
      </div>
    </DashboardLayout>
  );
}
