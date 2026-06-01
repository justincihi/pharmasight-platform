import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import DashboardLayout from '@/components/DashboardLayout';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Loader2, Plus, X, FlaskConical, AlertTriangle, CheckCircle2, AlertCircle, Info } from 'lucide-react';
import { toast } from 'sonner';

const SEVERITY_CONFIG: Record<string, { color: string; icon: typeof AlertTriangle }> = {
  major: { color: 'text-red-600 bg-red-50 border-red-200', icon: AlertTriangle },
  moderate: { color: 'text-yellow-600 bg-yellow-50 border-yellow-200', icon: AlertCircle },
  minor: { color: 'text-blue-600 bg-blue-50 border-blue-200', icon: Info },
  contraindicated: { color: 'text-red-800 bg-red-100 border-red-300', icon: AlertTriangle },
};

const RISK_CONFIG: Record<string, { color: string; label: string }> = {
  low: { color: 'text-green-600 bg-green-50 border-green-200', label: 'Low Risk' },
  moderate: { color: 'text-yellow-600 bg-yellow-50 border-yellow-200', label: 'Moderate Risk' },
  high: { color: 'text-red-600 bg-red-50 border-red-200', label: 'High Risk' },
  contraindicated: { color: 'text-red-800 bg-red-100 border-red-300', label: 'Contraindicated' },
};

const PRESET_COMBOS = [
  { label: 'SSRI + MAOI', compounds: [{ name: 'fluoxetine', dose: '20mg' }, { name: 'phenelzine', dose: '15mg' }] },
  { label: 'Ketamine + Benzo', compounds: [{ name: 'ketamine', dose: '0.5mg/kg' }, { name: 'diazepam', dose: '5mg' }] },
  { label: 'Psilocybin + SSRI', compounds: [{ name: 'psilocybin', dose: '25mg' }, { name: 'sertraline', dose: '100mg' }] },
  { label: 'Antipsychotic Stack', compounds: [{ name: 'risperidone', dose: '2mg' }, { name: 'haloperidol', dose: '5mg' }, { name: 'clozapine', dose: '100mg' }] },
];

interface Compound { name: string; dose: string; route: string; }

export default function CocktailAnalyzer() {
  const [compounds, setCompounds] = useState<Compound[]>([
    { name: '', dose: '', route: 'oral' },
    { name: '', dose: '', route: 'oral' },
  ]);
  const [result, setResult] = useState<any>(null);

  const analyzeMutation = trpc.cocktail.analyze.useMutation({
    onSuccess: (data) => {
      setResult(data);
      toast.success('Interaction analysis complete');
    },
    onError: (e) => toast.error(`Analysis failed: ${e.message}`),
  });

  const addCompound = () => {
    if (compounds.length < 8) setCompounds(c => [...c, { name: '', dose: '', route: 'oral' }]);
  };

  const removeCompound = (i: number) => {
    if (compounds.length > 2) setCompounds(c => c.filter((_, idx) => idx !== i));
  };

  const updateCompound = (i: number, field: keyof Compound, value: string) => {
    setCompounds(c => c.map((comp, idx) => idx === i ? { ...comp, [field]: value } : comp));
  };

  const loadPreset = (preset: typeof PRESET_COMBOS[0]) => {
    setCompounds(preset.compounds.map(c => ({ ...c, route: 'oral' })));
    setResult(null);
  };

  const handleAnalyze = () => {
    const valid = compounds.filter(c => c.name.trim());
    if (valid.length < 2) {
      toast.error('Enter at least 2 compound names');
      return;
    }
    analyzeMutation.mutate({
      compounds: valid.map(c => ({ name: c.name.trim(), dose: c.dose || undefined, route: c.route || undefined })),
    });
  };

  return (
    <DashboardLayout>
      <div className="p-6 space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3">
          <FlaskConical className="h-7 w-7 text-red-500" />
          <div>
            <h1 className="text-2xl font-bold">Psychiatric Cocktail Analyzer</h1>
            <p className="text-sm text-muted-foreground">
              Polypharmacy interaction checker — serotonin syndrome, CYP interactions, QTc prolongation, and clinical recommendations
            </p>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Input panel */}
          <div className="space-y-4">
            {/* Presets */}
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium">Quick Presets</CardTitle>
              </CardHeader>
              <CardContent className="flex flex-wrap gap-2">
                {PRESET_COMBOS.map((p) => (
                  <Button key={p.label} variant="outline" size="sm" className="text-xs" onClick={() => loadPreset(p)}>
                    {p.label}
                  </Button>
                ))}
              </CardContent>
            </Card>

            {/* Compound inputs */}
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium">Compounds ({compounds.length}/8)</CardTitle>
                <CardDescription className="text-xs">Enter drug names, analogs, or research compounds</CardDescription>
              </CardHeader>
              <CardContent className="space-y-3">
                {compounds.map((c, i) => (
                  <div key={i} className="flex items-center gap-2">
                    <div className="flex-1 grid grid-cols-2 gap-2">
                      <Input
                        placeholder={`Compound ${i + 1} name`}
                        value={c.name}
                        onChange={(e) => updateCompound(i, 'name', e.target.value)}
                        className="text-sm"
                      />
                      <Input
                        placeholder="Dose (e.g. 20mg)"
                        value={c.dose}
                        onChange={(e) => updateCompound(i, 'dose', e.target.value)}
                        className="text-sm"
                      />
                    </div>
                    {compounds.length > 2 && (
                      <Button variant="ghost" size="sm" className="h-8 w-8 p-0 shrink-0" onClick={() => removeCompound(i)}>
                        <X className="h-3.5 w-3.5" />
                      </Button>
                    )}
                  </div>
                ))}
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" onClick={addCompound} disabled={compounds.length >= 8} className="flex-1">
                    <Plus className="h-3.5 w-3.5 mr-1" /> Add Compound
                  </Button>
                  <Button
                    onClick={handleAnalyze}
                    disabled={analyzeMutation.isPending}
                    className="flex-1 bg-red-600 hover:bg-red-700 text-white"
                  >
                    {analyzeMutation.isPending ? (
                      <><Loader2 className="h-4 w-4 mr-2 animate-spin" /> Analyzing...</>
                    ) : (
                      <><FlaskConical className="h-4 w-4 mr-2" /> Analyze Interactions</>
                    )}
                  </Button>
                </div>
              </CardContent>
            </Card>

            {/* Disclaimer */}
            <div className="p-3 bg-amber-50 dark:bg-amber-950/20 border border-amber-200 dark:border-amber-800 rounded-lg">
              <p className="text-xs text-amber-700 dark:text-amber-300">
                <strong>Research use only.</strong> This tool provides computational predictions for research purposes. Always consult clinical pharmacology resources and healthcare professionals for patient care decisions.
              </p>
            </div>
          </div>

          {/* Results panel */}
          <div className="space-y-4">
            {!result && !analyzeMutation.isPending && (
              <Card className="h-64 flex items-center justify-center">
                <CardContent className="text-center">
                  <FlaskConical className="h-10 w-10 text-muted-foreground mx-auto mb-2 opacity-30" />
                  <p className="text-muted-foreground text-sm">Enter compounds and click Analyze</p>
                </CardContent>
              </Card>
            )}

            {analyzeMutation.isPending && (
              <Card className="h-64 flex items-center justify-center">
                <CardContent className="text-center">
                  <Loader2 className="h-8 w-8 animate-spin text-red-500 mx-auto mb-2" />
                  <p className="text-muted-foreground text-sm">Analyzing interactions...</p>
                </CardContent>
              </Card>
            )}

            {result && (
              <>
                {/* Overall risk */}
                <Card className={`border-2 ${RISK_CONFIG[result.overallRisk]?.color ?? ''}`}>
                  <CardContent className="pt-4">
                    <div className="flex items-center gap-3">
                      {result.overallRisk === 'low' ? (
                        <CheckCircle2 className="h-8 w-8 text-green-500 shrink-0" />
                      ) : (
                        <AlertTriangle className="h-8 w-8 text-red-500 shrink-0" />
                      )}
                      <div>
                        <p className="font-bold text-lg">{RISK_CONFIG[result.overallRisk]?.label ?? result.overallRisk}</p>
                        <p className="text-sm text-muted-foreground">{result.summary}</p>
                      </div>
                    </div>
                  </CardContent>
                </Card>

                {/* Interactions */}
                {result.interactions?.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium">
                        Drug Interactions ({result.interactions.length})
                      </CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-3">
                      {result.interactions.map((int: any, i: number) => {
                        const cfg = SEVERITY_CONFIG[int.severity] ?? SEVERITY_CONFIG.minor;
                        const Icon = cfg.icon;
                        return (
                          <div key={i} className={`p-3 rounded-lg border ${cfg.color}`}>
                            <div className="flex items-start gap-2">
                              <Icon className="h-4 w-4 shrink-0 mt-0.5" />
                              <div className="flex-1 min-w-0">
                                <div className="flex items-center gap-2 flex-wrap">
                                  <span className="font-medium text-sm">{int.drugA} + {int.drugB}</span>
                                  <Badge variant="outline" className={`text-xs h-4 px-1 capitalize ${cfg.color}`}>
                                    {int.severity}
                                  </Badge>
                                  <Badge variant="outline" className="text-xs h-4 px-1 capitalize">
                                    {int.type}
                                  </Badge>
                                </div>
                                <p className="text-xs mt-1">{int.mechanism}</p>
                                {int.recommendation && (
                                  <p className="text-xs mt-1 font-medium">→ {int.recommendation}</p>
                                )}
                              </div>
                            </div>
                          </div>
                        );
                      })}
                    </CardContent>
                  </Card>
                )}

                {/* Monitoring */}
                {result.monitoring?.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <Info className="h-4 w-4 text-blue-500" />
                        Monitoring Recommendations
                      </CardTitle>
                    </CardHeader>
                    <CardContent>
                      <ul className="space-y-1">
                        {result.monitoring.map((m: string, i: number) => (
                          <li key={i} className="text-xs flex items-start gap-1.5">
                            <span className="text-blue-500 mt-0.5">•</span>
                            {m}
                          </li>
                        ))}
                      </ul>
                    </CardContent>
                  </Card>
                )}

                {/* No interactions */}
                {result.interactions?.length === 0 && (
                  <Card>
                    <CardContent className="pt-4 text-center">
                      <CheckCircle2 className="h-8 w-8 text-green-500 mx-auto mb-2" />
                      <p className="text-sm text-muted-foreground">No significant interactions detected for this combination.</p>
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
