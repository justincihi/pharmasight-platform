import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import { Button } from '@/components/ui/button';
import { Card, CardContent } from '@/components/ui/card';
import { Loader2, Beaker, AlertCircle, TreePine, List } from 'lucide-react';
import { MetabolitePathwayTree } from './MetabolitePathwayTree';

interface MetaboliteViewerProps {
  analogId: number;
  parentName?: string;
  parentSmiles?: string;
}

export function MetaboliteViewer({ analogId, parentName = 'Parent Compound', parentSmiles = '' }: MetaboliteViewerProps) {
  const [isPredicting, setIsPredicting] = useState(false);
  const [viewMode, setViewMode] = useState<'tree' | 'cards'>('tree');
  
  const { data: metabolites, isLoading, refetch } = trpc.analog.getMetabolites.useQuery({ analogId });
  const predictMutation = trpc.analog.predictMetabolites.useMutation();

  const handlePredict = async () => {
    setIsPredicting(true);
    try {
      await predictMutation.mutateAsync({ analogId });
      await refetch();
    } catch (error) {
      console.error('Metabolite prediction failed:', error);
      alert('Failed to predict metabolites. Please try again.');
    } finally {
      setIsPredicting(false);
    }
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center p-8">
        <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
      </div>
    );
  }

  // Map stored metabolites to the shape MetabolitePathwayTree expects
  const treeMetabolites = (metabolites ?? []).map((m: any, i: number) => ({
    id: m.id ?? i,
    analogId,
    metaboliteSmiles: m.smiles ?? m.metaboliteSmiles ?? '',
    metaboliteName: m.name ?? m.metaboliteName ?? null,
    reactionType: m.transformation ?? m.reactionType ?? null,
    enzyme: m.enzyme ?? null,
    phase: m.phase === 'Phase I' ? 1 : m.phase === 'Phase II' ? 2 : (m.phase ?? 1),
    molecularWeight: m.molecularWeight != null ? parseFloat(m.molecularWeight) : null,
    logP: m.logP != null ? parseFloat(m.logP) : null,
    predictedActivity: m.predictedActivity ?? null,
    confidenceScore: m.probability != null ? parseFloat(m.probability) : (m.confidenceScore ?? null),
    createdAt: m.createdAt ?? new Date(),
  }));

  return (
    <div className="space-y-4">
      {/* Header with Predict Button + view toggle */}
      <div className="flex items-center justify-between flex-wrap gap-2">
        <div>
          <h3 className="text-lg font-semibold">Metabolite Pathway</h3>
          <p className="text-sm text-muted-foreground">
            Phase I &amp; II metabolites predicted using RDKit SMARTS engine (10 reaction rules)
          </p>
        </div>
        <div className="flex items-center gap-2">
          {metabolites && metabolites.length > 0 && (
            <div className="flex rounded-md border overflow-hidden text-xs">
              <button
                onClick={() => setViewMode('tree')}
                className={`px-2 py-1.5 flex items-center gap-1 transition-colors ${
                  viewMode === 'tree' ? 'bg-primary text-primary-foreground' : 'hover:bg-muted'
                }`}
              >
                <TreePine className="h-3 w-3" />
                Pathway
              </button>
              <button
                onClick={() => setViewMode('cards')}
                className={`px-2 py-1.5 flex items-center gap-1 transition-colors ${
                  viewMode === 'cards' ? 'bg-primary text-primary-foreground' : 'hover:bg-muted'
                }`}
              >
                <List className="h-3 w-3" />
                Cards
              </button>
            </div>
          )}
          <Button onClick={handlePredict} disabled={isPredicting} variant="default" size="sm">
            {isPredicting ? (
              <><Loader2 className="mr-2 h-4 w-4 animate-spin" />Predicting...</>
            ) : (
              <><Beaker className="mr-2 h-4 w-4" />{metabolites && metabolites.length > 0 ? 'Re-predict' : 'Predict Metabolites'}</>
            )}
          </Button>
        </div>
      </div>

      {/* Metabolites display */}
      {metabolites && metabolites.length > 0 ? (
        viewMode === 'tree' ? (
          <MetabolitePathwayTree
            parentName={parentName}
            parentSmiles={parentSmiles}
            metabolites={treeMetabolites}
          />
        ) : (
          <div className="grid gap-4">
            {treeMetabolites.map((met, idx) => (
              <Card key={met.id}>
                <CardContent className="pt-4">
                  <div className="flex items-start justify-between mb-3">
                    <div>
                      <p className="font-semibold text-sm">{met.metaboliteName ?? `Metabolite ${idx + 1}`}</p>
                      <p className="font-mono text-xs text-muted-foreground mt-0.5">{met.metaboliteSmiles}</p>
                    </div>
                    <span className={`text-xs font-bold ${
                      met.phase === 1 ? 'text-blue-500' : 'text-purple-500'
                    }`}>Phase {met.phase}</span>
                  </div>
                  <div className="grid grid-cols-2 gap-3 text-sm">
                    <div><span className="text-muted-foreground text-xs">Transformation:</span><p className="font-medium text-xs">{met.reactionType ?? '—'}</p></div>
                    <div><span className="text-muted-foreground text-xs">Enzyme:</span><p className="font-medium text-xs">{met.enzyme ?? '—'}</p></div>
                    <div><span className="text-muted-foreground text-xs">Confidence:</span><p className="font-medium text-xs">{met.confidenceScore != null ? (met.confidenceScore * 100).toFixed(1) + '%' : '—'}</p></div>
                    <div><span className="text-muted-foreground text-xs">MW:</span><p className="font-medium text-xs">{met.molecularWeight?.toFixed(2) ?? '—'} g/mol</p></div>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )
      ) : (
        <Card>
          <CardContent className="flex flex-col items-center justify-center p-8 text-center">
            <AlertCircle className="h-12 w-12 text-muted-foreground mb-4" />
            <p className="text-muted-foreground">
              No metabolites predicted yet. Click "Predict Metabolites" to generate predictions.
            </p>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
