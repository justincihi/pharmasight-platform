import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Loader2, Beaker, AlertCircle } from 'lucide-react';

interface MetaboliteViewerProps {
  analogId: number;
}

export function MetaboliteViewer({ analogId }: MetaboliteViewerProps) {
  const [isPredicting, setIsPredicting] = useState(false);
  
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

  return (
    <div className="space-y-4">
      {/* Header with Predict Button */}
      <div className="flex items-center justify-between">
        <div>
          <h3 className="text-lg font-semibold">Metabolite Prediction</h3>
          <p className="text-sm text-muted-foreground">
            Phase I & II metabolites predicted using RDKit CYP450 rules
          </p>
        </div>
        <Button
          onClick={handlePredict}
          disabled={isPredicting}
          variant="default"
        >
          {isPredicting ? (
            <>
              <Loader2 className="mr-2 h-4 w-4 animate-spin" />
              Predicting...
            </>
          ) : (
            <>
              <Beaker className="mr-2 h-4 w-4" />
              {metabolites && metabolites.length > 0 ? 'Re-predict' : 'Predict Metabolites'}
            </>
          )}
        </Button>
      </div>

      {/* Metabolites List */}
      {metabolites && metabolites.length > 0 ? (
        <div className="grid gap-4">
          {metabolites.map((met, idx) => (
            <Card key={idx}>
              <CardHeader>
                <div className="flex items-start justify-between">
                  <div>
                    <CardTitle className="text-base">
                      Metabolite {idx + 1}
                    </CardTitle>
                    <CardDescription className="font-mono text-xs mt-1">
                      {met.smiles}
                    </CardDescription>
                  </div>
                  <Badge variant={met.phase === 'Phase I' ? 'default' : 'secondary'}>
                    {met.phase}
                  </Badge>
                </div>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-2 gap-4 text-sm">
                  <div>
                    <span className="text-muted-foreground">Transformation:</span>
                    <p className="font-medium">{met.transformation}</p>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Enzyme:</span>
                    <p className="font-medium">{met.enzyme}</p>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Probability:</span>
                    <p className="font-medium">{(parseFloat(met.probability) * 100).toFixed(1)}%</p>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Molecular Weight:</span>
                    <p className="font-medium">{parseFloat(met.molecularWeight || '0').toFixed(2)} g/mol</p>
                  </div>
                  <div>
                    <span className="text-muted-foreground">LogP:</span>
                    <p className="font-medium">{parseFloat(met.logP || '0').toFixed(2)}</p>
                  </div>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
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
