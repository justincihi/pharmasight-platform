import { useState } from 'react';
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from '@/components/ui/dialog';
import { Button } from '@/components/ui/button';
import { Progress } from '@/components/ui/progress';
import { Badge } from '@/components/ui/badge';
import { trpc } from '@/lib/trpc';
import { Loader2, CheckCircle2, XCircle, Download } from 'lucide-react';
import { toast } from 'sonner';

interface BatchAnalysisModalProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  selectedAnalogs: Array<{ id: number; compoundName: string; smiles: string }>;
}

interface AnalysisResult {
  analogId: number;
  compoundName: string;
  status: 'pending' | 'running' | 'success' | 'error';
  result?: any;
  error?: string;
}

export function BatchAnalysisModal({ open, onOpenChange, selectedAnalogs }: BatchAnalysisModalProps) {
  const [results, setResults] = useState<AnalysisResult[]>([]);
  const [isRunning, setIsRunning] = useState(false);
  const [currentIndex, setCurrentIndex] = useState(0);
  
  const runAnalysisMutation = trpc.advancedAnalysis.comprehensive.useMutation();

  const startBatchAnalysis = async () => {
    setIsRunning(true);
    const initialResults: AnalysisResult[] = selectedAnalogs.map(analog => ({
      analogId: analog.id,
      compoundName: analog.compoundName,
      status: 'pending' as const,
    }));
    setResults(initialResults);

    for (let i = 0; i < selectedAnalogs.length; i++) {
      const analog = selectedAnalogs[i];
      setCurrentIndex(i);
      
      // Update status to running
      setResults(prev => prev.map((r, idx) => 
        idx === i ? { ...r, status: 'running' as const } : r
      ));

      try {
        const result = await runAnalysisMutation.mutateAsync({ smiles: analog.smiles });
        
        // Update with success
        setResults(prev => prev.map((r, idx) => 
          idx === i ? { ...r, status: 'success' as const, result } : r
        ));
      } catch (error: any) {
        // Update with error
        setResults(prev => prev.map((r, idx) => 
          idx === i ? { ...r, status: 'error' as const, error: error.message } : r
        ));
      }
    }

    setIsRunning(false);
    toast.success(`Batch analysis complete: ${results.filter(r => r.status === 'success').length}/${selectedAnalogs.length} succeeded`);
  };

  const exportResults = () => {
    const successfulResults = results.filter(r => r.status === 'success');
    
    if (successfulResults.length === 0) {
      toast.error('No successful results to export');
      return;
    }

    // Create CSV content
    const headers = [
      'Compound Name',
      'Analog ID',
      'hERG Risk',
      'Hepatotoxicity Risk',
      'Mutagenicity Risk',
      'Carcinogenicity Risk',
      'SA Score',
      'SA Difficulty',
      'Optimization Suggestions Count'
    ];

    const rows = successfulResults.map(r => {
      const tox = r.result?.toxicity_profile || {};
      const sa = r.result?.synthetic_accessibility || {};
      const optCount = r.result?.optimization_suggestions?.length || 0;

      return [
        r.compoundName,
        r.analogId,
        tox.hERG?.risk_level || 'N/A',
        tox.hepatotoxicity?.risk_level || 'N/A',
        tox.mutagenicity?.risk_level || 'N/A',
        tox.carcinogenicity?.risk_level || 'N/A',
        sa.sa_score?.toFixed(2) || 'N/A',
        sa.difficulty || 'N/A',
        optCount
      ];
    });

    const csvContent = [
      headers.join(','),
      ...rows.map(row => row.map(cell => `"${cell}"`).join(','))
    ].join('\n');

    // Download CSV
    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `batch_analysis_${new Date().toISOString().split('T')[0]}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);

    toast.success('Results exported to CSV');
  };

  const progress = results.length > 0 
    ? (results.filter(r => r.status !== 'pending').length / results.length) * 100 
    : 0;

  const successCount = results.filter(r => r.status === 'success').length;
  const errorCount = results.filter(r => r.status === 'error').length;

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-w-3xl max-h-[80vh] overflow-y-auto">
        <DialogHeader>
          <DialogTitle>Batch Advanced Analysis</DialogTitle>
          <DialogDescription>
            Running advanced analysis on {selectedAnalogs.length} selected analog{selectedAnalogs.length > 1 ? 's' : ''}
          </DialogDescription>
        </DialogHeader>

        <div className="space-y-4">
          {/* Progress Bar */}
          {results.length > 0 && (
            <div className="space-y-2">
              <div className="flex items-center justify-between text-sm">
                <span className="text-muted-foreground">
                  Progress: {currentIndex + 1} / {selectedAnalogs.length}
                </span>
                <div className="flex gap-2">
                  <Badge variant="outline" className="bg-green-50 text-green-700 border-green-200">
                    <CheckCircle2 className="h-3 w-3 mr-1" />
                    {successCount} Success
                  </Badge>
                  {errorCount > 0 && (
                    <Badge variant="outline" className="bg-red-50 text-red-700 border-red-200">
                      <XCircle className="h-3 w-3 mr-1" />
                      {errorCount} Failed
                    </Badge>
                  )}
                </div>
              </div>
              <Progress value={progress} className="h-2" />
            </div>
          )}

          {/* Results List */}
          {results.length > 0 ? (
            <div className="border rounded-lg divide-y max-h-[400px] overflow-y-auto">
              {results.map((result, index) => (
                <div key={result.analogId} className="p-3 flex items-center justify-between">
                  <div className="flex items-center gap-3">
                    <span className="text-sm font-medium text-muted-foreground">#{index + 1}</span>
                    <div>
                      <div className="font-medium">{result.compoundName}</div>
                      <div className="text-xs text-muted-foreground">ID: {result.analogId}</div>
                    </div>
                  </div>
                  
                  <div>
                    {result.status === 'pending' && (
                      <Badge variant="outline">Pending</Badge>
                    )}
                    {result.status === 'running' && (
                      <Badge variant="outline" className="bg-blue-50 text-blue-700 border-blue-200">
                        <Loader2 className="h-3 w-3 mr-1 animate-spin" />
                        Running
                      </Badge>
                    )}
                    {result.status === 'success' && (
                      <Badge variant="outline" className="bg-green-50 text-green-700 border-green-200">
                        <CheckCircle2 className="h-3 w-3 mr-1" />
                        Complete
                      </Badge>
                    )}
                    {result.status === 'error' && (
                      <Badge variant="outline" className="bg-red-50 text-red-700 border-red-200">
                        <XCircle className="h-3 w-3 mr-1" />
                        Failed
                      </Badge>
                    )}
                  </div>
                </div>
              ))}
            </div>
          ) : (
            <div className="text-center py-8 text-muted-foreground">
              Click "Start Analysis" to begin batch processing
            </div>
          )}

          {/* Action Buttons */}
          <div className="flex justify-between pt-4">
            <Button
              variant="outline"
              onClick={exportResults}
              disabled={successCount === 0}
            >
              <Download className="h-4 w-4 mr-2" />
              Export Results
            </Button>
            
            <div className="flex gap-2">
              <Button
                variant="outline"
                onClick={() => onOpenChange(false)}
                disabled={isRunning}
              >
                {isRunning ? 'Running...' : 'Close'}
              </Button>
              
              {results.length === 0 && (
                <Button onClick={startBatchAnalysis} disabled={isRunning}>
                  {isRunning ? (
                    <><Loader2 className="h-4 w-4 mr-2 animate-spin" />Running...</>
                  ) : (
                    'Start Analysis'
                  )}
                </Button>
              )}
            </div>
          </div>
        </div>
      </DialogContent>
    </Dialog>
  );
}

export default BatchAnalysisModal;
