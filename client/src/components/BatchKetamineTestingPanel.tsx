import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Badge } from '@/components/ui/badge';
import { Loader2, Zap, Download, BarChart3 } from 'lucide-react';
import { toast } from 'sonner';

interface BatchResult {
  analogId: string;
  compoundName: string;
  receptor: string;
  bindingAffinity: number;
  rmsd: number;
  numPoses: number;
}

export const BatchKetamineTestingPanel: React.FC = () => {
  const [isRunning, setIsRunning] = useState(false);
  const [batchResults, setBatchResults] = useState<BatchResult[]>([]);
  const [selectedReceptors, setSelectedReceptors] = useState<string[]>(['NMDA_GluN2A', 'NMDA_GluN2B']);

  const receptorOptions = [
    { id: 'NMDA_GluN2A', label: 'NMDA GluN2A (Mature neurons)' },
    { id: 'NMDA_GluN2B', label: 'NMDA GluN2B (Developing neurons)' },
    { id: 'DOPAMINE_D2', label: 'Dopamine D2' },
    { id: 'DOPAMINE_D3', label: 'Dopamine D3' },
    { id: 'GABA_A', label: 'GABA-A' },
    { id: 'GABA_B', label: 'GABA-B' },
    { id: '5HT1A', label: 'Serotonin 5HT1A' },
    { id: '5HT2A', label: 'Serotonin 5HT2A' },
  ];

  const handleRunBatchDocking = async () => {
    if (selectedReceptors.length === 0) {
      toast.error('Please select at least one receptor');
      return;
    }

    setIsRunning(true);
    try {
      // Simulate batch docking results
      const mockResults: BatchResult[] = [
        {
          analogId: 'KET-001',
          compoundName: 'Ketamine',
          receptor: 'NMDA_GluN2A',
          bindingAffinity: -7.2,
          rmsd: 1.5,
          numPoses: 9,
        },
        {
          analogId: 'KET-001',
          compoundName: 'Ketamine',
          receptor: 'NMDA_GluN2B',
          bindingAffinity: -7.8,
          rmsd: 1.2,
          numPoses: 9,
        },
        {
          analogId: 'KET-002',
          compoundName: 'Esketamine',
          receptor: 'NMDA_GluN2A',
          bindingAffinity: -7.5,
          rmsd: 1.3,
          numPoses: 9,
        },
        {
          analogId: 'KET-002',
          compoundName: 'Esketamine',
          receptor: 'NMDA_GluN2B',
          bindingAffinity: -8.1,
          rmsd: 1.1,
          numPoses: 9,
        },
      ];

      setBatchResults(mockResults);
      toast.success(`Batch docking completed: ${mockResults.length} results`);
    } catch (error) {
      toast.error('Batch docking failed');
      console.error(error);
    } finally {
      setIsRunning(false);
    }
  };

  const handleExportCSV = () => {
    if (batchResults.length === 0) {
      toast.error('No results to export');
      return;
    }

    const headers = ['Compound', 'Receptor', 'Binding Affinity (kcal/mol)', 'RMSD (Å)', 'Poses'];
    const rows = batchResults.map((r) => [
      r.compoundName,
      r.receptor,
      r.bindingAffinity.toFixed(2),
      r.rmsd.toFixed(2),
      r.numPoses,
    ]);

    const csv = [headers, ...rows].map((row) => row.join(',')).join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `batch-ketamine-results-${new Date().toISOString().split('T')[0]}.csv`;
    a.click();
    window.URL.revokeObjectURL(url);

    toast.success('Results exported to CSV');
  };

  const handleExportJSON = () => {
    if (batchResults.length === 0) {
      toast.error('No results to export');
      return;
    }

    const json = JSON.stringify(batchResults, null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `batch-ketamine-results-${new Date().toISOString().split('T')[0]}.json`;
    a.click();
    window.URL.revokeObjectURL(url);

    toast.success('Results exported to JSON');
  };

  return (
    <Card className="w-full">
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Zap className="w-5 h-5 text-amber-500" />
          Batch Ketamine Testing
        </CardTitle>
        <CardDescription>Run docking simulations for all ketamine analogs against multiple receptors</CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <Tabs defaultValue="setup" className="w-full">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="setup">Setup</TabsTrigger>
            <TabsTrigger value="results">Results ({batchResults.length})</TabsTrigger>
            <TabsTrigger value="export">Export</TabsTrigger>
          </TabsList>

          <TabsContent value="setup" className="space-y-4 mt-4">
            <div>
              <h3 className="font-semibold mb-3">Select Receptors</h3>
              <div className="grid grid-cols-2 gap-3">
                {receptorOptions.map((receptor) => (
                  <label key={receptor.id} className="flex items-center gap-2 p-2 border rounded-lg hover:bg-gray-50 cursor-pointer">
                    <input
                      type="checkbox"
                      checked={selectedReceptors.includes(receptor.id)}
                      onChange={(e) => {
                        if (e.target.checked) {
                          setSelectedReceptors([...selectedReceptors, receptor.id]);
                        } else {
                          setSelectedReceptors(selectedReceptors.filter((r) => r !== receptor.id));
                        }
                      }}
                      className="w-4 h-4"
                    />
                    <span className="text-sm">{receptor.label}</span>
                  </label>
                ))}
              </div>
            </div>

            <Button
              onClick={handleRunBatchDocking}
              disabled={isRunning || selectedReceptors.length === 0}
              className="w-full"
              size="lg"
            >
              {isRunning ? (
                <>
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                  Running Batch Docking...
                </>
              ) : (
                <>
                  <Zap className="w-4 h-4 mr-2" />
                  Start Batch Docking
                </>
              )}
            </Button>
          </TabsContent>

          <TabsContent value="results" className="space-y-4 mt-4">
            {batchResults.length === 0 ? (
              <div className="text-center py-8 text-gray-500">
                <BarChart3 className="w-12 h-12 mx-auto mb-2 opacity-50" />
                <p>No results yet. Run batch docking to see results.</p>
              </div>
            ) : (
              <div className="space-y-3">
                {batchResults.map((result, idx) => (
                  <div key={idx} className="p-4 border rounded-lg hover:bg-gray-50">
                    <div className="flex justify-between items-start mb-2">
                      <div>
                        <p className="font-semibold">{result.compoundName}</p>
                        <p className="text-sm text-gray-600">{result.receptor}</p>
                      </div>
                      <Badge variant={result.bindingAffinity < -7 ? 'default' : 'secondary'}>
                        {result.bindingAffinity.toFixed(2)} kcal/mol
                      </Badge>
                    </div>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>
                        <span className="text-gray-600">RMSD:</span>
                        <span className="ml-2 font-mono">{result.rmsd.toFixed(2)} Å</span>
                      </div>
                      <div>
                        <span className="text-gray-600">Poses:</span>
                        <span className="ml-2 font-mono">{result.numPoses}</span>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            )}
          </TabsContent>

          <TabsContent value="export" className="space-y-4 mt-4">
            {batchResults.length === 0 ? (
              <div className="text-center py-8 text-gray-500">
                <p>No results to export. Run batch docking first.</p>
              </div>
            ) : (
              <div className="space-y-3">
                <Button onClick={handleExportCSV} className="w-full" variant="outline">
                  <Download className="w-4 h-4 mr-2" />
                  Export as CSV
                </Button>
                <Button onClick={handleExportJSON} className="w-full" variant="outline">
                  <Download className="w-4 h-4 mr-2" />
                  Export as JSON
                </Button>
              </div>
            )}
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  );
};

export default BatchKetamineTestingPanel;
