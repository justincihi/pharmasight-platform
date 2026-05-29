import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Badge } from '@/components/ui/badge';
import { Checkbox } from '@/components/ui/checkbox';
import { trpc } from '@/lib/trpc';
import { Loader2, Download, Zap } from 'lucide-react';

interface OptimizedLead {
  smiles: string;
  name: string;
  rank: number;
  predicted_potency: number;
  predicted_selectivity: number;
  predicted_admet_score: number;
  metabolic_stability: number;
  overall_score: number;
  rationale: string;
}

interface PipelineResult {
  parentADMET?: any;
  metabolites?: any;
  sar?: any;
  optimizedLeads?: {
    top_leads: OptimizedLead[];
    optimization_metrics: Record<string, any>;
  };
}

export const LeadOptimizationPanel: React.FC = () => {
  const [parentSmiles, setParentSmiles] = useState('');
  const [targetSmiles, setTargetSmiles] = useState('');
  const [numAnalogs, setNumAnalogs] = useState(20);
  const [topN, setTopN] = useState(10);
  const [includeMetabolites, setIncludeMetabolites] = useState(true);
  const [includeSAR, setIncludeSAR] = useState(true);
  const [activeTab, setActiveTab] = useState('config');
  const [pipelineResults, setPipelineResults] = useState<PipelineResult | null>(null);

  const runPipelineMutation = trpc.leadOptimization.runCompletePipeline.useMutation({
    onSuccess: (data) => {
      setPipelineResults(data.data);
      setActiveTab('results');
    },
    onError: (error) => {
      console.error('Pipeline error:', error);
      alert(`Pipeline failed: ${error.message}`);
    },
  });

  const handleRunPipeline = async () => {
    if (!parentSmiles.trim()) {
      alert('Please enter a parent SMILES string');
      return;
    }

    runPipelineMutation.mutate({
      parentSmiles,
      targetSmiles,
      numAnalogs,
      topN,
      includeMetabolites,
      includeSAR,
    });
  };

  const handleExportResults = () => {
    if (!pipelineResults) return;

    const exportData = {
      timestamp: new Date().toISOString(),
      parentSmiles,
      targetSmiles,
      results: pipelineResults,
    };

    const blob = new Blob([JSON.stringify(exportData, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `lead-optimization-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="w-full space-y-4">
      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-4">
          <TabsTrigger value="config">Configuration</TabsTrigger>
          <TabsTrigger value="results" disabled={!pipelineResults}>
            Results
          </TabsTrigger>
          <TabsTrigger value="leads" disabled={!pipelineResults?.optimizedLeads}>
            Optimized Leads
          </TabsTrigger>
          <TabsTrigger value="analysis" disabled={!pipelineResults}>
            Analysis
          </TabsTrigger>
        </TabsList>

        {/* Configuration Tab */}
        <TabsContent value="config" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Lead Optimization Pipeline</CardTitle>
              <CardDescription>
                Configure and run multi-objective lead optimization with metabolite prediction, ADMET analysis, and SAR
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              {/* Parent SMILES Input */}
              <div className="space-y-2">
                <label className="text-sm font-medium">Parent Compound SMILES</label>
                <Input
                  placeholder="Enter SMILES string (e.g., CC(C)Cc1ccc(cc1)C(C)C(O)=O)"
                  value={parentSmiles}
                  onChange={(e) => setParentSmiles(e.target.value)}
                  disabled={runPipelineMutation.isPending}
                />
              </div>

              {/* Target SMILES Input */}
              <div className="space-y-2">
                <label className="text-sm font-medium">Target Protein SMILES (Optional)</label>
                <Input
                  placeholder="Enter target SMILES for selectivity calculation"
                  value={targetSmiles}
                  onChange={(e) => setTargetSmiles(e.target.value)}
                  disabled={runPipelineMutation.isPending}
                />
              </div>

              {/* Parameters */}
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <label className="text-sm font-medium">Number of Analogs</label>
                  <Input
                    type="number"
                    min="5"
                    max="100"
                    value={numAnalogs}
                    onChange={(e) => setNumAnalogs(parseInt(e.target.value))}
                    disabled={runPipelineMutation.isPending}
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-sm font-medium">Top N Results</label>
                  <Input
                    type="number"
                    min="1"
                    max="50"
                    value={topN}
                    onChange={(e) => setTopN(parseInt(e.target.value))}
                    disabled={runPipelineMutation.isPending}
                  />
                </div>
              </div>

              {/* Options */}
              <div className="space-y-3">
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="metabolites"
                    checked={includeMetabolites}
                    onCheckedChange={(checked) => setIncludeMetabolites(checked as boolean)}
                    disabled={runPipelineMutation.isPending}
                  />
                  <label htmlFor="metabolites" className="text-sm font-medium cursor-pointer">
                    Include Metabolite Prediction
                  </label>
                </div>
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="sar"
                    checked={includeSAR}
                    onCheckedChange={(checked) => setIncludeSAR(checked as boolean)}
                    disabled={runPipelineMutation.isPending}
                  />
                  <label htmlFor="sar" className="text-sm font-medium cursor-pointer">
                    Include SAR Analysis
                  </label>
                </div>
              </div>

              {/* Run Button */}
              <Button
                onClick={handleRunPipeline}
                disabled={runPipelineMutation.isPending || !parentSmiles.trim()}
                className="w-full"
                size="lg"
              >
                {runPipelineMutation.isPending ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running Pipeline...
                  </>
                ) : (
                  <>
                    <Zap className="mr-2 h-4 w-4" />
                    Run Optimization Pipeline
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Results Tab */}
        <TabsContent value="results" className="space-y-4">
          {pipelineResults && (
            <Card>
              <CardHeader>
                <CardTitle>Pipeline Results</CardTitle>
                <CardDescription>Overview of all analysis steps</CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                {/* Parent ADMET */}
                {pipelineResults.parentADMET && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">Parent Compound ADMET</h3>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>MW: {pipelineResults.parentADMET.molecular_weight?.toFixed(1)}</div>
                      <div>LogP: {pipelineResults.parentADMET.logp?.toFixed(2)}</div>
                      <div>HBD: {pipelineResults.parentADMET.hbd}</div>
                      <div>HBA: {pipelineResults.parentADMET.hba}</div>
                    </div>
                  </div>
                )}

                {/* Metabolites Summary */}
                {pipelineResults.metabolites && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">Metabolite Prediction</h3>
                    <div className="text-sm">
                      Total Metabolites:{' '}
                      <Badge>{pipelineResults.metabolites.total_metabolites || 0}</Badge>
                    </div>
                  </div>
                )}

                {/* SAR Summary */}
                {pipelineResults.sar && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">SAR Analysis</h3>
                    <div className="text-sm">
                      Pharmacophore Features:{' '}
                      <Badge variant="outline">
                        {pipelineResults.sar.pharmacophore_features?.length || 0}
                      </Badge>
                    </div>
                  </div>
                )}

                {/* Export Button */}
                <Button onClick={handleExportResults} variant="outline" className="w-full">
                  <Download className="mr-2 h-4 w-4" />
                  Export Results
                </Button>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        {/* Optimized Leads Tab */}
        <TabsContent value="leads" className="space-y-4">
          {pipelineResults?.optimizedLeads && (
            <div className="space-y-4">
              {/* Metrics Summary */}
              <Card>
                <CardHeader>
                  <CardTitle>Optimization Metrics</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="grid grid-cols-2 gap-4 text-sm">
                    <div>
                      <div className="text-gray-500">Total Analogs</div>
                      <div className="text-lg font-semibold">
                        {pipelineResults.optimizedLeads.optimization_metrics.total_analogs_generated}
                      </div>
                    </div>
                    <div>
                      <div className="text-gray-500">Valid Analogs</div>
                      <div className="text-lg font-semibold">
                        {pipelineResults.optimizedLeads.optimization_metrics.valid_analogs}
                      </div>
                    </div>
                    <div>
                      <div className="text-gray-500">Avg Potency</div>
                      <div className="text-lg font-semibold">
                        {pipelineResults.optimizedLeads.optimization_metrics.avg_potency?.toFixed(1)}%
                      </div>
                    </div>
                    <div>
                      <div className="text-gray-500">Best Score</div>
                      <div className="text-lg font-semibold">
                        {pipelineResults.optimizedLeads.optimization_metrics.best_overall_score?.toFixed(1)}
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>

              {/* Top Leads */}
              <div className="space-y-2">
                <h3 className="font-semibold">Top Optimized Leads</h3>
                {pipelineResults.optimizedLeads.top_leads.map((lead: OptimizedLead, idx: number) => (
                  <Card key={idx}>
                    <CardContent className="pt-4">
                      <div className="space-y-2">
                        <div className="flex items-center justify-between">
                          <div className="font-semibold">{lead.name}</div>
                          <Badge variant="default">Rank #{lead.rank}</Badge>
                        </div>
                        <div className="text-xs font-mono text-gray-600 break-all">{lead.smiles}</div>
                        <div className="grid grid-cols-4 gap-2 text-sm">
                          <div>
                            <div className="text-gray-500">Potency</div>
                            <div className="font-semibold">{lead.predicted_potency.toFixed(1)}%</div>
                          </div>
                          <div>
                            <div className="text-gray-500">Selectivity</div>
                            <div className="font-semibold">{lead.predicted_selectivity.toFixed(1)}%</div>
                          </div>
                          <div>
                            <div className="text-gray-500">ADMET</div>
                            <div className="font-semibold">{lead.predicted_admet_score.toFixed(1)}%</div>
                          </div>
                          <div>
                            <div className="text-gray-500">Overall</div>
                            <div className="font-semibold text-blue-600">{lead.overall_score.toFixed(1)}</div>
                          </div>
                        </div>
                        <div className="text-xs text-gray-600">{lead.rationale}</div>
                      </div>
                    </CardContent>
                  </Card>
                ))}
              </div>
            </div>
          )}
        </TabsContent>

        {/* Analysis Tab */}
        <TabsContent value="analysis" className="space-y-4">
          {pipelineResults && (
            <Card>
              <CardHeader>
                <CardTitle>Detailed Analysis</CardTitle>
              </CardHeader>
              <CardContent>
                <pre className="bg-gray-100 p-4 rounded text-xs overflow-auto max-h-96">
                  {JSON.stringify(pipelineResults, null, 2)}
                </pre>
              </CardContent>
            </Card>
          )}
        </TabsContent>
      </Tabs>
    </div>
  );
};
