import React, { useState } from 'react';
import { trpc } from '@/lib/trpc';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Loader2, AlertCircle, CheckCircle2, Download } from 'lucide-react';

interface AnalogResult {
  smiles: string;
  tanimoto: number;
  mw?: number;
  cid?: number;
  name?: string;
  patent_free?: boolean;
  patents?: string[];
  flag?: string;
}

interface PipelineResult {
  success: boolean;
  error?: string;
  data?: {
    parent_smiles?: string;
    parent_cid?: number;
    parent_name?: string;
    hits?: AnalogResult[];
    total_hits?: number;
    analogs?: AnalogResult[];
    total_generated?: number;
    master_list?: AnalogResult[];
    total_candidates?: number;
    patent_free_count?: number;
    canonical_smiles?: string;
    mw?: number;
    num_atoms?: number;
    num_bonds?: number;
  };
}

export function CheminformaticsPipeline() {
  const [activeTab, setActiveTab] = useState('similarity');
  const [inputSmiles, setInputSmiles] = useState('COc1ccc2[nH]cc(CCN(C)C)c2c1'); // 5-MeO-DMT example
  const [threshold, setThreshold] = useState(0.70);
  const [maxHits, setMaxHits] = useState(25);
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState<PipelineResult | null>(null);

  // tRPC queries
  const validateSmilesQuery = trpc.cheminformatics.validateSmiles.useQuery(
    { smiles: inputSmiles },
    { enabled: false }
  );

  const similarityQuery = trpc.cheminformatics.confirmAndFetchSimilars.useQuery(
    { nameOrSmiles: inputSmiles, threshold, maxHits },
    { enabled: false }
  );

  const bricsQuery = trpc.cheminformatics.generateBricsAnalogs.useQuery(
    { smiles: inputSmiles, n: maxHits },
    { enabled: false }
  );

  const fullPipelineQuery = trpc.cheminformatics.fullAnalogPipeline.useQuery(
    { inputSmiles, threshold, maxHits },
    { enabled: false }
  );

  const handleValidateSmiles = async () => {
    setIsLoading(true);
    try {
      const result = await validateSmilesQuery.refetch();
      if (result.data) {
        setResults(result.data as any);
        console.log('SMILES Validated:', result.data.data?.canonical_smiles);
      }
    } catch (error) {
        console.error('Validation Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const handleSimilaritySearch = async () => {
    setIsLoading(true);
    try {
      const result = await similarityQuery.refetch();
      if (result.data) {
        setResults(result.data as any);
        console.log('Similarity Search Complete:', result.data.data?.total_hits);
      }
    } catch (error) {
        console.error('Search Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const handleBricsGeneration = async () => {
    setIsLoading(true);
    try {
      const result = await bricsQuery.refetch();
      if (result.data) {
        setResults(result.data as any);
        console.log('BRICS Generation Complete:', result.data.data?.total_generated);
      }
    } catch (error) {
        console.error('Generation Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const handleFullPipeline = async () => {
    setIsLoading(true);
    try {
      const result = await fullPipelineQuery.refetch();
      if (result.data) {
        setResults(result.data as any);
        console.log('Full Pipeline Complete:', result.data.data?.total_candidates);
      }
    } catch (error) {
        console.error('Pipeline Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const downloadResults = () => {
    if (!results?.data) return;
    const json = JSON.stringify(results.data, null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `cheminformatics-results-${Date.now()}.json`;
    a.click();
  };

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>Cheminformatics Pipeline</CardTitle>
          <CardDescription>
            PubChem similarity screening, patent checking, and analog generation
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="smiles">SMILES String</Label>
            <Input
              id="smiles"
              placeholder="Enter compound SMILES (e.g., COc1ccc2[nH]cc(CCN(C)C)c2c1)"
              value={inputSmiles}
              onChange={(e) => setInputSmiles(e.target.value)}
              disabled={isLoading}
            />
          </div>

          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="threshold">Similarity Threshold</Label>
              <Input
                id="threshold"
                type="number"
                min="0"
                max="1"
                step="0.05"
                value={threshold}
                onChange={(e) => setThreshold(parseFloat(e.target.value))}
                disabled={isLoading}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="maxHits">Max Results</Label>
              <Input
                id="maxHits"
                type="number"
                min="1"
                max="100"
                value={maxHits}
                onChange={(e) => setMaxHits(parseInt(e.target.value))}
                disabled={isLoading}
              />
            </div>
          </div>
        </CardContent>
      </Card>

      <Tabs value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="grid w-full grid-cols-4">
          <TabsTrigger value="similarity">Similarity</TabsTrigger>
          <TabsTrigger value="brics">BRICS</TabsTrigger>
          <TabsTrigger value="validate">Validate</TabsTrigger>
          <TabsTrigger value="full">Full Pipeline</TabsTrigger>
        </TabsList>

        <TabsContent value="similarity" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>PubChem Similarity Search</CardTitle>
              <CardDescription>
                Find structurally similar compounds using Tanimoto similarity
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <Button
                onClick={handleSimilaritySearch}
                disabled={isLoading || !inputSmiles}
                className="w-full"
              >
                {isLoading && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                Search Similar Compounds
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="brics" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>BRICS Analog Generation</CardTitle>
              <CardDescription>
                Generate novel analogs via BRICS fragmentation and reassembly
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <Button
                onClick={handleBricsGeneration}
                disabled={isLoading || !inputSmiles}
                className="w-full"
              >
                {isLoading && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                Generate BRICS Analogs
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="validate" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>SMILES Validation</CardTitle>
              <CardDescription>Validate and canonicalize SMILES strings</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <Button
                onClick={handleValidateSmiles}
                disabled={isLoading || !inputSmiles}
                className="w-full"
              >
                {isLoading && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                Validate SMILES
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="full" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Full Integrated Pipeline</CardTitle>
              <CardDescription>
                Complete workflow: canonicalize → generate → screen → patent check
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <Button
                onClick={handleFullPipeline}
                disabled={isLoading || !inputSmiles}
                className="w-full"
              >
                {isLoading && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                Run Full Pipeline
              </Button>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>

      {results && (
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0">
            <div>
              <CardTitle>Results</CardTitle>
              <CardDescription>
                {results.success ? 'Pipeline executed successfully' : 'Error during execution'}
              </CardDescription>
            </div>
            {results.success && (
              <Button
                variant="outline"
                size="sm"
                onClick={downloadResults}
                className="gap-2"
              >
                <Download className="h-4 w-4" />
                Download JSON
              </Button>
            )}
          </CardHeader>
          <CardContent className="space-y-4">
            {results.success ? (
              <div className="space-y-4">
                {results.data?.canonical_smiles && (
                  <div className="rounded-lg bg-green-50 p-4 dark:bg-green-950">
                    <div className="flex gap-2">
                      <CheckCircle2 className="h-5 w-5 text-green-600 dark:text-green-400" />
                      <div>
                        <p className="font-semibold text-green-900 dark:text-green-100">
                          Canonical SMILES
                        </p>
                        <p className="text-sm text-green-800 dark:text-green-200 font-mono break-all">
                          {results.data.canonical_smiles}
                        </p>
                      </div>
                    </div>
                  </div>
                )}

                {results.data?.hits && results.data.hits.length > 0 && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">Similar Compounds ({results.data.hits.length})</h3>
                    <div className="max-h-96 overflow-y-auto space-y-2">
                      {results.data.hits.map((hit: any, idx: number) => (
                        <div key={idx} className="rounded border p-3 text-sm">
                          <div className="flex justify-between">
                            <span className="font-mono text-xs">{hit.smiles}</span>
                            <span className="font-semibold">{(hit.tanimoto * 100).toFixed(1)}%</span>
                          </div>
                          <div className="mt-1 text-xs text-gray-600 dark:text-gray-400">
                            CID: {hit.cid} | MW: {hit.mw?.toFixed(2)} | {hit.flag}
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {results.data?.analogs && results.data.analogs.length > 0 && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">Generated Analogs ({results.data.analogs.length})</h3>
                    <div className="max-h-96 overflow-y-auto space-y-2">
                      {results.data.analogs.map((analog: any, idx: number) => (
                        <div key={idx} className="rounded border p-3 text-sm">
                          <div className="flex justify-between">
                            <span className="font-mono text-xs">{analog.smiles}</span>
                            <span className="font-semibold">{(analog.tanimoto * 100).toFixed(1)}%</span>
                          </div>
                          <div className="mt-1 text-xs text-gray-600 dark:text-gray-400">
                            MW: {analog.mw?.toFixed(2)}
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {results.data?.master_list && results.data.master_list.length > 0 && (
                  <div className="space-y-2">
                    <h3 className="font-semibold">
                      Patent-Free Candidates ({results.data.master_list.length})
                    </h3>
                    <div className="max-h-96 overflow-y-auto space-y-2">
                      {results.data.master_list.map((candidate: any, idx: number) => (
                        <div key={idx} className="rounded border p-3 text-sm">
                          <div className="flex justify-between items-start">
                            <span className="font-mono text-xs flex-1">{candidate.smiles}</span>
                            <span className="ml-2 whitespace-nowrap text-xs font-semibold">
                              {candidate.flag}
                            </span>
                          </div>
                          <div className="mt-1 text-xs text-gray-600 dark:text-gray-400">
                            Similarity: {(candidate.tanimoto * 100).toFixed(1)}% | MW: {candidate.mw?.toFixed(2)}
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            ) : (
              <div className="rounded-lg bg-red-50 p-4 dark:bg-red-950">
                <div className="flex gap-2">
                  <AlertCircle className="h-5 w-5 text-red-600 dark:text-red-400" />
                  <div>
                    <p className="font-semibold text-red-900 dark:text-red-100">Error</p>
                    <p className="text-sm text-red-800 dark:text-red-200">{results.error}</p>
                  </div>
                </div>
              </div>
            )}
          </CardContent>
        </Card>
      )}
    </div>
  );
}
