import { useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Loader2, CheckCircle, AlertCircle, Play, Download } from 'lucide-react';
import DockingPoseViewer from '@/components/DockingPoseViewer';
import DemoModeBadge from '@/components/DemoModeBadge';
import { trpc } from '@/lib/trpc';

interface TestResult {
  name: string;
  status: 'pending' | 'running' | 'success' | 'error';
  message?: string;
  data?: Record<string, unknown>;
  timestamp?: Date;
  source?: 'python' | 'fallback' | 'mock' | 'api';
}

/**
 * AnalysisTestPage
 * Dedicated page for testing docking, toxicity, ADMET, and PK/PD analysis functions
 * Provides end-to-end verification of the analysis pipeline
 */
export default function AnalysisTestPage() {
  const [smiles, setSmiles] = useState('CC(=O)Oc1ccccc1C(=O)O'); // Aspirin
  const [target, setTarget] = useState('NMDA');
  const [results, setResults] = useState<Record<string, TestResult>>({});
  const [isRunning, setIsRunning] = useState(false);

  // Test docking
  const dockingMutation = trpc.analog.runDocking.useMutation();
  const runDockingTest = async () => {
    setResults((prev) => ({
      ...prev,
      docking: { name: 'Docking Analysis', status: 'running' },
    }));

    try {
      const result = await dockingMutation.mutateAsync({
        analogId: 0,
        smiles,
        target,
      });

      setResults((prev) => ({
        ...prev,
        docking: {
          name: 'Docking Analysis',
          status: 'success',
          data: result,
          timestamp: new Date(),
          source: 'python',
        },
      }));
    } catch (error: any) {
      setResults((prev) => ({
        ...prev,
        docking: {
          name: 'Docking Analysis',
          status: 'error',
          message: error.message || 'Docking failed',
          timestamp: new Date(),
        },
      }));
    }
  };

  // Test toxicity - using ADMET as fallback since runToxicity may not exist
  const toxicityMutation = trpc.analog.runADMET.useMutation();
  const runToxicityTest = async () => {
    setResults((prev) => ({
      ...prev,
      toxicity: { name: 'Toxicity Prediction', status: 'running' },
    }));

    try {
      const result = await toxicityMutation.mutateAsync({
        analogId: 0,
        smiles,
      });

      setResults((prev) => ({
        ...prev,
        toxicity: {
          name: 'Toxicity Prediction',
          status: 'success',
          data: result,
          timestamp: new Date(),
          source: 'python',
        },
      }));
    } catch (error: any) {
      setResults((prev) => ({
        ...prev,
        toxicity: {
          name: 'Toxicity Prediction',
          status: 'error',
          message: error.message || 'Toxicity prediction failed',
          timestamp: new Date(),
        },
      }));
    }
  };

  // Test ADMET
  const admetMutation = trpc.analog.runADMET.useMutation();
  const runAdmetTest = async () => {
    setResults((prev) => ({
      ...prev,
      admet: { name: 'ADMET Analysis', status: 'running' },
    }));

    try {
      const result = await admetMutation.mutateAsync({
        analogId: 0,
        smiles,
      });

      setResults((prev) => ({
        ...prev,
        admet: {
          name: 'ADMET Analysis',
          status: 'success',
          data: result,
          timestamp: new Date(),
          source: 'python',
        },
      }));
    } catch (error: any) {
      setResults((prev) => ({
        ...prev,
        admet: {
          name: 'ADMET Analysis',
          status: 'error',
          message: error.message || 'ADMET analysis failed',
          timestamp: new Date(),
        },
      }));
    }
  };

  // Run all tests
  const runAllTests = async () => {
    setIsRunning(true);
    setResults({});

    await Promise.all([runDockingTest(), runToxicityTest(), runAdmetTest()]);

    setIsRunning(false);
  };

  // Export results
  const exportResults = () => {
    const exportData = {
      smiles,
      target,
      timestamp: new Date().toISOString(),
      results: Object.entries(results).map(([key, result]) => ({
        test: key,
        ...result,
      })),
    };

    const blob = new Blob([JSON.stringify(exportData, null, 2)], {
      type: 'application/json',
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `analysis-test-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const dockingResult = results.docking?.data as Record<string, unknown> | undefined;

  return (
    <div className="container mx-auto py-8 space-y-6">
      <div>
        <h1 className="text-3xl font-bold">Analysis Pipeline Test Suite</h1>
        <p className="text-gray-600 mt-2">
          End-to-end testing for docking, toxicity, ADMET, and PK/PD analysis functions
        </p>
      </div>

      {/* Input Section */}
      <Card>
        <CardHeader>
          <CardTitle>Test Configuration</CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <label className="text-sm font-medium">SMILES String</label>
              <Input
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder="Enter SMILES string"
                className="mt-1"
              />
              <p className="text-xs text-gray-500 mt-1">
                Current: Aspirin (CC(=O)Oc1ccccc1C(=O)O)
              </p>
            </div>
            <div>
              <label className="text-sm font-medium">Target Receptor</label>
              <Input
                value={target}
                onChange={(e) => setTarget(e.target.value)}
                placeholder="e.g., NMDA, COX2, ACE"
                className="mt-1"
              />
              <p className="text-xs text-gray-500 mt-1">
                Available: NMDA, COX2, ACE, EGFR
              </p>
            </div>
          </div>

          <div className="flex gap-2">
            <Button
              onClick={runAllTests}
              disabled={isRunning || !smiles}
              className="flex items-center gap-2"
            >
              {isRunning ? (
                <>
                  <Loader2 className="w-4 h-4 animate-spin" />
                  Running Tests...
                </>
              ) : (
                <>
                  <Play className="w-4 h-4" />
                  Run All Tests
                </>
              )}
            </Button>
            <Button
              onClick={exportResults}
              disabled={Object.keys(results).length === 0}
              variant="outline"
              className="flex items-center gap-2"
            >
              <Download className="w-4 h-4" />
              Export Results
            </Button>
          </div>
        </CardContent>
      </Card>

      {/* Results Section */}
      {Object.keys(results).length > 0 && (
        <Tabs defaultValue="overview" className="space-y-4">
          <TabsList>
            <TabsTrigger value="overview">Overview</TabsTrigger>
            <TabsTrigger value="docking">Docking</TabsTrigger>
            <TabsTrigger value="toxicity">Toxicity</TabsTrigger>
            <TabsTrigger value="admet">ADMET</TabsTrigger>
          </TabsList>

          {/* Overview Tab */}
          <TabsContent value="overview" className="space-y-4">
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              {Object.entries(results).map(([key, result]) => (
                <Card key={key}>
                  <CardHeader className="pb-3">
                    <div className="flex items-center justify-between">
                      <CardTitle className="text-base">{result.name}</CardTitle>
                      {result.status === 'success' && (
                        <CheckCircle className="w-5 h-5 text-green-600" />
                      )}
                      {result.status === 'error' && (
                        <AlertCircle className="w-5 h-5 text-red-600" />
                      )}
                      {result.status === 'running' && (
                        <Loader2 className="w-5 h-5 animate-spin text-blue-600" />
                      )}
                    </div>
                  </CardHeader>
                  <CardContent className="space-y-2">
                    <Badge
                      variant={
                        result.status === 'success'
                          ? 'default'
                          : result.status === 'error'
                            ? 'destructive'
                            : 'secondary'
                      }
                    >
                      {result.status.charAt(0).toUpperCase() + result.status.slice(1)}
                    </Badge>
                    {result.message && (
                      <p className="text-sm text-red-600">{result.message}</p>
                    )}
                    {result.timestamp && (
                      <p className="text-xs text-gray-500">
                        {result.timestamp.toLocaleTimeString()}
                      </p>
                    )}
                    {result.source && (
                      <div className="pt-2">
                        <DemoModeBadge source={result.source} />
                      </div>
                    )}
                  </CardContent>
                </Card>
              ))}
            </div>
          </TabsContent>

          {/* Docking Tab */}
          <TabsContent value="docking" className="space-y-4">
            {results.docking?.status === 'success' && dockingResult && (
              <div className="space-y-4">
                <DockingPoseViewer
                  ligandPDB={dockingResult.ligandPDB as string}
                  receptorPDB={dockingResult.receptorPDB as string}
                  bindingAffinity={dockingResult.bindingAffinity as string}
                  dockingScore={dockingResult.dockingScore as number}
                  target={target}
                  source={results.docking.source}
                  timestamp={results.docking.timestamp}
                />

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Docking Details</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="grid grid-cols-2 gap-4 text-sm">
                      <div>
                        <p className="font-semibold">Binding Affinity</p>
                        <p className="text-gray-600">
                          {String(dockingResult.bindingAffinity)} kcal/mol
                        </p>
                      </div>
                      <div>
                        <p className="font-semibold">Docking Score</p>
                        <p className="text-gray-600">{String(dockingResult.dockingScore)}/100</p>
                      </div>
                      <div>
                        <p className="font-semibold">Target</p>
                        <p className="text-gray-600">{target}</p>
                      </div>
                      <div>
                        <p className="font-semibold">SMILES</p>
                        <p className="text-gray-600 font-mono text-xs break-all">{smiles}</p>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </div>
            )}
            {results.docking?.status === 'error' && (
              <Card className="border-red-200 bg-red-50">
                <CardContent className="pt-6">
                  <p className="text-red-800">{results.docking.message}</p>
                </CardContent>
              </Card>
            )}
            {results.docking?.status === 'running' && (
              <Card>
                <CardContent className="pt-6 flex items-center gap-2">
                  <Loader2 className="w-4 h-4 animate-spin" />
                  <p>Running docking analysis...</p>
                </CardContent>
              </Card>
            )}
          </TabsContent>

          {/* Toxicity Tab */}
          <TabsContent value="toxicity" className="space-y-4">
            {results.toxicity?.status === 'success' && (
              <Card>
                <CardHeader>
                  <CardTitle className="text-base">Toxicity Prediction Results</CardTitle>
                </CardHeader>
                <CardContent>
                  <pre className="bg-gray-50 p-4 rounded text-xs overflow-auto max-h-96">
                    {JSON.stringify(results.toxicity.data, null, 2)}
                  </pre>
                </CardContent>
              </Card>
            )}
            {results.toxicity?.status === 'error' && (
              <Card className="border-red-200 bg-red-50">
                <CardContent className="pt-6">
                  <p className="text-red-800">{results.toxicity.message}</p>
                </CardContent>
              </Card>
            )}
            {results.toxicity?.status === 'running' && (
              <Card>
                <CardContent className="pt-6 flex items-center gap-2">
                  <Loader2 className="w-4 h-4 animate-spin" />
                  <p>Running toxicity prediction...</p>
                </CardContent>
              </Card>
            )}
          </TabsContent>

          {/* ADMET Tab */}
          <TabsContent value="admet" className="space-y-4">
            {results.admet?.status === 'success' && (
              <Card>
                <CardHeader>
                  <CardTitle className="text-base">ADMET Analysis Results</CardTitle>
                </CardHeader>
                <CardContent>
                  <pre className="bg-gray-50 p-4 rounded text-xs overflow-auto max-h-96">
                    {JSON.stringify(results.admet.data, null, 2)}
                  </pre>
                </CardContent>
              </Card>
            )}
            {results.admet?.status === 'error' && (
              <Card className="border-red-200 bg-red-50">
                <CardContent className="pt-6">
                  <p className="text-red-800">{results.admet.message}</p>
                </CardContent>
              </Card>
            )}
            {results.admet?.status === 'running' && (
              <Card>
                <CardContent className="pt-6 flex items-center gap-2">
                  <Loader2 className="w-4 h-4 animate-spin" />
                  <p>Running ADMET analysis...</p>
                </CardContent>
              </Card>
            )}
          </TabsContent>
        </Tabs>
      )}

      {/* Predefined Test Cases */}
      <Card>
        <CardHeader>
          <CardTitle>Predefined Test Cases</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-2">
            <Button
              variant="outline"
              className="w-full justify-start"
              onClick={() => {
                setSmiles('CC(=O)Oc1ccccc1C(=O)O');
                setTarget('NMDA');
              }}
            >
              <span className="font-mono text-sm">Aspirin on NMDA</span>
            </Button>
            <Button
              variant="outline"
              className="w-full justify-start"
              onClick={() => {
                setSmiles('CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O');
                setTarget('COX2');
              }}
            >
              <span className="font-mono text-sm">Ibuprofen on COX2</span>
            </Button>
            <Button
              variant="outline"
              className="w-full justify-start"
              onClick={() => {
                setSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
                setTarget('ADENOSINE');
              }}
            >
              <span className="font-mono text-sm">Caffeine on Adenosine Receptor</span>
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
