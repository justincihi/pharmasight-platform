import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Checkbox } from "@/components/ui/checkbox";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Download, Play, CheckCircle2, XCircle, Loader2, FlaskConical } from "lucide-react";
import { toast } from "sonner";
import { BatchAdmetPanel } from "@/components/BatchAdmetPanel";

interface BatchTestResult {
  analogId: number;
  compoundName: string;
  status: "pending" | "running" | "completed" | "failed";
  results?: {
    admet?: any;
    docking?: any;
    toxicity?: any;
    pkpd?: any;
  };
  error?: string;
}

export default function BatchAnalysis() {
  const [selectedAnalogs, setSelectedAnalogs] = useState<number[]>([]);
  const [selectedTests, setSelectedTests] = useState<string[]>([]);
  const [batchResults, setBatchResults] = useState<BatchTestResult[]>([]);
  const [isRunning, setIsRunning] = useState(false);

  const { data: analogs, isLoading } = trpc.analog.list.useQuery({ limit: 1000, offset: 0 });
  const runBatchMutation = trpc.cheminformatics.runBatch.useMutation();

  const testOptions = [
    { id: "admet", label: "ADMET Prediction", description: "Drug metabolism and toxicity" },
    { id: "docking", label: "Molecular Docking", description: "Protein-ligand binding" },
    { id: "toxicity", label: "Toxicity Assessment", description: "Safety profile analysis" },
    { id: "pkpd", label: "PK/PD Simulation", description: "Pharmacokinetics modeling" },
  ];

  const toggleAnalog = (id: number) => {
    setSelectedAnalogs((prev) =>
      prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]
    );
  };

  const toggleTest = (testId: string) => {
    setSelectedTests((prev) =>
      prev.includes(testId) ? prev.filter((x) => x !== testId) : [...prev, testId]
    );
  };

  const selectAll = () => {
    if (analogs) {
      setSelectedAnalogs(analogs.map((a: any) => a.id));
    }
  };

  const deselectAll = () => {
    setSelectedAnalogs([]);
  };

  const runBatchAnalysis = async () => {
    if (selectedAnalogs.length === 0 || selectedTests.length === 0) {
      toast.error("Please select at least one analog and one test");
      return;
    }

    setIsRunning(true);
    
    // Initialize results
    const initialResults: BatchTestResult[] = selectedAnalogs.map((id) => ({
      analogId: id,
      compoundName: analogs?.find((a: any) => a.id === id)?.compoundName || "",
      status: "pending",
    }));
    setBatchResults(initialResults);

    try {
      const result = await runBatchMutation.mutateAsync({
        analogIds: selectedAnalogs,
        tests: selectedTests,
      });

      // Update results with completed data
      setBatchResults(result.results);
      toast.success(`Batch analysis complete: ${result.completed}/${result.total} successful`);
    } catch (error) {
      toast.error("Batch analysis failed");
      console.error(error);
    } finally {
      setIsRunning(false);
    }
  };

  const downloadResults = () => {
    const csv = generateCSV(batchResults);
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `batch_analysis_${new Date().toISOString()}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
    toast.success("Results downloaded");
  };

  const generateCSV = (results: BatchTestResult[]) => {
    const headers = ["Compound", "Status", "ADMET", "Docking", "Toxicity", "PK/PD"];
    const rows = results.map((r) => [
      r.compoundName,
      r.status,
      r.results?.admet ? "✓" : "-",
      r.results?.docking ? "✓" : "-",
      r.results?.toxicity ? "✓" : "-",
      r.results?.pkpd ? "✓" : "-",
    ]);
    return [headers, ...rows].map((row) => row.join(",")).join("\\n");
  };

  const completedCount = batchResults.filter((r) => r.status === "completed").length;
  const failedCount = batchResults.filter((r) => r.status === "failed").length;
  const progress = batchResults.length > 0 
    ? ((completedCount + failedCount) / batchResults.length) * 100 
    : 0;

  if (isLoading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <Loader2 className="w-8 h-8 animate-spin text-blue-600" />
      </div>
    );
  }

  return (
    <div className="container mx-auto py-8">
      <div className="mb-6">
        <h1 className="text-3xl font-bold mb-2">Batch Analysis</h1>
        <p className="text-muted-foreground">
          Run cheminformatics analyses on multiple analogs simultaneously
        </p>
      </div>

      <Tabs defaultValue="ml-admet" className="space-y-6">
        <TabsList>
          <TabsTrigger value="ml-admet" className="gap-2">
            <FlaskConical className="h-4 w-4" />
            ML ADMET Screening
          </TabsTrigger>
          <TabsTrigger value="classic">Classic Batch</TabsTrigger>
        </TabsList>

        <TabsContent value="ml-admet">
          <BatchAdmetPanel />
        </TabsContent>

        <TabsContent value="classic">
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Selection Panel */}
        <div className="lg:col-span-2 space-y-6">
          {/* Analog Selection */}
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <CardTitle>Select Analogs ({selectedAnalogs.length} selected)</CardTitle>
                <div className="flex gap-2">
                  <Button size="sm" variant="outline" onClick={selectAll}>
                    Select All
                  </Button>
                  <Button size="sm" variant="outline" onClick={deselectAll}>
                    Clear
                  </Button>
                </div>
              </div>
            </CardHeader>
            <CardContent>
              <div className="max-h-96 overflow-y-auto space-y-2">
                {analogs?.map((analog: any) => (
                  <div
                    key={analog.id}
                    className="flex items-center space-x-3 p-3 border rounded hover:bg-gray-50 cursor-pointer"
                    onClick={() => toggleAnalog(analog.id)}
                  >
                    <Checkbox
                      checked={selectedAnalogs.includes(analog.id)}
                      onCheckedChange={() => toggleAnalog(analog.id)}
                    />
                    <div className="flex-1">
                      <p className="font-medium text-sm">{analog.compoundName}</p>
                      <p className="text-xs text-gray-500">{analog.parentCompound}</p>
                    </div>
                    <Badge variant="secondary">{analog.confidenceScore}%</Badge>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>

          {/* Test Selection */}
          <Card>
            <CardHeader>
              <CardTitle>Select Tests ({selectedTests.length} selected)</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                {testOptions.map((test) => (
                  <div
                    key={test.id}
                    className="flex items-start space-x-3 p-3 border rounded hover:bg-gray-50 cursor-pointer"
                    onClick={() => toggleTest(test.id)}
                  >
                    <Checkbox
                      checked={selectedTests.includes(test.id)}
                      onCheckedChange={() => toggleTest(test.id)}
                    />
                    <div>
                      <p className="font-medium text-sm">{test.label}</p>
                      <p className="text-xs text-gray-500">{test.description}</p>
                    </div>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Control Panel */}
        <div className="space-y-6">
          <Card>
            <CardHeader>
              <CardTitle>Batch Control</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <p className="text-sm text-gray-600">
                  <strong>{selectedAnalogs.length}</strong> analogs selected
                </p>
                <p className="text-sm text-gray-600">
                  <strong>{selectedTests.length}</strong> tests selected
                </p>
                <p className="text-sm text-gray-600">
                  Total operations: <strong>{selectedAnalogs.length * selectedTests.length}</strong>
                </p>
              </div>

              <Button
                onClick={runBatchAnalysis}
                disabled={isRunning || selectedAnalogs.length === 0 || selectedTests.length === 0}
                className="w-full"
              >
                {isRunning ? (
                  <>
                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <Play className="w-4 h-4 mr-2" />
                    Run Batch Analysis
                  </>
                )}
              </Button>

              {batchResults.length > 0 && (
                <Button
                  onClick={downloadResults}
                  variant="outline"
                  className="w-full"
                >
                  <Download className="w-4 h-4 mr-2" />
                  Download Results (CSV)
                </Button>
              )}
            </CardContent>
          </Card>

          {/* Progress Panel */}
          {batchResults.length > 0 && (
            <Card>
              <CardHeader>
                <CardTitle>Progress</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <Progress value={progress} className="w-full" />
                <div className="space-y-2">
                  <div className="flex items-center justify-between text-sm">
                    <span className="flex items-center gap-2">
                      <CheckCircle2 className="w-4 h-4 text-green-600" />
                      Completed
                    </span>
                    <span className="font-medium">{completedCount}</span>
                  </div>
                  <div className="flex items-center justify-between text-sm">
                    <span className="flex items-center gap-2">
                      <XCircle className="w-4 h-4 text-red-600" />
                      Failed
                    </span>
                    <span className="font-medium">{failedCount}</span>
                  </div>
                  <div className="flex items-center justify-between text-sm">
                    <span className="flex items-center gap-2">
                      <Loader2 className="w-4 h-4 text-blue-600 animate-spin" />
                      Pending
                    </span>
                    <span className="font-medium">
                      {batchResults.length - completedCount - failedCount}
                    </span>
                  </div>
                </div>
              </CardContent>
            </Card>
          )}
        </div>
      </div>

      {/* Results Table */}
      {batchResults.length > 0 && (
        <Card className="mt-6">
          <CardHeader>
            <CardTitle>Results</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead>
                  <tr className="border-b">
                    <th className="text-left p-3 text-sm font-medium">Compound</th>
                    <th className="text-left p-3 text-sm font-medium">Status</th>
                    <th className="text-center p-3 text-sm font-medium">ADMET</th>
                    <th className="text-center p-3 text-sm font-medium">Docking</th>
                    <th className="text-center p-3 text-sm font-medium">Toxicity</th>
                    <th className="text-center p-3 text-sm font-medium">PK/PD</th>
                  </tr>
                </thead>
                <tbody>
                  {batchResults.map((result) => (
                    <tr key={result.analogId} className="border-b hover:bg-gray-50">
                      <td className="p-3 text-sm">{result.compoundName}</td>
                      <td className="p-3">
                        <Badge
                          variant={
                            result.status === "completed"
                              ? "default"
                              : result.status === "failed"
                              ? "destructive"
                              : "secondary"
                          }
                        >
                          {result.status}
                        </Badge>
                      </td>
                      <td className="text-center p-3">
                        {result.results?.admet ? (
                          <CheckCircle2 className="w-5 h-5 text-green-600 mx-auto" />
                        ) : (
                          <span className="text-gray-400">-</span>
                        )}
                      </td>
                      <td className="text-center p-3">
                        {result.results?.docking ? (
                          <CheckCircle2 className="w-5 h-5 text-green-600 mx-auto" />
                        ) : (
                          <span className="text-gray-400">-</span>
                        )}
                      </td>
                      <td className="text-center p-3">
                        {result.results?.toxicity ? (
                          <CheckCircle2 className="w-5 h-5 text-green-600 mx-auto" />
                        ) : (
                          <span className="text-gray-400">-</span>
                        )}
                      </td>
                      <td className="text-center p-3">
                        {result.results?.pkpd ? (
                          <CheckCircle2 className="w-5 h-5 text-green-600 mx-auto" />
                        ) : (
                          <span className="text-gray-400">-</span>
                        )}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>
      )}
        </TabsContent>
      </Tabs>
    </div>
  );
}
