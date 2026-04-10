import React, { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Badge } from "@/components/ui/badge";
import { Loader2, Beaker, Activity, Skull, Pill } from "lucide-react";
import { toast } from "sonner";
import { DockingParametersPanel } from "@/components/DockingParametersPanel";
import { PDBUploadDialog } from "@/components/PDBUploadDialog";
import { ReceptorSelector, type SelectedReceptor } from "@/components/ReceptorSelector";
import { DockingResultsPanel } from "@/components/DockingResultsPanel";
import { DockingPoseViewer3D } from "@/components/DockingPoseViewer3D";
import { BatchKetamineTestingPanel } from "@/components/BatchKetamineTestingPanel";
import { ReceptorSelectivityPanel } from "@/components/ReceptorSelectivityPanel";
import { DockingResultsExportPanel } from "@/components/DockingResultsExportPanel";

export default function CompoundTesting() {
  const [selectedAnalog, setSelectedAnalog] = useState<number | null>(null);
  const [activeTest, setActiveTest] = useState<string | null>(null);
  const [pdbDialogOpen, setPdbDialogOpen] = useState(false);
  const [showDockingParams, setShowDockingParams] = useState(false);
  const [selectedReceptor, setSelectedReceptor] = useState<SelectedReceptor | null>(null);
  const [showReceptorSelector, setShowReceptorSelector] = useState(false);
  const [uploadedPdbFile, setUploadedPdbFile] = useState<File | null>(null);

  const { data: analogs, isLoading: analogsLoading } = trpc.analog.list.useQuery({
    limit: 100,
    offset: 0,
  });

  // Fetch test results for selected analog
  const { data: testResults, refetch: refetchResults } = trpc.analog.getTestResults.useQuery(
    { analogId: selectedAnalog! },
    { enabled: !!selectedAnalog }
  );

  const admetMutation = trpc.analog.runADMET.useMutation({
    onSuccess: () => {
      toast.success("ADMET analysis completed successfully");
      setActiveTest(null);
      refetchResults();
    },
    onError: (error: any) => {
      toast.error(`ADMET analysis failed: ${error.message}`);
      setActiveTest(null);
    },
  });

  const dockingMutation = trpc.analog.runDocking.useMutation({
    onSuccess: () => {
      toast.success("Molecular docking completed successfully");
      setActiveTest(null);
    },
    onError: (error: any) => {
      toast.error(`Docking failed: ${error.message}`);
      setActiveTest(null);
    },
  });

  const toxicityMutation = trpc.advancedAnalysis.toxicity.useMutation({
    onSuccess: () => {
      toast.success("Toxicity prediction completed successfully");
      setActiveTest(null);
    },
    onError: (error: any) => {
      toast.error(`Toxicity prediction failed: ${error.message}`);
      setActiveTest(null);
    },
  });

  const pkpdMutation = trpc.cheminformatics.simulatePKPD.useMutation({
    onSuccess: () => {
      toast.success("PK/PD simulation completed successfully");
      setActiveTest(null);
    },
    onError: (error: any) => {
      toast.error(`PK/PD simulation failed: ${error.message}`);
      setActiveTest(null);
    },
  });

  const runADMET = () => {
    if (!selectedAnalog) {
      toast.error("Please select a compound first");
      return;
    }

    const analog = analogs?.find((a: any) => a.id === selectedAnalog);
    if (!analog) return;

    setActiveTest("admet");
    admetMutation.mutate({
      analogId: selectedAnalog,
      smiles: analog.smiles,
    });
  };

  const runDocking = () => {
    if (!selectedAnalog) {
      toast.error("Please select a compound first");
      return;
    }

    const analog = analogs?.find((a: any) => a.id === selectedAnalog);
    if (!analog) return;

    setActiveTest("docking");
    dockingMutation.mutate({
      analogId: selectedAnalog,
      smiles: analog.smiles,
      target: "5HT2A", // Default receptor
    });
  };

  const runToxicity = () => {
    if (!selectedAnalog) {
      toast.error("Please select a compound first");
      return;
    }

    const analog = analogs?.find((a: any) => a.id === selectedAnalog);
    if (!analog) return;

    setActiveTest("toxicity");
    toxicityMutation.mutate({
      smiles: analog.smiles,
    });
  };

  const runPKPD = () => {
    if (!selectedAnalog) {
      toast.error("Please select a compound first");
      return;
    }

    const analog = analogs?.find((a: any) => a.id === selectedAnalog);
    if (!analog) return;

    setActiveTest("pkpd");
    pkpdMutation.mutate({
      smiles: analog.smiles,
      dose: 100,
      route: "oral",
    });
  };

  const selectedAnalogData = analogs?.find((a: any) => a.id === selectedAnalog);

  return (
    <div className="container mx-auto py-8">
      <div className="mb-8">
        <h1 className="text-3xl font-bold mb-2">Compound Testing Interface</h1>
        <p className="text-muted-foreground">
          Run cheminformatics analyses on discovered analogs
        </p>
      </div>

      <div className="grid gap-6 md:grid-cols-[300px_1fr]">
        {/* Compound Selector */}
        <Card>
          <CardHeader>
            <CardTitle>Select Compound</CardTitle>
            <CardDescription>Choose an analog to test</CardDescription>
          </CardHeader>
          <CardContent>
            {analogsLoading ? (
              <div className="flex items-center justify-center py-8">
                <Loader2 className="h-6 w-6 animate-spin" />
              </div>
            ) : (
              <Select
                value={selectedAnalog?.toString() || ""}
                onValueChange={(value) => setSelectedAnalog(parseInt(value))}
              >
                <SelectTrigger>
                  <SelectValue placeholder="Select compound..." />
                </SelectTrigger>
                <SelectContent>
                  {analogs?.map((analog: any) => (
                    <SelectItem key={analog.id} value={analog.id.toString()}>
                      {analog.compoundName}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            )}

            {selectedAnalogData && (
              <div className="mt-4 space-y-2">
                <div className="text-sm">
                  <span className="font-medium">Parent:</span>{" "}
                  {selectedAnalogData.parentCompound}
                </div>
                <div className="text-sm">
                  <span className="font-medium">Confidence:</span>{" "}
                  <Badge variant="secondary">
                    {selectedAnalogData.confidenceScore}%
                  </Badge>
                </div>
                <div className="text-sm">
                  <span className="font-medium">SMILES:</span>
                  <code className="block mt-1 p-2 bg-muted rounded text-xs break-all">
                    {selectedAnalogData.smiles}
                  </code>
                </div>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Testing Interface */}
        <Card>
          <CardHeader>
            <CardTitle>Cheminformatics Analyses</CardTitle>
            <CardDescription>
              Run various tests on the selected compound
            </CardDescription>
          </CardHeader>
          <CardContent>
            <Tabs defaultValue="admet" className="w-full">
              <TabsList className="grid w-full grid-cols-7">
                <TabsTrigger value="admet">ADMET</TabsTrigger>
                <TabsTrigger value="docking">Docking</TabsTrigger>
                <TabsTrigger value="toxicity">Toxicity</TabsTrigger>
                <TabsTrigger value="pkpd">PK/PD</TabsTrigger>
                <TabsTrigger value="batch">Batch Test</TabsTrigger>
                <TabsTrigger value="selectivity">Selectivity</TabsTrigger>
                <TabsTrigger value="export">Export</TabsTrigger>
              </TabsList>

              <TabsContent value="admet" className="space-y-4">
                <div className="flex items-start gap-4">
                  <Beaker className="h-8 w-8 text-blue-500 mt-1" />
                  <div className="flex-1">
                    <h3 className="font-semibold mb-2">ADMET Prediction</h3>
                    <p className="text-sm text-muted-foreground mb-4">
                      Predict Absorption, Distribution, Metabolism, Excretion, and
                      Toxicity properties using advanced ML models.
                    </p>
                    <Button
                      onClick={runADMET}
                      disabled={!selectedAnalog || activeTest === "admet"}
                    >
                      {activeTest === "admet" ? (
                        <>
                          <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                          Running Analysis...
                        </>
                      ) : (
                        "Run ADMET Analysis"
                      )}
                    </Button>

                    {/* Display ADMET Results */}
                    {testResults && testResults.filter((r: any) => r.testType === 'admet').length > 0 && (
                      <div className="mt-6 space-y-4">
                        {testResults
                          .filter((r: any) => r.testType === 'admet')
                          .sort((a: any, b: any) => new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime())
                          .slice(0, 1)
                          .map((result: any) => {
                            const data = JSON.parse(result.results || '{}');
                            const toxicity = data.toxicity_profile || {};
                            const sa = data.synthetic_accessibility || {};
                            
                            return (
                              <Card key={result.id} className="mt-4">
                                <CardHeader>
                                  <CardTitle className="text-lg">Latest ADMET Results</CardTitle>
                                  <CardDescription>
                                    Analyzed {new Date(result.createdAt).toLocaleString()}
                                  </CardDescription>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                  {/* Toxicity Profile */}
                                  <div>
                                    <h4 className="font-semibold mb-2">Toxicity Profile</h4>
                                    <div className="grid grid-cols-2 gap-3">
                                      {toxicity.hERG && (
                                        <div className="p-3 border rounded-lg">
                                          <div className="text-sm font-medium">hERG Risk</div>
                                          <div className="text-2xl font-bold">{toxicity.hERG.risk_score?.toFixed(1)}</div>
                                          <Badge variant={toxicity.hERG.risk_level === 'Low' ? 'default' : 'destructive'}>
                                            {toxicity.hERG.risk_level}
                                          </Badge>
                                        </div>
                                      )}
                                      {toxicity.hepatotoxicity && (
                                        <div className="p-3 border rounded-lg">
                                          <div className="text-sm font-medium">Hepatotoxicity</div>
                                          <div className="text-2xl font-bold">{toxicity.hepatotoxicity.risk_score?.toFixed(1)}</div>
                                          <Badge variant={toxicity.hepatotoxicity.risk_level === 'Low' ? 'default' : 'destructive'}>
                                            {toxicity.hepatotoxicity.risk_level}
                                          </Badge>
                                        </div>
                                      )}
                                      {toxicity.mutagenicity && (
                                        <div className="p-3 border rounded-lg">
                                          <div className="text-sm font-medium">Mutagenicity</div>
                                          <Badge variant={toxicity.mutagenicity.risk_level === 'Low' ? 'default' : 'destructive'}>
                                            {toxicity.mutagenicity.prediction}
                                          </Badge>
                                        </div>
                                      )}
                                      {toxicity.carcinogenicity && (
                                        <div className="p-3 border rounded-lg">
                                          <div className="text-sm font-medium">Carcinogenicity</div>
                                          <Badge variant={toxicity.carcinogenicity.risk_level === 'Low' ? 'default' : 'destructive'}>
                                            {toxicity.carcinogenicity.prediction}
                                          </Badge>
                                        </div>
                                      )}
                                    </div>
                                  </div>

                                  {/* Synthetic Accessibility */}
                                  {sa.sa_score && (
                                    <div>
                                      <h4 className="font-semibold mb-2">Synthetic Accessibility</h4>
                                      <div className="p-3 border rounded-lg">
                                        <div className="flex justify-between items-center">
                                          <div>
                                            <div className="text-sm text-muted-foreground">SA Score</div>
                                            <div className="text-2xl font-bold">{sa.sa_score?.toFixed(1)}</div>
                                          </div>
                                          <div className="text-right">
                                            <Badge>{sa.difficulty}</Badge>
                                            <div className="text-sm text-muted-foreground mt-1">{sa.estimated_steps}</div>
                                          </div>
                                        </div>
                                        <p className="text-sm text-muted-foreground mt-2">{sa.recommendation}</p>
                                      </div>
                                    </div>
                                  )}

                                  {/* Optimization Suggestions */}
                                  {data.optimization_suggestions && data.optimization_suggestions.length > 0 && (
                                    <div>
                                      <h4 className="font-semibold mb-2">Optimization Suggestions</h4>
                                      {data.optimization_suggestions.map((sug: any, idx: number) => (
                                        <div key={idx} className="p-3 border rounded-lg mb-2">
                                          <div className="font-medium">{sug.modification}</div>
                                          <p className="text-sm text-muted-foreground mt-1">{sug.rationale}</p>
                                        </div>
                                      ))}
                                    </div>
                                  )}
                                </CardContent>
                              </Card>
                            );
                          })}
                      </div>
                    )}
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="docking" className="space-y-4">
                <div className="flex items-start gap-4">
                  <Activity className="h-8 w-8 text-green-500 mt-1" />
                  <div className="flex-1">
                    <h3 className="font-semibold mb-2">Molecular Docking</h3>
                    <p className="text-sm text-muted-foreground mb-4">
                      Simulate binding interactions with target receptors using
                      AutoDock Vina to predict binding affinity and pose.
                    </p>
                    <div className="flex gap-2 mb-4 flex-wrap">
                      <Button
                        onClick={runDocking}
                        disabled={!selectedAnalog || activeTest === "docking"}
                      >
                        {activeTest === "docking" ? (
                          <>
                            <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                            Running Docking...
                          </>
                        ) : (
                          "Run Docking Simulation"
                        )}
                      </Button>
                      <Button
                        variant="outline"
                        onClick={() => setShowReceptorSelector(!showReceptorSelector)}
                      >
                        {showReceptorSelector ? "Hide" : "Select"} Receptor
                      </Button>
                      <Button
                        variant="outline"
                        onClick={() => setShowDockingParams(!showDockingParams)}
                      >
                        {showDockingParams ? "Hide" : "Show"} Parameters
                      </Button>
                      <Button
                        variant="outline"
                        onClick={() => setPdbDialogOpen(true)}
                      >
                        Upload PDB
                      </Button>
                    </div>
                    {selectedReceptor && (
                      <div className="mb-4 p-3 bg-blue-50 border border-blue-200 rounded">
                        <p className="text-sm font-medium text-blue-900">Receptor Selected:</p>
                        <p className="text-sm text-blue-800 mt-1">
                          {selectedReceptor.family} - {selectedReceptor.subtype} ({selectedReceptor.species})
                        </p>
                      </div>
                    )}
                    {showReceptorSelector && (
                      <div className="mt-6 border-t pt-6">
                        <ReceptorSelector
                          onReceptorSelect={(receptor) => {
                            setSelectedReceptor(receptor);
                            setShowReceptorSelector(false);
                            toast.success(`Receptor selected: ${receptor.family}`);
                          }}
                          onPdbUpload={(file) => {
                            setUploadedPdbFile(file);
                            toast.success(`PDB file loaded: ${file.name}`);
                          }}
                          allowPdbUpload={true}
                        />
                      </div>
                    )}
                    {showDockingParams && (
                      <div className="mt-6 border-t pt-6">
                        <DockingParametersPanel />
                      </div>
                    )}
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="toxicity" className="space-y-4">
                <div className="flex items-start gap-4">
                  <Skull className="h-8 w-8 text-red-500 mt-1" />
                  <div className="flex-1">
                    <h3 className="font-semibold mb-2">Toxicity Prediction</h3>
                    <p className="text-sm text-muted-foreground mb-4">
                      Assess potential toxic effects including acute toxicity,
                      genotoxicity, and carcinogenicity using ML models.
                    </p>
                    <Button
                      onClick={runToxicity}
                      disabled={!selectedAnalog || activeTest === "toxicity"}
                    >
                      {activeTest === "toxicity" ? (
                        <>
                          <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                          Predicting Toxicity...
                        </>
                      ) : (
                        "Run Toxicity Prediction"
                      )}
                    </Button>
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="pkpd" className="space-y-4">
                <div className="flex items-start gap-4">
                  <Pill className="h-8 w-8 text-purple-500 mt-1" />
                  <div className="flex-1">
                    <h3 className="font-semibold mb-2">PK/PD Simulation</h3>
                    <p className="text-sm text-muted-foreground mb-4">
                      Simulate pharmacokinetic and pharmacodynamic profiles to
                      predict drug concentration over time and therapeutic effects.
                    </p>
                    <Button
                      onClick={runPKPD}
                      disabled={!selectedAnalog || activeTest === "pkpd"}
                    >
                      {activeTest === "pkpd" ? (
                        <>
                          <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                          Running Simulation...
                        </>
                      ) : (
                        "Run PK/PD Simulation"
                      )}
                    </Button>
                    
                    {/* PK/PD Results Display */}
                    {testResults?.pkpd && testResults.pkpd.length > 0 && (
                      <div className="mt-6 space-y-4">
                        <h4 className="font-semibold">Latest PK/PD Results</h4>
                        {testResults.pkpd.map((result: any, idx: number) => (
                          <Card key={idx}>
                            <CardHeader>
                              <CardTitle className="text-lg">
                                PK/PD Profile - {result.dose}mg {result.route}
                              </CardTitle>
                              <CardDescription>
                                Analyzed {new Date(result.analyzedAt).toLocaleString()}
                              </CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-4">
                              {/* Pharmacokinetic Parameters */}
                              {result.pharmacokinetics && (
                                <div>
                                  <h4 className="font-semibold mb-3">Pharmacokinetic Parameters</h4>
                                  <div className="grid grid-cols-2 gap-3">
                                    {result.pharmacokinetics.absorption && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Absorption (Tmax)</div>
                                        <div className="text-lg font-bold">{result.pharmacokinetics.absorption.tmax?.toFixed(2)} h</div>
                                      </div>
                                    )}
                                    {result.pharmacokinetics.distribution && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Volume of Distribution</div>
                                        <div className="text-lg font-bold">{result.pharmacokinetics.distribution.vd?.toFixed(2)} L/kg</div>
                                      </div>
                                    )}
                                    {result.pharmacokinetics.elimination && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Half-life</div>
                                        <div className="text-lg font-bold">{result.pharmacokinetics.elimination.t_half?.toFixed(2)} h</div>
                                      </div>
                                    )}
                                    {result.pharmacokinetics.clearance && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Clearance</div>
                                        <div className="text-lg font-bold">{result.pharmacokinetics.clearance.cl?.toFixed(2)} L/h</div>
                                      </div>
                                    )}
                                  </div>
                                </div>
                              )}

                              {/* Pharmacodynamic Parameters */}
                              {result.pharmacodynamics && (
                                <div>
                                  <h4 className="font-semibold mb-3">Pharmacodynamic Parameters</h4>
                                  <div className="grid grid-cols-2 gap-3">
                                    {result.pharmacodynamics.potency && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Potency (EC50)</div>
                                        <div className="text-lg font-bold">{result.pharmacodynamics.potency.ec50?.toFixed(2)} nM</div>
                                      </div>
                                    )}
                                    {result.pharmacodynamics.efficacy && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Efficacy (Emax)</div>
                                        <div className="text-lg font-bold">{(result.pharmacodynamics.efficacy.emax * 100).toFixed(1)}%</div>
                                      </div>
                                    )}
                                    {result.pharmacodynamics.onset && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Onset Time</div>
                                        <div className="text-lg font-bold">{result.pharmacodynamics.onset.tmax?.toFixed(2)} h</div>
                                      </div>
                                    )}
                                    {result.pharmacodynamics.duration && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm text-muted-foreground">Duration</div>
                                        <div className="text-lg font-bold">{result.pharmacodynamics.duration.duration?.toFixed(2)} h</div>
                                      </div>
                                    )}
                                  </div>
                                </div>
                              )}

                              {/* Therapeutic Window */}
                              {result.therapeuticWindow && (
                                <div>
                                  <h4 className="font-semibold mb-3">Therapeutic Window</h4>
                                  <div className="p-3 border rounded-lg">
                                    <div className="flex justify-between items-center mb-2">
                                      <span className="text-sm text-muted-foreground">Therapeutic Index</span>
                                      <Badge variant={result.therapeuticWindow.index > 2 ? 'default' : 'destructive'}>
                                        {result.therapeuticWindow.index?.toFixed(2)}
                                      </Badge>
                                    </div>
                                    <p className="text-sm text-muted-foreground">{result.therapeuticWindow.interpretation}</p>
                                  </div>
                                </div>
                              )}

                              {/* Safety Assessment */}
                              {result.safetyAssessment && (
                                <div>
                                  <h4 className="font-semibold mb-3">Safety Assessment</h4>
                                  <div className="space-y-2">
                                    {result.safetyAssessment.hepatotoxicity && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm font-medium">Hepatotoxicity Risk</div>
                                        <Badge variant={result.safetyAssessment.hepatotoxicity.risk === 'Low' ? 'default' : 'destructive'}>
                                          {result.safetyAssessment.hepatotoxicity.risk}
                                        </Badge>
                                      </div>
                                    )}
                                    {result.safetyAssessment.nephrotoxicity && (
                                      <div className="p-3 border rounded-lg">
                                        <div className="text-sm font-medium">Nephrotoxicity Risk</div>
                                        <Badge variant={result.safetyAssessment.nephrotoxicity.risk === 'Low' ? 'default' : 'destructive'}>
                                          {result.safetyAssessment.nephrotoxicity.risk}
                                        </Badge>
                                      </div>
                                    )}
                                  </div>
                                </div>
                              )}
                            </CardContent>
                          </Card>
                        ))}
                      </div>
                    )}
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="batch" className="space-y-4">
                <BatchKetamineTestingPanel />
              </TabsContent>

              <TabsContent value="selectivity" className="space-y-4">
                <ReceptorSelectivityPanel results={[]} />
              </TabsContent>

              <TabsContent value="export" className="space-y-4">
                <DockingResultsExportPanel results={[]} />
              </TabsContent>
            </Tabs>
          </CardContent>
        </Card>
      </div>

      {/* PDB Upload Dialog */}
      <PDBUploadDialog
        open={pdbDialogOpen}
        onOpenChange={setPdbDialogOpen}
        onSuccess={() => {
          alert('✅ PDB file uploaded successfully');
        }}
      />
    </div>
  );
}
