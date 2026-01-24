import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Badge } from "@/components/ui/badge";
import { Loader2, Beaker, Activity, Skull, Pill } from "lucide-react";
import { toast } from "sonner";

export default function CompoundTesting() {
  const [selectedAnalog, setSelectedAnalog] = useState<number | null>(null);
  const [activeTest, setActiveTest] = useState<string | null>(null);

  const { data: analogs, isLoading: analogsLoading } = trpc.analog.list.useQuery({
    limit: 100,
    offset: 0,
  });

  const admetMutation = trpc.analog.runADMET.useMutation({
    onSuccess: () => {
      toast.success("ADMET analysis completed successfully");
      setActiveTest(null);
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

  // Toxicity mutation temporarily disabled
  // const toxicityMutation = trpc.analog.runToxicity.useMutation({
  //   onSuccess: () => {
  //     toast.success("Toxicity prediction completed successfully");
  //     setActiveTest(null);
  //   },
  //   onError: (error: any) => {
  //     toast.error(`Toxicity prediction failed: ${error.message}`);
  //     setActiveTest(null);
  //   },
  // });

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
    // toxicityMutation.mutate({
    //   analogId: selectedAnalog,
    // });
    toast.info("Toxicity testing temporarily disabled");
    setActiveTest(null);
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
              <TabsList className="grid w-full grid-cols-4">
                <TabsTrigger value="admet">ADMET</TabsTrigger>
                <TabsTrigger value="docking">Docking</TabsTrigger>
                <TabsTrigger value="toxicity">Toxicity</TabsTrigger>
                <TabsTrigger value="pkpd">PK/PD</TabsTrigger>
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
                  </div>
                </div>
              </TabsContent>
            </Tabs>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
