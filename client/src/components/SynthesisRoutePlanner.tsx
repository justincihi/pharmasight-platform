import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { Alert, AlertDescription } from "@/components/ui/alert";
import { Loader2, FlaskConical, DollarSign, Clock, TrendingUp, ChevronRight, AlertCircle, GitCompare } from "lucide-react";

interface SynthesisStep {
  stepNumber: number;
  reaction: string;
  reagents: string[];
  conditions: string;
  yield: string;
  difficulty: "easy" | "moderate" | "difficult";
  estimatedCost: number;
  notes: string;
}

interface SynthesisRoute {
  routeId: string;
  targetSmiles: string;
  targetName: string;
  totalSteps: number;
  overallYield: string;
  totalCost: number;
  feasibilityScore: number;
  difficulty: "easy" | "moderate" | "difficult";
  estimatedTime: string;
  steps: SynthesisStep[];
  startingMaterials: Array<{
    name: string;
    smiles: string;
    availability: "commercial" | "synthesize";
    cost: number;
  }>;
  summary: string;
}

interface SynthesisRoutePlannerProps {
  smiles: string;
  compoundName: string;
}

const difficultyColors = {
  easy: "bg-green-100 text-green-800 border-green-300",
  moderate: "bg-yellow-100 text-yellow-800 border-yellow-300",
  difficult: "bg-red-100 text-red-800 border-red-300",
};

const difficultyLabels = {
  easy: "Easy",
  moderate: "Moderate",
  difficult: "Difficult",
};

export default function SynthesisRoutePlanner({ smiles, compoundName }: SynthesisRoutePlannerProps) {
  const [routes, setRoutes] = useState<SynthesisRoute[]>([]);
  const [selectedRoute, setSelectedRoute] = useState<number | null>(null);
  const [compareMode, setCompareMode] = useState(false);
  const [selectedRoutes, setSelectedRoutes] = useState<number[]>([]);

  const generateMutation = trpc.synthesis.generateRoutes.useMutation({
    onSuccess: (data) => {
      setRoutes(data);
      if (data.length > 0) {
        setSelectedRoute(0);
      }
    },
  });

  const handleGenerate = () => {
    generateMutation.mutate({
      smiles,
      compoundName,
      numRoutes: 2,
    });
  };

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-2xl font-bold">Synthesis Route Planner</h2>
          <p className="text-muted-foreground">
            AI-powered retrosynthesis for {compoundName}
          </p>
        </div>
        <div className="flex gap-2">
          {routes.length >= 2 && (
            <Button
              onClick={() => {
                setCompareMode(!compareMode);
                if (!compareMode) {
                  setSelectedRoutes([0, 1]);
                } else {
                  setSelectedRoutes([]);
                }
              }}
              variant={compareMode ? "default" : "outline"}
              size="lg"
            >
              <GitCompare className="mr-2 h-4 w-4" />
              {compareMode ? "Exit Compare" : "Compare Routes"}
            </Button>
          )}
          <Button
            onClick={handleGenerate}
            disabled={generateMutation.isPending}
            size="lg"
          >
            {generateMutation.isPending ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                Generating Routes...
              </>
            ) : (
              <>
                <FlaskConical className="mr-2 h-4 w-4" />
                Generate Synthesis Routes
              </>
            )}
          </Button>
        </div>
      </div>

      {generateMutation.isError && (
        <Alert variant="destructive">
          <AlertCircle className="h-4 w-4" />
          <AlertDescription>
            Failed to generate synthesis routes. Please try again.
          </AlertDescription>
        </Alert>
      )}

      {routes.length > 0 && (
        compareMode ? (
          /* Comparison View */
          <div className="space-y-6">
            <Card>
              <CardHeader>
                <CardTitle>Route Comparison</CardTitle>
                <CardDescription>Side-by-side comparison of synthesis routes</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-2 gap-6">
                  {selectedRoutes.map((routeIndex) => {
                    const route = routes[routeIndex];
                    return (
                      <div key={routeIndex} className="space-y-4">
                        <div className="flex items-center justify-between">
                          <h3 className="text-lg font-semibold">Route {routeIndex + 1}</h3>
                          <Badge className={difficultyColors[route.difficulty]}>
                            {difficultyLabels[route.difficulty]}
                          </Badge>
                        </div>
                        
                        {/* Metrics Comparison */}
                        <div className="space-y-3 p-4 bg-muted/30 rounded-lg">
                          <div className="flex justify-between items-center">
                            <span className="text-sm text-muted-foreground">Feasibility Score</span>
                            <span className="font-bold text-lg">{route.feasibilityScore}/100</span>
                          </div>
                          <Separator />
                          <div className="flex justify-between">
                            <span className="text-sm text-muted-foreground">Total Steps</span>
                            <span className="font-medium">{route.totalSteps}</span>
                          </div>
                          <div className="flex justify-between">
                            <span className="text-sm text-muted-foreground">Total Cost</span>
                            <span className="font-medium">${route.totalCost.toFixed(0)}</span>
                          </div>
                          <div className="flex justify-between">
                            <span className="text-sm text-muted-foreground">Overall Yield</span>
                            <span className="font-medium">{route.overallYield}</span>
                          </div>
                          <div className="flex justify-between">
                            <span className="text-sm text-muted-foreground">Estimated Time</span>
                            <span className="font-medium">{route.estimatedTime}</span>
                          </div>
                        </div>

                        {/* Starting Materials */}
                        <div>
                          <h4 className="font-semibold mb-2">Starting Materials</h4>
                          <div className="space-y-2">
                            {route.startingMaterials.map((material, idx) => (
                              <div key={idx} className="flex justify-between items-center text-sm p-2 bg-background rounded">
                                <span>{material.name}</span>
                                <Badge variant={material.availability === "commercial" ? "default" : "outline"}>
                                  ${material.cost}
                                </Badge>
                              </div>
                            ))}
                          </div>
                        </div>

                        {/* Steps Summary */}
                        <div>
                          <h4 className="font-semibold mb-2">Synthesis Steps</h4>
                          <div className="space-y-2">
                            {route.steps.map((step) => (
                              <div key={step.stepNumber} className="p-3 bg-background rounded border">
                                <div className="flex items-center justify-between mb-1">
                                  <span className="font-medium text-sm">Step {step.stepNumber}: {step.reaction}</span>
                                  <Badge className={difficultyColors[step.difficulty]} variant="outline">
                                    {difficultyLabels[step.difficulty]}
                                  </Badge>
                                </div>
                                <div className="grid grid-cols-2 gap-2 text-xs text-muted-foreground">
                                  <div>Yield: {step.yield}</div>
                                  <div>Cost: ${step.estimatedCost}</div>
                                </div>
                              </div>
                            ))}
                          </div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              </CardContent>
            </Card>
          </div>
        ) : (
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* Route Selection Sidebar */}
            <div className="space-y-4">
              <h3 className="font-semibold text-lg">Available Routes</h3>
              {routes.map((route, index) => (
                <Card
                  key={route.routeId}
                  className={`cursor-pointer transition-all ${
                    selectedRoute === index
                      ? "ring-2 ring-primary"
                      : "hover:shadow-md"
                  }`}
                  onClick={() => setSelectedRoute(index)}
                >
                  <CardHeader className="pb-3">
                    <CardTitle className="text-base">
                      Route {index + 1}
                    </CardTitle>
                    <CardDescription className="flex items-center gap-2">
                      <Badge className={difficultyColors[route.difficulty]}>
                        {difficultyLabels[route.difficulty]}
                      </Badge>
                      <span className="text-sm">
                        Score: {route.feasibilityScore}/100
                      </span>
                    </CardDescription>
                  </CardHeader>
                  <CardContent className="space-y-2">
                    <div className="flex items-center justify-between text-sm">
                      <span className="text-muted-foreground">Steps:</span>
                      <span className="font-medium">{route.totalSteps}</span>
                    </div>
                    <div className="flex items-center justify-between text-sm">
                      <span className="text-muted-foreground">Cost:</span>
                      <span className="font-medium">${route.totalCost.toFixed(0)}</span>
                    </div>
                    <div className="flex items-center justify-between text-sm">
                      <span className="text-muted-foreground">Yield:</span>
                      <span className="font-medium">{route.overallYield}</span>
                    </div>
                    <div className="flex items-center justify-between text-sm">
                      <span className="text-muted-foreground">Time:</span>
                      <span className="font-medium">{route.estimatedTime}</span>
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>

          {/* Route Details */}
          {selectedRoute !== null && routes[selectedRoute] && (
            <div className="lg:col-span-2 space-y-6">
              {/* Route Summary */}
              <Card>
                <CardHeader>
                  <CardTitle>Route {selectedRoute + 1} Details</CardTitle>
                  <CardDescription>
                    {routes[selectedRoute].summary}
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div className="flex items-center gap-2">
                      <TrendingUp className="h-4 w-4 text-muted-foreground" />
                      <div>
                        <p className="text-xs text-muted-foreground">Feasibility</p>
                        <p className="text-lg font-semibold">
                          {routes[selectedRoute].feasibilityScore}/100
                        </p>
                      </div>
                    </div>
                    <div className="flex items-center gap-2">
                      <DollarSign className="h-4 w-4 text-muted-foreground" />
                      <div>
                        <p className="text-xs text-muted-foreground">Total Cost</p>
                        <p className="text-lg font-semibold">
                          ${routes[selectedRoute].totalCost.toFixed(0)}
                        </p>
                      </div>
                    </div>
                    <div className="flex items-center gap-2">
                      <Clock className="h-4 w-4 text-muted-foreground" />
                      <div>
                        <p className="text-xs text-muted-foreground">Time</p>
                        <p className="text-lg font-semibold">
                          {routes[selectedRoute].estimatedTime}
                        </p>
                      </div>
                    </div>
                    <div className="flex items-center gap-2">
                      <FlaskConical className="h-4 w-4 text-muted-foreground" />
                      <div>
                        <p className="text-xs text-muted-foreground">Yield</p>
                        <p className="text-lg font-semibold">
                          {routes[selectedRoute].overallYield}
                        </p>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>

              {/* Starting Materials */}
              <Card>
                <CardHeader>
                  <CardTitle>Starting Materials</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="space-y-3">
                    {routes[selectedRoute].startingMaterials.map((material, idx) => (
                      <div
                        key={idx}
                        className="flex items-center justify-between p-3 border rounded-lg"
                      >
                        <div>
                          <p className="font-medium">{material.name}</p>
                          <p className="text-xs text-muted-foreground font-mono">
                            {material.smiles}
                          </p>
                        </div>
                        <div className="flex items-center gap-3">
                          <Badge
                            variant={
                              material.availability === "commercial"
                                ? "default"
                                : "secondary"
                            }
                          >
                            {material.availability === "commercial"
                              ? "Commercial"
                              : "Synthesize"}
                          </Badge>
                          <span className="text-sm font-medium">
                            ${material.cost.toFixed(0)}
                          </span>
                        </div>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>

              {/* Synthesis Steps */}
              <Card>
                <CardHeader>
                  <CardTitle>Synthesis Steps</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="space-y-6">
                    {routes[selectedRoute].steps.map((step, idx) => (
                      <div key={step.stepNumber}>
                        <div className="flex items-start gap-4">
                          <div className="flex-shrink-0 w-8 h-8 rounded-full bg-primary text-primary-foreground flex items-center justify-center font-semibold">
                            {step.stepNumber}
                          </div>
                          <div className="flex-1 space-y-2">
                            <div className="flex items-center justify-between">
                              <h4 className="font-semibold">{step.reaction}</h4>
                              <Badge className={difficultyColors[step.difficulty]}>
                                {difficultyLabels[step.difficulty]}
                              </Badge>
                            </div>
                            <div className="grid grid-cols-2 gap-4 text-sm">
                              <div>
                                <p className="text-muted-foreground">Reagents:</p>
                                <p className="font-medium">
                                  {step.reagents.join(", ")}
                                </p>
                              </div>
                              <div>
                                <p className="text-muted-foreground">Conditions:</p>
                                <p className="font-medium">{step.conditions}</p>
                              </div>
                              <div>
                                <p className="text-muted-foreground">Expected Yield:</p>
                                <p className="font-medium">{step.yield}</p>
                              </div>
                              <div>
                                <p className="text-muted-foreground">Estimated Cost:</p>
                                <p className="font-medium">
                                  ${step.estimatedCost.toFixed(0)}
                                </p>
                              </div>
                            </div>
                            {step.notes && (
                              <Alert>
                                <AlertDescription className="text-sm">
                                  {step.notes}
                                </AlertDescription>
                              </Alert>
                            )}
                          </div>
                        </div>
                        {idx < routes[selectedRoute].steps.length - 1 && (
                          <div className="flex justify-center my-4">
                            <ChevronRight className="h-6 w-6 text-muted-foreground rotate-90" />
                          </div>
                        )}
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            </div>
          )}
          </div>
        )
      )}
    </div>
  );
}
