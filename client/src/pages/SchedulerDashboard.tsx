import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Checkbox } from "@/components/ui/checkbox";
import { trpc } from "@/lib/trpc";
import { Activity, Calendar, CheckCircle2, Clock, PlayCircle, RefreshCw, AlertCircle, Target, Plus, X, Sparkles, Filter } from "lucide-react";
import { toast } from "sonner";
import { useState } from "react";

export default function SchedulerDashboard() {
  const { data: schedulerStatus, isLoading, refetch } = trpc.scheduler.status.useQuery();
  const { data: researchGoals, refetch: refetchGoals } = trpc.scheduler.getResearchGoals.useQuery();
  const { data: medicalTrends } = trpc.scheduler.getMedicalTrends.useQuery();
  
  const [newGoal, setNewGoal] = useState("");
  const [isAddingGoal, setIsAddingGoal] = useState(false);
  
  const therapeuticAreas = [
    "CNS Disorders",
    "Oncology",
    "Cardiovascular",
    "Metabolic Diseases",
    "Infectious Diseases",
    "Immunology",
    "Rare Diseases",
    "Respiratory",
    "Pain Management",
    "Psychedelic Therapy",
  ];
  
  const [selectedAreas, setSelectedAreas] = useState<string[]>(therapeuticAreas); // All selected by default
  
  const runNow = trpc.scheduler.runNow.useMutation({
    onSuccess: () => {
      toast.success("Scheduler triggered successfully");
      refetch();
    },
    onError: (error) => {
      toast.error(`Failed to run scheduler: ${error.message}`);
    },
  });
  
  const saveGoals = trpc.scheduler.saveResearchGoals.useMutation({
    onSuccess: () => {
      toast.success("Research goals saved successfully");
      refetchGoals();
    },
    onError: (error) => {
      toast.error(`Failed to save goals: ${error.message}`);
    },
  });
  
  const refreshTrends = trpc.scheduler.refreshMedicalTrends.useMutation({
    onSuccess: () => {
      toast.success("Medical trends refreshed");
      refetch();
    },
    onError: (error) => {
      toast.error(`Failed to refresh trends: ${error.message}`);
    },
  });

  if (isLoading) {
    return (
      <div className="container mx-auto py-8">
        <div className="flex items-center justify-center h-64">
          <RefreshCw className="h-8 w-8 animate-spin text-muted-foreground" />
        </div>
      </div>
    );
  }

  const getStatusColor = (running: boolean) => {
    return running
      ? "bg-green-100 text-green-800 border-green-200"
      : "bg-blue-100 text-blue-800 border-blue-200";
  };

  const formatDate = (dateString: string | null) => {
    if (!dateString) return "Never";
    return new Date(dateString).toLocaleString();
  };

  return (
    <div className="container mx-auto py-8 space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold">Scheduler Dashboard</h1>
          <p className="text-muted-foreground mt-1">
            Monitor autonomous research engine and scheduled tasks
          </p>
        </div>
        <Button
          onClick={() => runNow.mutate()}
          disabled={runNow.isPending || schedulerStatus?.running}
        >
          <PlayCircle className="h-4 w-4 mr-2" />
          Run Now
        </Button>
      </div>

      {/* Status Overview */}
      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Status</CardTitle>
            <Activity className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <Badge className={getStatusColor(schedulerStatus?.running || false)}>
              {schedulerStatus?.running ? "RUNNING" : "IDLE"}
            </Badge>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Schedule</CardTitle>
            <Clock className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-sm font-medium">
              {schedulerStatus?.config?.cronSchedule || "Not configured"}
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Next Run</CardTitle>
            <Calendar className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-sm font-medium">
              {formatDate(schedulerStatus?.nextRun || null)}
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Enabled</CardTitle>
            <CheckCircle2 className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <Badge variant={schedulerStatus?.config?.enabled ? "default" : "outline"}>
              {schedulerStatus?.config?.enabled ? "Yes" : "No"}
            </Badge>
          </CardContent>
        </Card>
      </div>

      {/* Scheduler Configuration */}
      <Card>
        <CardHeader>
          <CardTitle>Scheduler Configuration</CardTitle>
          <CardDescription>
            Autonomous research engine configuration and status
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            <div className="flex items-center justify-between p-4 border rounded-lg">
              <div className="flex items-center gap-4">
                <div className="flex h-10 w-10 items-center justify-center rounded-full bg-primary/10">
                  <Activity className="h-5 w-5 text-primary" />
                </div>
                <div>
                  <h4 className="font-semibold">Autonomous Discovery Engine</h4>
                  <p className="text-sm text-muted-foreground">
                    Imports analog discoveries from research engine
                  </p>
                </div>
              </div>
              <div className="text-right">
                <Badge className={getStatusColor(schedulerStatus?.running || false)}>
                  {schedulerStatus?.running ? "RUNNING" : "IDLE"}
                </Badge>
                <p className="text-xs text-muted-foreground mt-1">
                  Next: {formatDate(schedulerStatus?.nextRun || null)}
                </p>
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Research Goals */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle className="flex items-center gap-2">
                <Target className="h-5 w-5" />
                Research Goals
              </CardTitle>
              <CardDescription>
                Define custom research targets for the autonomous engine
              </CardDescription>
            </div>
            <Button
              variant="outline"
              size="sm"
              onClick={() => setIsAddingGoal(!isAddingGoal)}
            >
              <Plus className="h-4 w-4 mr-2" />
              Add Goal
            </Button>
          </div>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            {isAddingGoal && (
              <div className="flex gap-2">
                <Input
                  placeholder="e.g., psychedelics, nootropics, anxiolytics"
                  value={newGoal}
                  onChange={(e) => setNewGoal(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === 'Enter' && newGoal.trim()) {
                      const currentGoals = researchGoals?.goals || [];
                      saveGoals.mutate({ goals: [...currentGoals, newGoal.trim()] });
                      setNewGoal("");
                      setIsAddingGoal(false);
                    }
                  }}
                />
                <Button
                  size="sm"
                  onClick={() => {
                    if (newGoal.trim()) {
                      const currentGoals = researchGoals?.goals || [];
                      saveGoals.mutate({ goals: [...currentGoals, newGoal.trim()] });
                      setNewGoal("");
                      setIsAddingGoal(false);
                    }
                  }}
                >
                  Add
                </Button>
              </div>
            )}
            
            <div className="flex flex-wrap gap-2">
              {researchGoals?.goals?.map((goal: string, index: number) => (
                <Badge key={index} variant="secondary" className="gap-2">
                  {goal}
                  <button
                    onClick={() => {
                      const updatedGoals = researchGoals.goals.filter((_: string, i: number) => i !== index);
                      saveGoals.mutate({ goals: updatedGoals });
                    }}
                    className="hover:text-destructive"
                  >
                    <X className="h-3 w-3" />
                  </button>
                </Badge>
              ))}
              {(!researchGoals?.goals || researchGoals.goals.length === 0) && (
                <p className="text-sm text-muted-foreground">No research goals defined yet</p>
              )}
            </div>
          </div>
        </CardContent>
      </Card>

      {/* AI-Powered Medical Trends */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle className="flex items-center gap-2">
                <Sparkles className="h-5 w-5 text-primary" />
                AI-Discovered Medical Trends
              </CardTitle>
              <CardDescription>
                High-value research targets identified from recent breakthroughs
              </CardDescription>
            </div>
            <Button
              variant="outline"
              size="sm"
              onClick={() => refreshTrends.mutate()}
              disabled={refreshTrends.isPending}
            >
              <RefreshCw className={`h-4 w-4 mr-2 ${refreshTrends.isPending ? 'animate-spin' : ''}`} />
              Refresh
            </Button>
          </div>
        </CardHeader>
        <CardContent>
          <div className="space-y-3">
            {medicalTrends?.trends?.map((trend: any, index: number) => (
              <div key={index} className="p-3 border rounded-lg space-y-2">
                <div className="flex items-start justify-between">
                  <div className="flex-1">
                    <h4 className="font-semibold text-sm">{trend.title}</h4>
                    <p className="text-xs text-muted-foreground mt-1">{trend.description}</p>
                  </div>
                  <Badge variant="outline" className="ml-2">
                    {trend.priority}
                  </Badge>
                </div>
                <div className="flex gap-2 flex-wrap">
                  {trend.keywords?.map((keyword: string, i: number) => (
                    <Badge key={i} variant="secondary" className="text-xs">
                      {keyword}
                    </Badge>
                  ))}
                </div>
              </div>
            ))}
            {(!medicalTrends?.trends || medicalTrends.trends.length === 0) && (
              <p className="text-sm text-muted-foreground">No trends discovered yet. Click Refresh to analyze recent breakthroughs.</p>
            )}
          </div>
        </CardContent>
      </Card>

      {/* Therapeutic Area Filtering */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle className="flex items-center gap-2">
                <Filter className="h-5 w-5" />
                Therapeutic Area Filtering
              </CardTitle>
              <CardDescription>
                Select focus areas for autonomous discovery engine
              </CardDescription>
            </div>
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                onClick={() => setSelectedAreas(therapeuticAreas)}
              >
                Select All
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={() => setSelectedAreas([])}
              >
                Clear All
              </Button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-3 gap-4">
            {therapeuticAreas.map((area) => (
              <div key={area} className="flex items-center space-x-2">
                <Checkbox
                  id={area}
                  checked={selectedAreas.includes(area)}
                  onCheckedChange={(checked) => {
                    if (checked) {
                      setSelectedAreas([...selectedAreas, area]);
                    } else {
                      setSelectedAreas(selectedAreas.filter(a => a !== area));
                    }
                  }}
                />
                <label
                  htmlFor={area}
                  className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70 cursor-pointer"
                >
                  {area}
                </label>
              </div>
            ))}
          </div>
          <div className="mt-4 p-3 bg-muted rounded-lg">
            <p className="text-xs text-muted-foreground">
              <strong>Note:</strong> Therapeutic area filtering influences the discovery engine's focus but doesn't guarantee exclusivity. 
              The engine will prioritize selected areas while maintaining diversity across the pharmaceutical landscape.
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Configuration Details */}
      <Card>
        <CardHeader>
          <CardTitle>Configuration Details</CardTitle>
          <CardDescription>Scheduler settings and parameters</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid gap-4 md:grid-cols-2">
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Cron Schedule</p>
              <p className="text-sm font-mono">{schedulerStatus?.config?.cronSchedule || "N/A"}</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Timezone</p>
              <p className="text-sm">America/New_York (EST)</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Enabled</p>
              <p className="text-sm">{schedulerStatus?.config?.enabled ? "Yes" : "No"}</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Next Execution</p>
              <p className="text-sm">{formatDate(schedulerStatus?.nextRun || null)}</p>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
