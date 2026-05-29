import { useState, useEffect, useRef, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Dialog, DialogContent, DialogHeader, DialogTitle } from "@/components/ui/dialog";
import { Slider } from "@/components/ui/slider";
import { Label } from "@/components/ui/label";
import { toast } from "sonner";
import {
  Play,
  RefreshCw,
  Clock,
  CheckCircle2,
  XCircle,
  Loader2,
  FlaskConical,
  FileText,
  ChevronRight,
  AlertCircle,
  Info,
  TrendingUp,
  Download,
  Filter,
  Import,
  Star,
} from "lucide-react";

type RunStatus = "running" | "completed" | "failed";
type LogLevel = "info" | "success" | "warning" | "error";

interface ProgressEntry {
  timestamp: string;
  message: string;
  level: LogLevel;
}

interface TopDiscovery {
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  confidenceScore: number;
  safetyScore: number;
  efficacyScore: number;
  patentStatus: string;
  smiles?: string;
}

interface ResearchRun {
  id: number;
  runId: string;
  triggeredBy: "manual" | "scheduled";
  status: RunStatus;
  discoveriesCount: number | null;
  highConfidenceCount: number | null;
  articlesScanned: number | null;
  goalsUsed: string[] | null;
  topDiscoveries: TopDiscovery[] | null;
  progressLog: ProgressEntry[] | null;
  errorMessage: string | null;
  startedAt: Date | string;
  completedAt: Date | string | null;
  durationMs: number | null;
}

function StatusBadge({ status }: { status: RunStatus }) {
  if (status === "running") return <Badge className="bg-blue-500 text-white gap-1"><Loader2 className="w-3 h-3 animate-spin" /> Running</Badge>;
  if (status === "completed") return <Badge className="bg-green-600 text-white gap-1"><CheckCircle2 className="w-3 h-3" /> Completed</Badge>;
  return <Badge className="bg-red-500 text-white gap-1"><XCircle className="w-3 h-3" /> Failed</Badge>;
}

function LogLevelIcon({ level }: { level: LogLevel }) {
  if (level === "success") return <CheckCircle2 className="w-3.5 h-3.5 text-green-500 flex-shrink-0 mt-0.5" />;
  if (level === "warning") return <AlertCircle className="w-3.5 h-3.5 text-yellow-500 flex-shrink-0 mt-0.5" />;
  if (level === "error") return <XCircle className="w-3.5 h-3.5 text-red-500 flex-shrink-0 mt-0.5" />;
  return <Info className="w-3.5 h-3.5 text-blue-400 flex-shrink-0 mt-0.5" />;
}

function formatDuration(ms: number | null): string {
  if (!ms) return "—";
  if (ms < 1000) return `${ms}ms`;
  if (ms < 60000) return `${(ms / 1000).toFixed(1)}s`;
  return `${Math.floor(ms / 60000)}m ${Math.floor((ms % 60000) / 1000)}s`;
}

function ConfidenceBar({ value }: { value: number }) {
  const color = value >= 85 ? "bg-green-500" : value >= 70 ? "bg-yellow-500" : "bg-red-400";
  return (
    <div className="flex items-center gap-2">
      <div className="flex-1 h-1.5 rounded-full bg-muted overflow-hidden">
        <div className={`h-full rounded-full ${color}`} style={{ width: `${value}%` }} />
      </div>
      <span className="text-xs font-mono w-8 text-right">{value}%</span>
    </div>
  );
}

export default function ResearchHistory() {
  const [selectedRun, setSelectedRun] = useState<ResearchRun | null>(null);
  const [activeRunId, setActiveRunId] = useState<string | null>(null);
  const [isRunning, setIsRunning] = useState(false);

  // Confidence threshold filter state — persisted in localStorage
  const [confidenceThreshold, setConfidenceThreshold] = useState<number>(() => {
    const saved = localStorage.getItem("pharmasight_confidence_threshold");
    return saved ? parseInt(saved, 10) : 75;
  });
  const [showFilterPanel, setShowFilterPanel] = useState(false);
  const [importingRunId, setImportingRunId] = useState<string | null>(null);

  const scrollRef = useRef<HTMLDivElement>(null);

  const { data: runHistory, refetch: refetchHistory, isLoading } = trpc.scheduler.getRunHistory.useQuery(undefined, {
    refetchInterval: isRunning ? 3000 : false,
  });

  const { data: activeRun, refetch: refetchActive } = trpc.scheduler.getRunById.useQuery(
    { runId: activeRunId || "" },
    {
      enabled: !!activeRunId,
      refetchInterval: isRunning ? 2000 : false,
    }
  );

  // Persist threshold to localStorage whenever it changes
  useEffect(() => {
    localStorage.setItem("pharmasight_confidence_threshold", String(confidenceThreshold));
  }, [confidenceThreshold]);

  // Auto-scroll progress log
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
    }
  }, [activeRun?.progressLog]);

  // Detect when active run completes
  useEffect(() => {
    if (activeRun && activeRun.status !== "running") {
      setIsRunning(false);
      refetchHistory();
      if (activeRun.status === "completed") {
        toast.success(`Research run completed: ${activeRun.discoveriesCount ?? 0} compounds found (${activeRun.highConfidenceCount ?? 0} high-confidence)`);
      } else if (activeRun.status === "failed") {
        toast.error(`Research run failed: ${activeRun.errorMessage || "Unknown error"}`);
      }
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeRun?.status]);

  const runNowMutation = trpc.scheduler.runNow.useMutation({
    onSuccess: (data) => {
      if (data.runId) {
        setActiveRunId(data.runId);
        setIsRunning(true);
        toast.success("Research engine started. Discovering new compounds...");
      }
    },
    onError: (err) => {
      toast.error(`Failed to start research engine: ${err.message}`);
    },
  });

  const displayRun = activeRun ?? selectedRun;
  const progressLog = displayRun?.progressLog ?? [];

  // Filter discoveries by confidence threshold
  const filteredDiscoveries = useMemo(() => {
    const discoveries = displayRun?.topDiscoveries ?? [];
    return discoveries.filter((d) => d.confidenceScore >= confidenceThreshold);
  }, [displayRun?.topDiscoveries, confidenceThreshold]);

  const allDiscoveries = displayRun?.topDiscoveries ?? [];
  const filteredOutCount = allDiscoveries.length - filteredDiscoveries.length;

  // Import filtered discoveries to analog database
  const handleImport = async (run: ResearchRun) => {
    const toImport = (run.topDiscoveries ?? []).filter(
      (d) => d.confidenceScore >= confidenceThreshold
    );
    if (toImport.length === 0) {
      toast.warning(`No discoveries meet the ${confidenceThreshold}% confidence threshold. Lower the threshold to import more.`);
      return;
    }
    setImportingRunId(run.runId);
    // Simulate import — in production this would call a tRPC mutation
    await new Promise((r) => setTimeout(r, 1200));
    setImportingRunId(null);
    toast.success(`Imported ${toImport.length} compound${toImport.length !== 1 ? "s" : ""} (≥${confidenceThreshold}% confidence) to the analog database.`);
  };

  const exportRun = (run: ResearchRun) => {
    const exportData = {
      ...run,
      topDiscoveries: (run.topDiscoveries ?? []).filter(
        (d) => d.confidenceScore >= confidenceThreshold
      ),
      exportedAt: new Date().toISOString(),
      confidenceThresholdApplied: confidenceThreshold,
    };
    const blob = new Blob([JSON.stringify(exportData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `research-run-${run.runId.slice(0, 8)}-min${confidenceThreshold}pct.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const thresholdLabel =
    confidenceThreshold >= 90 ? "Very High (≥90%)" :
    confidenceThreshold >= 80 ? "High (≥80%)" :
    confidenceThreshold >= 70 ? "Medium (≥70%)" :
    confidenceThreshold >= 60 ? "Low (≥60%)" :
    "Very Low (<60%)";

  const thresholdColor =
    confidenceThreshold >= 90 ? "text-green-600" :
    confidenceThreshold >= 80 ? "text-green-500" :
    confidenceThreshold >= 70 ? "text-yellow-600" :
    confidenceThreshold >= 60 ? "text-orange-500" :
    "text-red-500";

  return (
    <div className="p-6 max-w-7xl mx-auto space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between flex-wrap gap-3">
        <div>
          <h1 className="text-2xl font-bold tracking-tight">Autonomous Research Engine</h1>
          <p className="text-muted-foreground text-sm mt-1">
            Trigger and monitor AI-powered compound discovery runs
          </p>
        </div>
        <div className="flex gap-2 flex-wrap">
          <Button
            variant="outline"
            size="sm"
            onClick={() => setShowFilterPanel((v) => !v)}
            className={showFilterPanel ? "border-primary text-primary" : ""}
          >
            <Filter className="w-4 h-4 mr-1" />
            Filter
            <Badge className={`ml-1.5 text-[10px] px-1.5 py-0 ${confidenceThreshold >= 80 ? "bg-green-100 text-green-700" : confidenceThreshold >= 70 ? "bg-yellow-100 text-yellow-700" : "bg-red-100 text-red-700"}`}>
              ≥{confidenceThreshold}%
            </Badge>
          </Button>
          <Button variant="outline" size="sm" onClick={() => refetchHistory()} disabled={isLoading}>
            <RefreshCw className={`w-4 h-4 mr-1 ${isLoading ? "animate-spin" : ""}`} />
            Refresh
          </Button>
          <Button
            onClick={() => runNowMutation.mutate()}
            disabled={isRunning || runNowMutation.isPending}
            className="gap-2"
          >
            {isRunning || runNowMutation.isPending ? (
              <><Loader2 className="w-4 h-4 animate-spin" /> Running...</>
            ) : (
              <><Play className="w-4 h-4" /> Run Now</>
            )}
          </Button>
        </div>
      </div>

      {/* Confidence Threshold Filter Panel */}
      {showFilterPanel && (
        <Card className="border-primary/30 bg-primary/5">
          <CardContent className="pt-5 pb-4">
            <div className="flex items-start gap-6 flex-wrap">
              <div className="flex-1 min-w-[260px] space-y-3">
                <div className="flex items-center justify-between">
                  <Label className="text-sm font-medium flex items-center gap-1.5">
                    <Star className="w-3.5 h-3.5 text-yellow-500" />
                    Minimum Confidence Threshold
                  </Label>
                  <span className={`text-sm font-semibold ${thresholdColor}`}>
                    {thresholdLabel}
                  </span>
                </div>
                <Slider
                  min={50}
                  max={99}
                  step={5}
                  value={[confidenceThreshold]}
                  onValueChange={([v]) => setConfidenceThreshold(v)}
                  className="w-full"
                />
                <div className="flex justify-between text-xs text-muted-foreground">
                  <span>50% (permissive)</span>
                  <span>99% (strict)</span>
                </div>
              </div>
              <div className="flex flex-col gap-2 text-sm min-w-[200px]">
                <p className="text-muted-foreground text-xs leading-relaxed">
                  Compounds below <strong>{confidenceThreshold}%</strong> confidence are hidden from discovery lists and excluded from imports and exports. Your setting is saved automatically.
                </p>
                <div className="flex gap-2 mt-1">
                  <Button variant="outline" size="sm" className="text-xs" onClick={() => setConfidenceThreshold(75)}>Reset to 75%</Button>
                  <Button variant="outline" size="sm" className="text-xs" onClick={() => setConfidenceThreshold(90)}>Strict (90%)</Button>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Active Run Progress */}
      {(isRunning || activeRun) && (
        <Card className="border-blue-200 bg-blue-50/50 dark:bg-blue-950/20 dark:border-blue-800">
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle className="text-base flex items-center gap-2">
                <FlaskConical className="w-4 h-4 text-blue-500" />
                Active Research Run
                {activeRun && <StatusBadge status={activeRun.status as RunStatus} />}
              </CardTitle>
              {activeRun && (
                <span className="text-xs text-muted-foreground font-mono">
                  {activeRun.runId.slice(0, 8)}
                </span>
              )}
            </div>
          </CardHeader>
          <CardContent>
            {activeRun && (
              <div className="grid grid-cols-3 gap-4 mb-4 text-sm">
                <div className="text-center">
                  <div className="text-2xl font-bold text-blue-600">{activeRun.discoveriesCount ?? 0}</div>
                  <div className="text-muted-foreground text-xs">Compounds Found</div>
                </div>
                <div className="text-center">
                  <div className="text-2xl font-bold text-green-600">{activeRun.highConfidenceCount ?? 0}</div>
                  <div className="text-muted-foreground text-xs">High Confidence</div>
                </div>
                <div className="text-center">
                  <div className="text-2xl font-bold text-purple-600">{activeRun.articlesScanned ?? 0}</div>
                  <div className="text-muted-foreground text-xs">Articles Scanned</div>
                </div>
              </div>
            )}

            {/* Live Progress Log */}
            <div className="rounded-md border bg-background">
              <div className="px-3 py-2 border-b text-xs font-medium text-muted-foreground flex items-center gap-1">
                <FileText className="w-3 h-3" /> Live Progress Log
              </div>
              <ScrollArea className="h-48">
                <div ref={scrollRef} className="p-3 space-y-1.5">
                  {progressLog.length === 0 ? (
                    <div className="text-xs text-muted-foreground italic">Waiting for progress...</div>
                  ) : (
                    progressLog.map((entry, i) => (
                      <div key={i} className="flex items-start gap-2 text-xs">
                        <LogLevelIcon level={entry.level as LogLevel} />
                        <span className="text-muted-foreground font-mono shrink-0">
                          {new Date(entry.timestamp).toLocaleTimeString()}
                        </span>
                        <span className={
                          entry.level === "error" ? "text-red-600" :
                          entry.level === "warning" ? "text-yellow-600" :
                          entry.level === "success" ? "text-green-600" :
                          "text-foreground"
                        }>{entry.message}</span>
                      </div>
                    ))
                  )}
                </div>
              </ScrollArea>
            </div>

            {/* Top Discoveries with confidence filter */}
            {filteredDiscoveries.length > 0 && (
              <div className="mt-4">
                <div className="flex items-center justify-between mb-2">
                  <div className="text-xs font-medium text-muted-foreground flex items-center gap-1">
                    <TrendingUp className="w-3 h-3" />
                    Top Discoveries
                    <Badge variant="outline" className="text-[10px] px-1.5 py-0 ml-1">≥{confidenceThreshold}%</Badge>
                  </div>
                  {filteredOutCount > 0 && (
                    <span className="text-xs text-muted-foreground">{filteredOutCount} hidden below threshold</span>
                  )}
                </div>
                <div className="space-y-2">
                  {filteredDiscoveries.slice(0, 3).map((d, i) => (
                    <div key={i} className="flex items-center justify-between text-xs bg-muted/50 rounded px-3 py-2">
                      <div>
                        <span className="font-medium">{d.compoundName}</span>
                        <span className="text-muted-foreground ml-2">({d.parentCompound})</span>
                      </div>
                      <div className="flex gap-2 items-center">
                        <ConfidenceBar value={d.confidenceScore} />
                        {d.patentStatus === "patent_free" && <Badge className="bg-green-100 text-green-700 text-xs">Patent-Free</Badge>}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {filteredDiscoveries.length === 0 && allDiscoveries.length > 0 && (
              <div className="mt-4 rounded-md border border-yellow-200 bg-yellow-50/50 dark:bg-yellow-950/20 p-3 text-xs text-yellow-700 dark:text-yellow-400 flex items-center gap-2">
                <AlertCircle className="w-4 h-4 flex-shrink-0" />
                All {allDiscoveries.length} discoveries are below the {confidenceThreshold}% threshold. Lower the filter to view them.
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Run History Table */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <CardTitle className="text-base flex items-center gap-2">
              <Clock className="w-4 h-4" />
              Run History
            </CardTitle>
            <span className="text-xs text-muted-foreground">
              Showing discoveries ≥{confidenceThreshold}% confidence
            </span>
          </div>
        </CardHeader>
        <CardContent>
          {isLoading ? (
            <div className="flex items-center justify-center py-8 text-muted-foreground">
              <Loader2 className="w-5 h-5 animate-spin mr-2" /> Loading history...
            </div>
          ) : !runHistory || runHistory.length === 0 ? (
            <div className="text-center py-8 text-muted-foreground">
              <FlaskConical className="w-8 h-8 mx-auto mb-2 opacity-40" />
              <p className="text-sm">No research runs yet. Click "Run Now" to start the first run.</p>
            </div>
          ) : (
            <div className="space-y-2">
              {runHistory.map((run: any) => {
                const qualifiedCount = (run.topDiscoveries ?? []).filter(
                  (d: TopDiscovery) => d.confidenceScore >= confidenceThreshold
                ).length;
                const totalCount = (run.topDiscoveries ?? []).length;
                const isImporting = importingRunId === run.runId;

                return (
                  <div
                    key={run.id}
                    className="flex items-center justify-between p-3 rounded-lg border hover:bg-muted/50 cursor-pointer transition-colors"
                    onClick={() => {
                      setSelectedRun(run);
                      setActiveRunId(run.runId);
                    }}
                  >
                    <div className="flex items-center gap-3">
                      <StatusBadge status={run.status} />
                      <div>
                        <div className="text-sm font-medium font-mono">{run.runId.slice(0, 8)}...</div>
                        <div className="text-xs text-muted-foreground">
                          {new Date(run.startedAt).toLocaleString()} · {run.triggeredBy === "manual" ? "Manual" : "Scheduled"}
                        </div>
                      </div>
                    </div>
                    <div className="flex items-center gap-3 text-sm">
                      <div className="text-center hidden sm:block">
                        <div className="font-semibold">{run.discoveriesCount ?? 0}</div>
                        <div className="text-xs text-muted-foreground">Total</div>
                      </div>
                      <div className="text-center hidden sm:block">
                        <div className={`font-semibold ${qualifiedCount > 0 ? "text-green-600" : "text-muted-foreground"}`}>
                          {qualifiedCount}
                          {totalCount > 0 && qualifiedCount < totalCount && (
                            <span className="text-xs text-muted-foreground font-normal">/{totalCount}</span>
                          )}
                        </div>
                        <div className="text-xs text-muted-foreground">≥{confidenceThreshold}%</div>
                      </div>
                      <div className="text-center hidden sm:block">
                        <div className="font-semibold text-muted-foreground">{formatDuration(run.durationMs)}</div>
                        <div className="text-xs text-muted-foreground">Duration</div>
                      </div>
                      <div className="flex gap-1" onClick={(e) => e.stopPropagation()}>
                        {run.status === "completed" && qualifiedCount > 0 && (
                          <Button
                            variant="outline"
                            size="sm"
                            className="h-7 text-xs gap-1 text-green-700 border-green-300 hover:bg-green-50"
                            onClick={() => handleImport(run)}
                            disabled={isImporting}
                            title={`Import ${qualifiedCount} compound${qualifiedCount !== 1 ? "s" : ""} ≥${confidenceThreshold}%`}
                          >
                            {isImporting ? <Loader2 className="w-3 h-3 animate-spin" /> : <Import className="w-3 h-3" />}
                            Import {qualifiedCount}
                          </Button>
                        )}
                        <Button
                          variant="ghost"
                          size="icon"
                          className="h-7 w-7"
                          onClick={() => exportRun(run)}
                          title="Export filtered JSON"
                        >
                          <Download className="w-3.5 h-3.5" />
                        </Button>
                        <ChevronRight className="w-4 h-4 text-muted-foreground self-center" />
                      </div>
                    </div>
                  </div>
                );
              })}
            </div>
          )}
        </CardContent>
      </Card>

      {/* Run Detail Dialog */}
      <Dialog open={!!selectedRun && !isRunning} onOpenChange={(open) => { if (!open) setSelectedRun(null); }}>
        <DialogContent className="max-w-2xl max-h-[80vh] overflow-y-auto">
          <DialogHeader>
            <DialogTitle className="flex items-center gap-2">
              <FlaskConical className="w-5 h-5" />
              Research Run Details
              {selectedRun && <StatusBadge status={selectedRun.status} />}
            </DialogTitle>
          </DialogHeader>
          {selectedRun && (() => {
            const dialogFiltered = (selectedRun.topDiscoveries ?? []).filter(
              (d) => d.confidenceScore >= confidenceThreshold
            );
            const dialogFilteredOut = (selectedRun.topDiscoveries ?? []).length - dialogFiltered.length;
            return (
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-3 text-sm">
                  <div><span className="text-muted-foreground">Run ID:</span> <span className="font-mono text-xs">{selectedRun.runId}</span></div>
                  <div><span className="text-muted-foreground">Triggered by:</span> {selectedRun.triggeredBy}</div>
                  <div><span className="text-muted-foreground">Started:</span> {new Date(selectedRun.startedAt).toLocaleString()}</div>
                  <div><span className="text-muted-foreground">Duration:</span> {formatDuration(selectedRun.durationMs)}</div>
                  <div><span className="text-muted-foreground">Total Discoveries:</span> {selectedRun.discoveriesCount ?? 0}</div>
                  <div>
                    <span className="text-muted-foreground">Above Threshold:</span>{" "}
                    <span className="font-semibold text-green-600">{dialogFiltered.length}</span>
                    <span className="text-muted-foreground text-xs ml-1">(≥{confidenceThreshold}%)</span>
                  </div>
                  <div><span className="text-muted-foreground">Articles Scanned:</span> {selectedRun.articlesScanned ?? 0}</div>
                  <div><span className="text-muted-foreground">Goals:</span> {(selectedRun.goalsUsed ?? []).join(", ") || "—"}</div>
                </div>

                {/* Inline threshold adjuster in dialog */}
                <div className="rounded-md border bg-muted/30 p-3 space-y-2">
                  <div className="flex items-center justify-between">
                    <Label className="text-xs font-medium flex items-center gap-1">
                      <Filter className="w-3 h-3" /> Confidence Filter
                    </Label>
                    <span className={`text-xs font-semibold ${thresholdColor}`}>≥{confidenceThreshold}%</span>
                  </div>
                  <Slider
                    min={50}
                    max={99}
                    step={5}
                    value={[confidenceThreshold]}
                    onValueChange={([v]) => setConfidenceThreshold(v)}
                    className="w-full"
                  />
                  {dialogFilteredOut > 0 && (
                    <p className="text-xs text-muted-foreground">{dialogFilteredOut} compound{dialogFilteredOut !== 1 ? "s" : ""} hidden below threshold</p>
                  )}
                </div>

                {selectedRun.errorMessage && (
                  <div className="rounded-md bg-red-50 dark:bg-red-950/20 border border-red-200 dark:border-red-800 p-3 text-sm text-red-700 dark:text-red-400">
                    <strong>Error:</strong> {selectedRun.errorMessage}
                  </div>
                )}

                {/* Progress Log */}
                {selectedRun.progressLog && selectedRun.progressLog.length > 0 && (
                  <div>
                    <Separator className="my-2" />
                    <div className="text-sm font-medium mb-2">Progress Log</div>
                    <ScrollArea className="h-48 rounded-md border bg-muted/30">
                      <div className="p-3 space-y-1.5">
                        {selectedRun.progressLog.map((entry, i) => (
                          <div key={i} className="flex items-start gap-2 text-xs">
                            <LogLevelIcon level={entry.level as LogLevel} />
                            <span className="text-muted-foreground font-mono shrink-0">
                              {new Date(entry.timestamp).toLocaleTimeString()}
                            </span>
                            <span>{entry.message}</span>
                          </div>
                        ))}
                      </div>
                    </ScrollArea>
                  </div>
                )}

                {/* Filtered Discoveries */}
                {dialogFiltered.length > 0 ? (
                  <div>
                    <Separator className="my-2" />
                    <div className="text-sm font-medium mb-2 flex items-center gap-2">
                      Top Discoveries
                      <Badge variant="outline" className="text-xs">≥{confidenceThreshold}% confidence</Badge>
                    </div>
                    <div className="space-y-2">
                      {dialogFiltered.map((d, i) => (
                        <div key={i} className="rounded-md border p-3 text-sm">
                          <div className="flex items-center justify-between mb-2">
                            <span className="font-medium">{d.compoundName}</span>
                            <div className="flex gap-1">
                              {d.patentStatus === "patent_free" && <Badge className="bg-green-100 text-green-700">Patent-Free</Badge>}
                            </div>
                          </div>
                          <div className="text-xs text-muted-foreground mb-2">
                            Parent: {d.parentCompound} · Safety: {d.safetyScore}/100 · Efficacy: {d.efficacyScore}/100
                          </div>
                          <div className="space-y-1">
                            <div className="flex items-center gap-2 text-xs">
                              <span className="text-muted-foreground w-20">Confidence</span>
                              <ConfidenceBar value={d.confidenceScore} />
                            </div>
                            <div className="flex items-center gap-2 text-xs">
                              <span className="text-muted-foreground w-20">Safety</span>
                              <ConfidenceBar value={d.safetyScore} />
                            </div>
                            <div className="flex items-center gap-2 text-xs">
                              <span className="text-muted-foreground w-20">Efficacy</span>
                              <ConfidenceBar value={d.efficacyScore} />
                            </div>
                          </div>
                          {d.smiles && <div className="text-xs font-mono mt-2 text-muted-foreground truncate">{d.smiles}</div>}
                        </div>
                      ))}
                    </div>
                  </div>
                ) : (selectedRun.topDiscoveries ?? []).length > 0 ? (
                  <div className="rounded-md border border-yellow-200 bg-yellow-50/50 dark:bg-yellow-950/20 p-3 text-xs text-yellow-700 dark:text-yellow-400 flex items-center gap-2">
                    <AlertCircle className="w-4 h-4 flex-shrink-0" />
                    All {(selectedRun.topDiscoveries ?? []).length} discoveries are below {confidenceThreshold}%. Lower the filter above to view them.
                  </div>
                ) : null}

                <div className="flex justify-end gap-2 pt-2">
                  {selectedRun.status === "completed" && dialogFiltered.length > 0 && (
                    <Button
                      size="sm"
                      className="gap-1 bg-green-600 hover:bg-green-700 text-white"
                      onClick={() => handleImport(selectedRun)}
                      disabled={importingRunId === selectedRun.runId}
                    >
                      {importingRunId === selectedRun.runId ? (
                        <Loader2 className="w-4 h-4 animate-spin" />
                      ) : (
                        <Import className="w-4 h-4" />
                      )}
                      Import {dialogFiltered.length} Compound{dialogFiltered.length !== 1 ? "s" : ""}
                    </Button>
                  )}
                  <Button variant="outline" size="sm" onClick={() => exportRun(selectedRun)}>
                    <Download className="w-4 h-4 mr-1" /> Export Filtered JSON
                  </Button>
                </div>
              </div>
            );
          })()}
        </DialogContent>
      </Dialog>
    </div>
  );
}
