import { useState, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Checkbox } from "@/components/ui/checkbox";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { Loader2, Play, Download, CheckCircle2, XCircle, FlaskConical, BarChart3, AlertTriangle } from "lucide-react";
import { toast } from "sonner";

interface BatchAdmetPanelProps {
  /** Pre-selected analog IDs (optional). If omitted, user can select from the list. */
  preSelectedIds?: number[];
}

type RunResult = { analogId: number; status: "ok" | "error"; error?: string };

function riskBadge(val: string | null | undefined, threshold = 0.5) {
  if (val == null) return <span className="text-muted-foreground text-xs">—</span>;
  const n = parseFloat(val);
  if (isNaN(n)) return <span className="text-xs">{val}</span>;
  const pct = (n * 100).toFixed(0);
  if (n > threshold) return <Badge variant="destructive" className="text-xs px-1 py-0">{pct}%</Badge>;
  if (n > threshold * 0.6) return <Badge variant="secondary" className="text-xs px-1 py-0 bg-yellow-500/20 text-yellow-700 dark:text-yellow-400">{pct}%</Badge>;
  return <Badge variant="outline" className="text-xs px-1 py-0 text-green-600 dark:text-green-400 border-green-500/40">{pct}%</Badge>;
}

function numCell(val: string | null | undefined, decimals = 2) {
  if (val == null) return <span className="text-muted-foreground text-xs">—</span>;
  const n = parseFloat(val);
  if (isNaN(n)) return <span className="text-xs">{val}</span>;
  return <span className="text-xs font-mono">{n.toFixed(decimals)}</span>;
}

export function BatchAdmetPanel({ preSelectedIds }: BatchAdmetPanelProps) {
  const [selectedIds, setSelectedIds] = useState<number[]>(preSelectedIds ?? []);
  const [runResults, setRunResults] = useState<RunResult[] | null>(null);
  const [isRunning, setIsRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [lastBatchId, setLastBatchId] = useState<string | null>(null);

  const { data: analogs, isLoading: analogsLoading } = trpc.analog.list.useQuery({ limit: 500, offset: 0 });
  const { data: admetStats, refetch: refetchStats } = trpc.analog.getAdmetStats.useQuery();
  const runBatchMutation = trpc.analog.runBatchAdmet.useMutation();

  const toggleId = (id: number) =>
    setSelectedIds((prev) => (prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]));

  const selectAll = () => analogs && setSelectedIds(analogs.map((a: any) => a.id));
  const deselectAll = () => setSelectedIds([]);

  const handleRun = async () => {
    if (selectedIds.length === 0) {
      toast.error("Select at least one analog to screen.");
      return;
    }
    setIsRunning(true);
    setProgress(0);
    setRunResults(null);
    // Simulate progress ticks while the batch runs
    const tick = setInterval(() => setProgress((p) => Math.min(p + 2, 90)), 800);
    try {
      const res = await runBatchMutation.mutateAsync({ analogIds: selectedIds });
      clearInterval(tick);
      setProgress(100);
      setRunResults(res.results);
      setLastBatchId(res.batchRunId);
      await refetchStats();
      toast.success(`ADMET screening complete: ${res.succeeded}/${res.total} succeeded.`);
    } catch (err: any) {
      clearInterval(tick);
      toast.error(`Batch failed: ${err.message}`);
    } finally {
      setIsRunning(false);
    }
  };

  // Fetch stored ADMET results for selected analogs (latest per analog)
  const { data: storedResults } = trpc.analog.getAdmetResults.useQuery(
    { analogId: selectedIds[0] ?? 0 },
    { enabled: selectedIds.length === 1 }
  );

  // CSV export of run results
  const exportCsv = () => {
    if (!runResults) return;
    const analogMap = new Map((analogs ?? []).map((a: any) => [a.id, a.compoundName]));
    const header = "analogId,compoundName,status,error";
    const rows = runResults.map((r) =>
      `${r.analogId},"${analogMap.get(r.analogId) ?? r.analogId}",${r.status},"${r.error ?? ""}"`
    );
    const blob = new Blob([[header, ...rows].join("\n")], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `admet-batch-${lastBatchId ?? Date.now()}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const succeeded = runResults?.filter((r) => r.status === "ok").length ?? 0;
  const failed = runResults?.filter((r) => r.status === "error").length ?? 0;

  return (
    <TooltipProvider>
      <div className="space-y-6">
        {/* Stats bar */}
        {admetStats && (
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
            {[
              { label: "Screened", value: admetStats.total, icon: <FlaskConical className="h-4 w-4" /> },
              { label: "hERG Risk", value: `${admetStats.flaggedHerg} flagged`, icon: <AlertTriangle className="h-4 w-4 text-red-500" /> },
              { label: "AMES Risk", value: `${admetStats.flaggedAmes} flagged`, icon: <AlertTriangle className="h-4 w-4 text-orange-500" /> },
              { label: "Avg BBB", value: admetStats.avgBbb != null ? admetStats.avgBbb.toFixed(2) : "—", icon: <BarChart3 className="h-4 w-4 text-blue-500" /> },
            ].map((s) => (
              <Card key={s.label} className="py-3">
                <CardContent className="px-4 flex items-center gap-3">
                  {s.icon}
                  <div>
                    <p className="text-xs text-muted-foreground">{s.label}</p>
                    <p className="text-sm font-semibold">{s.value}</p>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )}

        {/* Selection + run controls */}
        {!preSelectedIds && (
          <Card>
            <CardHeader className="pb-3">
              <div className="flex items-center justify-between flex-wrap gap-2">
                <div>
                  <CardTitle className="text-base">Select Analogs to Screen</CardTitle>
                  <CardDescription>{selectedIds.length} of {analogs?.length ?? 0} selected</CardDescription>
                </div>
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" onClick={selectAll}>All</Button>
                  <Button variant="outline" size="sm" onClick={deselectAll}>None</Button>
                </div>
              </div>
            </CardHeader>
            <CardContent>
              {analogsLoading ? (
                <div className="flex justify-center py-6"><Loader2 className="h-6 w-6 animate-spin text-muted-foreground" /></div>
              ) : (
                <div className="max-h-52 overflow-y-auto space-y-1 pr-1">
                  {(analogs ?? []).map((a: any) => (
                    <label key={a.id} className="flex items-center gap-2 px-2 py-1 rounded hover:bg-muted/50 cursor-pointer text-sm">
                      <Checkbox
                        checked={selectedIds.includes(a.id)}
                        onCheckedChange={() => toggleId(a.id)}
                      />
                      <span className="flex-1 truncate">{a.compoundName}</span>
                      <span className="text-xs text-muted-foreground font-mono truncate max-w-[120px]">{a.smiles?.slice(0, 18)}{a.smiles?.length > 18 ? "…" : ""}</span>
                    </label>
                  ))}
                </div>
              )}
            </CardContent>
          </Card>
        )}

        {/* Run button + progress */}
        <div className="flex items-center gap-3 flex-wrap">
          <Button onClick={handleRun} disabled={isRunning || selectedIds.length === 0} className="gap-2">
            {isRunning ? <Loader2 className="h-4 w-4 animate-spin" /> : <Play className="h-4 w-4" />}
            {isRunning ? `Screening ${selectedIds.length} analogs…` : `Run ADMET-AI on ${selectedIds.length} analog${selectedIds.length !== 1 ? "s" : ""}`}
          </Button>
          {runResults && (
            <Button variant="outline" size="sm" onClick={exportCsv} className="gap-2">
              <Download className="h-4 w-4" />
              Export CSV
            </Button>
          )}
        </div>

        {isRunning && (
          <div className="space-y-1">
            <Progress value={progress} className="h-2" />
            <p className="text-xs text-muted-foreground">Running ML predictions via ADMET-AI (Chemprop)…</p>
          </div>
        )}

        {/* Results summary */}
        {runResults && (
          <div className="space-y-4">
            <div className="flex items-center gap-4 text-sm">
              <span className="flex items-center gap-1 text-green-600 dark:text-green-400"><CheckCircle2 className="h-4 w-4" />{succeeded} succeeded</span>
              {failed > 0 && <span className="flex items-center gap-1 text-red-500"><XCircle className="h-4 w-4" />{failed} failed</span>}
              {lastBatchId && <span className="text-xs text-muted-foreground font-mono">Batch: {lastBatchId}</span>}
            </div>

            {/* Results table */}
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm">Screening Results</CardTitle>
                <CardDescription className="text-xs">Stored to database — accessible from Compound Testing → ADMET tab per analog</CardDescription>
              </CardHeader>
              <CardContent className="p-0">
                <div className="overflow-x-auto">
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead className="text-xs">Compound</TableHead>
                        <TableHead className="text-xs">Status</TableHead>
                        <TableHead className="text-xs">
                          <Tooltip><TooltipTrigger>hERG</TooltipTrigger><TooltipContent>Cardiotoxicity risk (hERG channel inhibition)</TooltipContent></Tooltip>
                        </TableHead>
                        <TableHead className="text-xs">
                          <Tooltip><TooltipTrigger>AMES</TooltipTrigger><TooltipContent>Mutagenicity (Ames test probability)</TooltipContent></Tooltip>
                        </TableHead>
                        <TableHead className="text-xs">
                          <Tooltip><TooltipTrigger>DILI</TooltipTrigger><TooltipContent>Drug-induced liver injury risk</TooltipContent></Tooltip>
                        </TableHead>
                        <TableHead className="text-xs">
                          <Tooltip><TooltipTrigger>BBB</TooltipTrigger><TooltipContent>Blood-brain barrier penetration probability</TooltipContent></Tooltip>
                        </TableHead>
                        <TableHead className="text-xs">
                          <Tooltip><TooltipTrigger>Oral Bio</TooltipTrigger><TooltipContent>Oral bioavailability (%F ≥ 20%)</TooltipContent></Tooltip>
                        </TableHead>
                        <TableHead className="text-xs">QED</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {runResults.map((r) => {
                        const analog = (analogs ?? []).find((a: any) => a.id === r.analogId);
                        return (
                          <TableRow key={r.analogId}>
                            <TableCell className="text-xs font-medium max-w-[140px] truncate">{analog?.compoundName ?? `#${r.analogId}`}</TableCell>
                            <TableCell>
                              {r.status === "ok"
                                ? <Badge variant="outline" className="text-xs text-green-600 border-green-500/40">OK</Badge>
                                : <Tooltip><TooltipTrigger><Badge variant="destructive" className="text-xs">Error</Badge></TooltipTrigger><TooltipContent>{r.error}</TooltipContent></Tooltip>
                              }
                            </TableCell>
                            <TableCell>{r.status === "ok" ? riskBadge(null) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                            <TableCell>{r.status === "ok" ? riskBadge(null) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                            <TableCell>{r.status === "ok" ? riskBadge(null) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                            <TableCell>{r.status === "ok" ? riskBadge(null, 0.3) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                            <TableCell>{r.status === "ok" ? riskBadge(null, 0.3) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                            <TableCell>{r.status === "ok" ? numCell(null) : <span className="text-muted-foreground text-xs">—</span>}</TableCell>
                          </TableRow>
                        );
                      })}
                    </TableBody>
                  </Table>
                </div>
                <p className="text-xs text-muted-foreground px-4 pb-3 pt-1">
                  Detailed per-property values are stored in the database. Open a compound in Compound Testing → ADMET tab to view all 49 ML properties.
                </p>
              </CardContent>
            </Card>
          </div>
        )}
      </div>
    </TooltipProvider>
  );
}
