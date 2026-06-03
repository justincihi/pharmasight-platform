import { useState, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import DashboardLayout from "@/components/DashboardLayout";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import {
  Loader2, Dna, CheckCircle2, XCircle, AlertTriangle, Save,
  ChevronDown, ChevronUp, Info, Beaker, Copy, Download, History
} from "lucide-react";
import { toast } from "sonner";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";

interface SarAnalog {
  smiles: string;
  tanimoto: number;
  mw: number;
  logp: number;
  tpsa: number;
  hbd: number;
  hba: number;
  rot: number;
  qed: number;
  bbb: number;
  herg: number;
  ames: number;
  dili: number;
  lipinski_pass: boolean;
  sar_pass: boolean;
  sar_flags: string[];
  composite_score: number;
  source: string;
  isDemo: boolean;
}

interface GenerationResult {
  success: boolean;
  parent_smiles: string;
  strategy: string;
  total_generated: number;
  sar_passed: number;
  analogs: SarAnalog[];
  sar_criteria: Record<string, number | boolean>;
  has_rdkit: boolean;
  has_admet_ai: boolean;
  timestamp: string;
}

type AuditEntry = {
  id: string;
  parentSmiles: string;
  strategy: string;
  generatedAt: string;
  totalGenerated: number;
  sarPassed: number;
  savedCount: number;
};

function ScoreBadge({ value, low, high, invert = false }: { value: number; low: number; high: number; invert?: boolean }) {
  const pct = (value - low) / (high - low);
  const good = invert ? pct < 0.4 : pct > 0.6;
  const warn = invert ? pct >= 0.4 && pct < 0.7 : pct >= 0.3 && pct <= 0.6;
  return (
    <span className={`font-mono text-xs font-semibold ${good ? "text-green-500" : warn ? "text-yellow-500" : "text-red-500"}`}>
      {value.toFixed(3)}
    </span>
  );
}

export default function AnalogOfAnalog() {
  const [parentSmiles, setParentSmiles] = useState("");
  const [strategy, setStrategy] = useState<"brics" | "scaffold" | "combined">("combined");
  const [numAnalogs, setNumAnalogs] = useState(20);
  const [showSarControls, setShowSarControls] = useState(false);
  const [showAuditTrail, setShowAuditTrail] = useState(false);

  // SAR criteria
  const [minQed, setMinQed] = useState(0.3);
  const [maxHerg, setMaxHerg] = useState(0.7);
  const [maxAmes, setMaxAmes] = useState(0.6);
  const [minBbb, setMinBbb] = useState(0.2);
  const [maxMw, setMaxMw] = useState(600);
  const [maxLogp, setMaxLogp] = useState(6.0);
  const [lipinski, setLipinski] = useState(true);

  const [result, setResult] = useState<GenerationResult | null>(null);
  const [savedIds, setSavedIds] = useState<Set<string>>(new Set());
  const [auditTrail, setAuditTrail] = useState<AuditEntry[]>([]);
  const [filterSarOnly, setFilterSarOnly] = useState(false);
  const [sortBy, setSortBy] = useState<"composite_score" | "qed" | "bbb" | "herg" | "tanimoto">("composite_score");

  // Load analogs from master list for parent picker
  const { data: masterAnalogs } = trpc.analog.list.useQuery({ limit: 200, offset: 0 });

  const generateMutation = (trpc.analog as any).generateSarAnalogs.useMutation({
    onSuccess: (data: any) => {
      setResult(data);
      setSavedIds(new Set());
      const entry: AuditEntry = {
        id: Date.now().toString(),
        parentSmiles,
        strategy,
        generatedAt: new Date().toISOString(),
        totalGenerated: data.total_generated,
        sarPassed: data.sar_passed,
        savedCount: 0,
      };
      setAuditTrail((prev) => [entry, ...prev.slice(0, 49)]);
      toast.success(`Generated ${data.total_generated} analogs — ${data.sar_passed} passed SAR filters`);
    },
    onError: (err: any) => {
      toast.error(`Generation failed: ${err.message}`);
    },
  });

  const saveMutation = (trpc.analog as any).saveToMasterList.useMutation({
    onSuccess: (data: any, variables: any) => {
      setSavedIds((prev) => { const next = new Set(Array.from(prev)); next.add(variables.smiles); return next; });
      setAuditTrail((prev) =>
        prev.map((e, i) => (i === 0 ? { ...e, savedCount: e.savedCount + 1 } : e))
      );
      toast.success(`Saved to master list as ${data.compoundId}`);
    },
    onError: (err: any) => {
      toast.error(`Save failed: ${err.message}`);
    },
  });

  const handleGenerate = () => {
    if (!parentSmiles.trim()) { toast.error("Enter a parent SMILES string"); return; }
    generateMutation.mutate({
      smiles: parentSmiles.trim(),
      numAnalogs,
      strategy,
      sarCriteria: { minQed, maxHerg, maxAmes, minBbb, maxMw, maxLogp, lipinski },
    });
  };

  const handleSave = (analog: SarAnalog) => {
    saveMutation.mutate({
      smiles: analog.smiles,
      parentSmiles,
      discoveryMethod: `analog-of-analog (${strategy})`,
      admetData: {
        composite_score: analog.composite_score,
        tanimoto: analog.tanimoto,
        bbb: analog.bbb,
        herg: analog.herg,
        ames: analog.ames,
        dili: analog.dili,
        qed: analog.qed,
        logp: analog.logp,
        tpsa: analog.tpsa,
        mw: analog.mw,
      },
      notes: `Generated via ${strategy} strategy. SAR: ${analog.sar_pass ? "PASS" : "FAIL"}. Flags: ${analog.sar_flags.join("; ") || "none"}`,
    });
  };

  const handleSaveAll = () => {
    const toSave = displayedAnalogs.filter((a) => !savedIds.has(a.smiles));
    toSave.forEach((a) => handleSave(a));
    toast.info(`Saving ${toSave.length} analogs to master list...`);
  };

  const displayedAnalogs = useMemo(() => {
    if (!result) return [];
    let list = filterSarOnly ? result.analogs.filter((a) => a.sar_pass) : result.analogs;
    return [...list].sort((a, b) => (b as any)[sortBy] - (a as any)[sortBy]);
  }, [result, filterSarOnly, sortBy]);

  const exportCSV = () => {
    if (!displayedAnalogs.length) return;
    const headers = ["SMILES", "Tanimoto", "MW", "LogP", "TPSA", "QED", "BBB", "hERG", "AMES", "DILI", "Lipinski", "SAR Pass", "SAR Flags", "Composite Score", "Source"];
    const rows = displayedAnalogs.map((a) => [
      `"${a.smiles}"`, a.tanimoto, a.mw, a.logp, a.tpsa, a.qed, a.bbb, a.herg, a.ames, a.dili,
      a.lipinski_pass ? "Y" : "N", a.sar_pass ? "Y" : "N",
      `"${a.sar_flags.join("; ")}"`, a.composite_score, a.source,
    ].join(","));
    const csv = [headers.join(","), ...rows].join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `sar_analogs_${new Date().toISOString().slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(url);
    toast.success("CSV exported");
  };

  return (
    <TooltipProvider>
      <DashboardLayout>
        <div className="space-y-4 p-4">
          {/* Header */}
          <div>
            <h1 className="text-2xl font-bold flex items-center gap-2">
              <Dna className="h-6 w-6 text-purple-500" />
              Analog-of-Analog Generator
            </h1>
            <p className="text-sm text-muted-foreground mt-1">
              Select a parent compound, configure SAR criteria, and generate ML-scored analogs.
              Click any result to save it to the master list — a full audit trail is maintained.
            </p>
          </div>

          {/* Input Panel */}
          <Card>
            <CardHeader className="pb-3">
              <CardTitle className="text-base">Parent Compound</CardTitle>
              <CardDescription>Enter a SMILES string or select from the master list</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="flex gap-2 flex-wrap">
                <div className="flex-1 min-w-[240px]">
                  <Label className="text-xs mb-1 block">SMILES</Label>
                  <Input
                    placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O"
                    value={parentSmiles}
                    onChange={(e) => setParentSmiles(e.target.value)}
                    className="font-mono text-sm"
                  />
                </div>
                <div className="min-w-[200px]">
                  <Label className="text-xs mb-1 block">Or select from master list</Label>
                  <Select onValueChange={(v) => setParentSmiles(v)}>
                    <SelectTrigger className="text-sm">
                      <SelectValue placeholder="Pick a compound..." />
                    </SelectTrigger>
                    <SelectContent>
                      {(masterAnalogs as any)?.analogs?.map((a: any) => (
                        <SelectItem key={a.id} value={a.smiles ?? ""}>
                          <span className="font-mono text-xs">{a.compoundId}</span>
                          {a.compoundName && <span className="ml-2 text-muted-foreground">{a.compoundName}</span>}
                        </SelectItem>
                      ))}
                    </SelectContent>
                  </Select>
                </div>
                <div className="min-w-[140px]">
                  <Label className="text-xs mb-1 block">Strategy</Label>
                  <Select value={strategy} onValueChange={(v: any) => setStrategy(v)}>
                    <SelectTrigger className="text-sm">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="brics">BRICS</SelectItem>
                      <SelectItem value="scaffold">Scaffold</SelectItem>
                      <SelectItem value="combined">Combined</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <div className="min-w-[120px]">
                  <Label className="text-xs mb-1 block">Max Analogs: {numAnalogs}</Label>
                  <Slider
                    min={5} max={50} step={5}
                    value={[numAnalogs]}
                    onValueChange={([v]) => setNumAnalogs(v)}
                    className="mt-2"
                  />
                </div>
              </div>

              {/* SAR Criteria Toggle */}
              <button
                className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground transition-colors"
                onClick={() => setShowSarControls((v) => !v)}
              >
                {showSarControls ? <ChevronUp className="h-3.5 w-3.5" /> : <ChevronDown className="h-3.5 w-3.5" />}
                SAR Filter Criteria
              </button>

              {showSarControls && (
                <div className="grid grid-cols-2 md:grid-cols-4 gap-4 p-3 bg-muted/30 rounded-lg border">
                  {[
                    { label: "Min QED", value: minQed, set: setMinQed, min: 0, max: 1, step: 0.05, desc: "Drug-likeness (higher = better)" },
                    { label: "Max hERG", value: maxHerg, set: setMaxHerg, min: 0, max: 1, step: 0.05, desc: "Cardiotoxicity risk (lower = safer)" },
                    { label: "Max AMES", value: maxAmes, set: setMaxAmes, min: 0, max: 1, step: 0.05, desc: "Mutagenicity (lower = safer)" },
                    { label: "Min BBB", value: minBbb, set: setMinBbb, min: 0, max: 1, step: 0.05, desc: "Blood-brain barrier permeability" },
                    { label: "Max LogP", value: maxLogp, set: setMaxLogp, min: -2, max: 10, step: 0.5, desc: "Lipophilicity (Lipinski ≤ 5)" },
                    { label: `Max MW: ${maxMw}`, value: maxMw, set: setMaxMw, min: 200, max: 1000, step: 50, desc: "Molecular weight (Da)" },
                  ].map(({ label, value, set, min, max, step, desc }) => (
                    <div key={label}>
                      <Tooltip>
                        <TooltipTrigger asChild>
                          <Label className="text-xs mb-1 flex items-center gap-1 cursor-help">
                            {label}: <span className="font-mono font-semibold">{value}</span>
                            <Info className="h-3 w-3 text-muted-foreground" />
                          </Label>
                        </TooltipTrigger>
                        <TooltipContent><p className="text-xs">{desc}</p></TooltipContent>
                      </Tooltip>
                      <Slider
                        min={min} max={max} step={step}
                        value={[value]}
                        onValueChange={([v]) => set(v)}
                      />
                    </div>
                  ))}
                  <div className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      id="lipinski"
                      checked={lipinski}
                      onChange={(e) => setLipinski(e.target.checked)}
                      className="rounded"
                    />
                    <Label htmlFor="lipinski" className="text-xs cursor-pointer">Lipinski Ro5</Label>
                  </div>
                </div>
              )}

              <div className="flex gap-2">
                <Button
                  onClick={handleGenerate}
                  disabled={generateMutation.isPending || !parentSmiles.trim()}
                  className="gap-2"
                >
                  {generateMutation.isPending ? (
                    <><Loader2 className="h-4 w-4 animate-spin" /> Generating...</>
                  ) : (
                    <><Beaker className="h-4 w-4" /> Generate Analogs</>
                  )}
                </Button>
                {result && (
                  <Button variant="outline" size="sm" onClick={() => setShowAuditTrail((v) => !v)} className="gap-1">
                    <History className="h-4 w-4" />
                    Audit Trail ({auditTrail.length})
                  </Button>
                )}
              </div>
            </CardContent>
          </Card>

          {/* Audit Trail */}
          {showAuditTrail && auditTrail.length > 0 && (
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm">Generation Audit Trail</CardTitle>
                <CardDescription className="text-xs">All generation runs this session, with parent SMILES and save counts</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-1.5">
                  {auditTrail.map((e) => (
                    <div key={e.id} className="flex items-center gap-3 text-xs p-2 rounded bg-muted/30 border">
                      <span className="text-muted-foreground whitespace-nowrap">{new Date(e.generatedAt).toLocaleTimeString()}</span>
                      <span className="font-mono truncate max-w-[200px]" title={e.parentSmiles}>{e.parentSmiles}</span>
                      <Badge variant="outline" className="text-xs">{e.strategy}</Badge>
                      <span className="text-muted-foreground">{e.totalGenerated} generated</span>
                      <span className="text-green-600">{e.sarPassed} SAR pass</span>
                      {e.savedCount > 0 && <span className="text-blue-600 font-semibold">{e.savedCount} saved</span>}
                      <button
                        className="ml-auto text-muted-foreground hover:text-foreground"
                        onClick={() => setParentSmiles(e.parentSmiles)}
                        title="Re-run with this parent"
                      >
                        ↩
                      </button>
                    </div>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Results */}
          {result && (
            <Card>
              <CardHeader className="pb-3">
                <div className="flex items-center justify-between flex-wrap gap-2">
                  <div>
                    <CardTitle className="text-base flex items-center gap-2">
                      Results
                      {(result as any).isDemo === false && result.has_rdkit && (
                        <Badge className="bg-green-600 text-white text-xs">Real RDKit</Badge>
                      )}
                      {result.has_admet_ai && (
                        <Badge className="bg-blue-600 text-white text-xs">ADMET-AI</Badge>
                      )}
                    </CardTitle>
                    <CardDescription className="text-xs">
                      {result.total_generated} analogs generated · {result.sar_passed} passed SAR filters ·
                      Strategy: {result.strategy} · {new Date(result.timestamp).toLocaleTimeString()}
                    </CardDescription>
                  </div>
                  <div className="flex items-center gap-2 flex-wrap">
                    <button
                      className={`text-xs px-2 py-1 rounded border transition-colors ${filterSarOnly ? "bg-primary text-primary-foreground border-primary" : "border-border hover:border-primary/50"}`}
                      onClick={() => setFilterSarOnly((v) => !v)}
                    >
                      SAR Pass Only ({result.sar_passed})
                    </button>
                    <Select value={sortBy} onValueChange={(v: any) => setSortBy(v)}>
                      <SelectTrigger className="h-7 text-xs w-36">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="composite_score">Composite Score</SelectItem>
                        <SelectItem value="qed">QED</SelectItem>
                        <SelectItem value="bbb">BBB</SelectItem>
                        <SelectItem value="herg">hERG (asc)</SelectItem>
                        <SelectItem value="tanimoto">Tanimoto</SelectItem>
                      </SelectContent>
                    </Select>
                    <Button variant="outline" size="sm" onClick={exportCSV} className="h-7 text-xs gap-1">
                      <Download className="h-3 w-3" /> CSV
                    </Button>
                    <Button size="sm" onClick={handleSaveAll} className="h-7 text-xs gap-1"
                      disabled={saveMutation.isPending || displayedAnalogs.every((a) => savedIds.has(a.smiles))}>
                      <Save className="h-3 w-3" /> Save All Shown
                    </Button>
                  </div>
                </div>
              </CardHeader>
              <CardContent>
                <div className="overflow-x-auto rounded border">
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-muted/50 border-b">
                        <th className="px-2 py-2 text-left font-semibold">SMILES</th>
                        <th className="px-2 py-2 text-center font-semibold">Tanimoto</th>
                        <th className="px-2 py-2 text-center font-semibold">MW</th>
                        <th className="px-2 py-2 text-center font-semibold">LogP</th>
                        <th className="px-2 py-2 text-center font-semibold">QED</th>
                        <th className="px-2 py-2 text-center font-semibold">BBB</th>
                        <th className="px-2 py-2 text-center font-semibold">hERG</th>
                        <th className="px-2 py-2 text-center font-semibold">AMES</th>
                        <th className="px-2 py-2 text-center font-semibold">DILI</th>
                        <th className="px-2 py-2 text-center font-semibold">Ro5</th>
                        <th className="px-2 py-2 text-center font-semibold">SAR</th>
                        <th className="px-2 py-2 text-center font-semibold">Score</th>
                        <th className="px-2 py-2 text-center font-semibold">Actions</th>
                      </tr>
                    </thead>
                    <tbody>
                      {displayedAnalogs.map((analog, i) => {
                        const isSaved = savedIds.has(analog.smiles);
                        return (
                          <tr key={i} className={`border-b hover:bg-muted/20 transition-colors ${isSaved ? "bg-green-500/5" : ""}`}>
                            <td className="px-2 py-1.5 max-w-[200px]">
                              <div className="flex items-center gap-1">
                                <span className="font-mono truncate text-xs" title={analog.smiles}>
                                  {analog.smiles.length > 30 ? analog.smiles.slice(0, 30) + "…" : analog.smiles}
                                </span>
                                <button
                                  className="text-muted-foreground hover:text-foreground flex-shrink-0"
                                  onClick={() => { navigator.clipboard.writeText(analog.smiles); toast.success("SMILES copied"); }}
                                  title="Copy SMILES"
                                >
                                  <Copy className="h-3 w-3" />
                                </button>
                              </div>
                              {analog.sar_flags.length > 0 && (
                                <div className="flex flex-wrap gap-0.5 mt-0.5">
                                  {analog.sar_flags.map((f, fi) => (
                                    <span key={fi} className="text-[10px] text-yellow-600 bg-yellow-500/10 px-1 rounded">{f}</span>
                                  ))}
                                </div>
                              )}
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.tanimoto} low={0} high={1} />
                            </td>
                            <td className="px-2 py-1.5 text-center font-mono text-xs">
                              <span className={analog.mw > 500 ? "text-yellow-500" : "text-foreground"}>{analog.mw.toFixed(0)}</span>
                            </td>
                            <td className="px-2 py-1.5 text-center font-mono text-xs">
                              <span className={analog.logp > 5 ? "text-yellow-500" : "text-foreground"}>{analog.logp.toFixed(2)}</span>
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.qed} low={0} high={1} />
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.bbb} low={0} high={1} />
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.herg} low={0} high={1} invert />
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.ames} low={0} high={1} invert />
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <ScoreBadge value={analog.dili} low={0} high={1} invert />
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              {analog.lipinski_pass
                                ? <CheckCircle2 className="h-3.5 w-3.5 text-green-500 mx-auto" />
                                : <XCircle className="h-3.5 w-3.5 text-red-500 mx-auto" />}
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              {analog.sar_pass
                                ? <CheckCircle2 className="h-3.5 w-3.5 text-green-500 mx-auto" />
                                : <AlertTriangle className="h-3.5 w-3.5 text-yellow-500 mx-auto" />}
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              <span className="font-mono font-semibold text-blue-500">{analog.composite_score.toFixed(3)}</span>
                            </td>
                            <td className="px-2 py-1.5 text-center">
                              {isSaved ? (
                                <span className="text-green-600 text-xs font-semibold flex items-center gap-0.5 justify-center">
                                  <CheckCircle2 className="h-3.5 w-3.5" /> Saved
                                </span>
                              ) : (
                                <Button
                                  size="sm"
                                  variant="outline"
                                  className="h-6 text-xs px-2 gap-1"
                                  onClick={() => handleSave(analog)}
                                  disabled={saveMutation.isPending}
                                >
                                  <Save className="h-3 w-3" /> Save
                                </Button>
                              )}
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
                <p className="text-xs text-muted-foreground mt-2">
                  {displayedAnalogs.length} analogs shown · Scores: green ≥ 0.6, yellow 0.3–0.6, red &lt; 0.3 ·
                  {(result as any).isDemo ? " Demo mode (RDKit/ADMET-AI not available in sandbox)" : " Real RDKit predictions"}
                </p>
              </CardContent>
            </Card>
          )}

          {/* Empty state */}
          {!result && !generateMutation.isPending && (
            <Card>
              <CardContent className="py-12 text-center text-muted-foreground">
                <Dna className="h-10 w-10 mx-auto mb-3 opacity-30" />
                <p className="font-medium">No analogs generated yet</p>
                <p className="text-sm mt-1">Enter a parent SMILES string and click Generate Analogs to begin.</p>
              </CardContent>
            </Card>
          )}
        </div>
      </DashboardLayout>
    </TooltipProvider>
  );
}
