import { useState, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import DashboardLayout from "@/components/DashboardLayout";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import {
  Loader2, Download, Search, FileJson, FileText, Table2,
  ArrowUpDown, ArrowUp, ArrowDown, CheckCircle2, Circle, Info
} from "lucide-react";
import { toast } from "sonner";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";

// Column definitions for the master sheet
const COLUMNS: Array<{
  key: string;
  label: string;
  group: string;
  description?: string;
  format?: (v: any) => string;
}> = [
  // Identity
  { key: "compoundId", label: "Compound ID", group: "Identity" },
  { key: "compoundName", label: "Name", group: "Identity" },
  { key: "smiles", label: "SMILES", group: "Identity" },
  { key: "parentCompound", label: "Parent", group: "Identity" },
  { key: "discoveryMethod", label: "Method", group: "Identity" },
  { key: "patentStatus", label: "Patent", group: "Identity" },
  { key: "fdaStatus", label: "FDA Status", group: "Identity" },
  // Scores
  { key: "confidenceScore", label: "Confidence", group: "Scores", format: (v) => v != null ? `${Number(v).toFixed(0)}%` : "—" },
  { key: "similarityScore", label: "Similarity", group: "Scores", format: (v) => v != null ? `${Number(v).toFixed(0)}%` : "—" },
  { key: "safetyScore", label: "Safety", group: "Scores", format: (v) => v != null ? `${Number(v).toFixed(0)}` : "—" },
  { key: "efficacyScore", label: "Efficacy", group: "Scores", format: (v) => v != null ? `${Number(v).toFixed(0)}` : "—" },
  { key: "drugLikenessScore", label: "Drug-like", group: "Scores", format: (v) => v != null ? `${Number(v).toFixed(0)}` : "—" },
  // ADMET
  { key: "bbbPermeability", label: "BBB", group: "ADMET", description: "Blood-brain barrier permeability (0–1)", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "oralBioavailability", label: "Oral BA", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "hia", label: "HIA", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "caco2", label: "Caco-2", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "pgp", label: "P-gp", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "ppbr", label: "PPBR (%)", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(1) : "—" },
  { key: "halfLife", label: "t½ (h)", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(1) : "—" },
  { key: "clearanceHepatocyte", label: "Clearance", group: "ADMET", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  // Toxicity
  { key: "herg", label: "hERG", group: "Toxicity", description: "hERG cardiotoxicity risk (0–1)", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "ames", label: "AMES", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "dili", label: "DILI", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "ld50", label: "LD50", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "clintox", label: "ClinTox", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "hepatotoxicity", label: "Hepatotox", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "mutagenicity", label: "Mutagenic", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "carcinogenicity", label: "Carcino", group: "Toxicity", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  // CYP
  { key: "cyp1a2", label: "CYP1A2", group: "CYP", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "cyp2c9", label: "CYP2C9", group: "CYP", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "cyp2c19", label: "CYP2C19", group: "CYP", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "cyp2d6", label: "CYP2D6", group: "CYP", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "cyp3a4", label: "CYP3A4", group: "CYP", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  // Physicochemical
  { key: "molecularWeight", label: "MW", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(1) : "—" },
  { key: "logp", label: "LogP", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(2) : "—" },
  { key: "tpsa", label: "TPSA (Å²)", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(1) : "—" },
  { key: "qed", label: "QED", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "solubility", label: "Solubility", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(3) : "—" },
  { key: "lipophilicity", label: "Lipophilicity", group: "Physchem", format: (v) => v != null ? Number(v).toFixed(2) : "—" },
  // Docking
  { key: "bestDockingScore", label: "Docking (kcal/mol)", group: "Docking", format: (v) => v != null ? Number(v).toFixed(2) : "—" },
  { key: "bestReceptor", label: "Best Receptor", group: "Docking" },
  // Audit
  { key: "admetRunAt", label: "ADMET Run", group: "Audit", format: (v) => v ? new Date(v).toLocaleDateString() : "—" },
  { key: "toxicityRunAt", label: "Tox Run", group: "Audit", format: (v) => v ? new Date(v).toLocaleDateString() : "—" },
  { key: "dockingRunAt", label: "Dock Run", group: "Audit", format: (v) => v ? new Date(v).toLocaleDateString() : "—" },
  { key: "discoveredAt", label: "Discovered", group: "Audit", format: (v) => v ? new Date(v).toLocaleDateString() : "—" },
];

const GROUPS = ["Identity", "Scores", "ADMET", "Toxicity", "CYP", "Physchem", "Docking", "Audit"];

type SortDir = "asc" | "desc";

function completionRate(row: any): number {
  const admetKeys = ["bbbPermeability", "herg", "ames", "dili", "qed", "tpsa", "logp"];
  const filled = admetKeys.filter((k) => row[k] != null).length;
  return Math.round((filled / admetKeys.length) * 100);
}

export default function MasterCompoundSheet() {
  const [sortKey, setSortKey] = useState("discoveredAt");
  const [sortDir, setSortDir] = useState<SortDir>("desc");
  const [search, setSearch] = useState("");
  const [visibleGroups, setVisibleGroups] = useState<Set<string>>(
    new Set(["Identity", "Scores", "ADMET", "Toxicity", "Docking"])
  );

  const { data: rows, isLoading } = trpc.analog.getMasterSheet.useQuery();

  const visibleCols = COLUMNS.filter((c) => visibleGroups.has(c.group));

  const filtered = useMemo(() => {
    if (!rows) return [];
    const q = search.toLowerCase();
    return rows.filter((r: any) =>
      !q ||
      (r.compoundId ?? "").toLowerCase().includes(q) ||
      (r.compoundName ?? "").toLowerCase().includes(q) ||
      (r.smiles ?? "").toLowerCase().includes(q) ||
      (r.parentCompound ?? "").toLowerCase().includes(q)
    );
  }, [rows, search]);

  const sorted = useMemo(() => {
    return [...filtered].sort((a: any, b: any) => {
      let av = a[sortKey];
      let bv = b[sortKey];
      if (av == null) av = sortDir === "asc" ? "\uffff" : "";
      if (bv == null) bv = sortDir === "asc" ? "\uffff" : "";
      const numA = parseFloat(av);
      const numB = parseFloat(bv);
      if (!isNaN(numA) && !isNaN(numB)) return sortDir === "asc" ? numA - numB : numB - numA;
      return sortDir === "asc" ? String(av).localeCompare(String(bv)) : String(bv).localeCompare(String(av));
    });
  }, [filtered, sortKey, sortDir]);

  const handleSort = (key: string) => {
    if (sortKey === key) setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    else { setSortKey(key); setSortDir("asc"); }
  };

  const SortIcon = ({ col }: { col: string }) =>
    sortKey !== col ? <ArrowUpDown className="h-3 w-3 opacity-40" /> :
    sortDir === "asc" ? <ArrowUp className="h-3 w-3" /> : <ArrowDown className="h-3 w-3" />;

  const exportCSV = () => {
    if (!sorted.length) return;
    const headers = visibleCols.map((c) => c.label);
    const csvRows = sorted.map((r: any) =>
      visibleCols.map((c) => {
        const v = r[c.key];
        const formatted = c.format ? c.format(v) : (v ?? "");
        return `"${String(formatted).replace(/"/g, '""')}"`;
      }).join(",")
    );
    const csv = [headers.join(","), ...csvRows].join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `master_compounds_${new Date().toISOString().slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(url);
    toast.success("CSV exported");
  };

  const exportJSON = () => {
    if (!sorted.length) return;
    const blob = new Blob([JSON.stringify(sorted, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `master_compounds_${new Date().toISOString().slice(0, 10)}.json`;
    a.click();
    URL.revokeObjectURL(url);
    toast.success("JSON exported");
  };

  const exportPDF = () => {
    const win = window.open("", "_blank");
    if (!win) { toast.error("Popup blocked — allow popups to export PDF"); return; }
    const rows_html = sorted.map((r: any) =>
      `<tr>${visibleCols.map((c) => {
        const v = r[c.key];
        const formatted = c.format ? c.format(v) : (v ?? "—");
        return `<td style="border:1px solid #ddd;padding:4px 6px;font-size:10px;white-space:nowrap">${formatted}</td>`;
      }).join("")}</tr>`
    ).join("");
    win.document.write(`<!DOCTYPE html><html><head><title>Master Compound Sheet</title>
      <style>body{font-family:sans-serif;font-size:11px}table{border-collapse:collapse;width:100%}th{background:#1e293b;color:#fff;padding:5px 6px;font-size:10px;text-align:left;white-space:nowrap}h1{font-size:16px;margin-bottom:4px}p{color:#666;font-size:11px;margin-bottom:12px}</style>
      </head><body>
      <h1>PharmaSight™ Master Compound Sheet</h1>
      <p>Generated ${new Date().toLocaleString()} · ${sorted.length} compounds</p>
      <table><thead><tr>${visibleCols.map((c) => `<th>${c.label}</th>`).join("")}</tr></thead>
      <tbody>${rows_html}</tbody></table>
      <script>window.onload=()=>{window.print();}</script></body></html>`);
    win.document.close();
    toast.success("PDF print dialog opened");
  };

  const toggleGroup = (g: string) => {
    setVisibleGroups((prev) => {
      const next = new Set(prev);
      if (next.has(g)) { if (next.size > 1) next.delete(g); }
      else next.add(g);
      return next;
    });
  };

  const totalFields = COLUMNS.length;
  const filledFields = rows ? rows.reduce((sum: number, r: any) => {
    return sum + COLUMNS.filter((c) => r[c.key] != null && r[c.key] !== "").length;
  }, 0) : 0;
  const totalPossible = (rows?.length ?? 0) * totalFields;
  const overallCompletion = totalPossible > 0 ? Math.round((filledFields / totalPossible) * 100) : 0;

  return (
    <TooltipProvider>
      <DashboardLayout>
        <div className="space-y-4 p-4">
          {/* Header */}
          <div className="flex items-start justify-between gap-4 flex-wrap">
            <div>
              <h1 className="text-2xl font-bold flex items-center gap-2">
                <Table2 className="h-6 w-6 text-blue-500" />
                Master Compound Sheet
              </h1>
              <p className="text-sm text-muted-foreground mt-1">
                Canonical data record for all analogs. Empty fields (—) indicate measurements not yet run.
                {rows && <span className="ml-2 font-medium">{rows.length} compounds · {overallCompletion}% data completeness</span>}
              </p>
            </div>
            <div className="flex items-center gap-2 flex-wrap">
              <Button variant="outline" size="sm" onClick={exportCSV} disabled={!sorted.length}>
                <Download className="h-4 w-4 mr-1" /> CSV
              </Button>
              <Button variant="outline" size="sm" onClick={exportJSON} disabled={!sorted.length}>
                <FileJson className="h-4 w-4 mr-1" /> JSON
              </Button>
              <Button variant="outline" size="sm" onClick={exportPDF} disabled={!sorted.length}>
                <FileText className="h-4 w-4 mr-1" /> PDF
              </Button>
            </div>
          </div>

          {/* Column group selector */}
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Column Groups</CardTitle>
              <CardDescription className="text-xs">Toggle groups to show/hide columns</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="flex flex-wrap gap-1.5">
                {GROUPS.map((g) => (
                  <button
                    key={g}
                    onClick={() => toggleGroup(g)}
                    className={`px-2 py-0.5 rounded text-xs font-medium border transition-colors ${
                      visibleGroups.has(g)
                        ? "bg-primary text-primary-foreground border-primary"
                        : "bg-transparent text-muted-foreground border-border hover:border-primary/50"
                    }`}
                  >
                    {g}
                  </button>
                ))}
              </div>
            </CardContent>
          </Card>

          {/* Search */}
          <div className="relative max-w-xs">
            <Search className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground" />
            <Input
              placeholder="Search compounds..."
              value={search}
              onChange={(e) => setSearch(e.target.value)}
              className="pl-8 h-8 text-sm"
            />
          </div>

          {/* Table */}
          {isLoading ? (
            <div className="flex items-center justify-center py-16">
              <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
              <span className="ml-2 text-muted-foreground">Loading master compound sheet...</span>
            </div>
          ) : sorted.length === 0 ? (
            <Card>
              <CardContent className="py-12 text-center text-muted-foreground">
                <Table2 className="h-10 w-10 mx-auto mb-3 opacity-30" />
                <p className="font-medium">No compounds yet</p>
                <p className="text-sm mt-1">Add analogs to the master list to see them here.</p>
              </CardContent>
            </Card>
          ) : (
            <div className="overflow-x-auto rounded-lg border">
              <table className="w-full text-xs border-collapse">
                <thead>
                  <tr className="bg-muted/50 border-b">
                    <th className="sticky left-0 z-10 bg-muted/80 px-3 py-2 text-left font-semibold whitespace-nowrap w-8">
                      Fill%
                    </th>
                    {visibleCols.map((c) => (
                      <th
                        key={c.key}
                        className="px-2 py-2 text-left font-semibold cursor-pointer whitespace-nowrap"
                        onClick={() => handleSort(c.key)}
                      >
                        <span className="flex items-center gap-1">
                          {c.label}
                          {c.description && (
                            <Tooltip>
                              <TooltipTrigger asChild>
                                <Info className="h-3 w-3 text-muted-foreground" />
                              </TooltipTrigger>
                              <TooltipContent><p className="text-xs">{c.description}</p></TooltipContent>
                            </Tooltip>
                          )}
                          <SortIcon col={c.key} />
                        </span>
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {sorted.map((row: any, i: number) => {
                    const pct = completionRate(row);
                    return (
                      <tr key={row.compoundId ?? i} className="border-b hover:bg-muted/20 transition-colors">
                        <td className="sticky left-0 z-10 bg-background px-2 py-1.5 text-center border-r">
                          <div className="flex flex-col items-center gap-0.5">
                            {pct === 100 ? (
                              <CheckCircle2 className="h-3.5 w-3.5 text-green-500" />
                            ) : pct > 0 ? (
                              <div className="relative h-3.5 w-3.5">
                                <Circle className="h-3.5 w-3.5 text-muted-foreground/30" />
                                <span className="absolute inset-0 flex items-center justify-center text-[7px] font-bold text-yellow-600">{pct}</span>
                              </div>
                            ) : (
                              <Circle className="h-3.5 w-3.5 text-muted-foreground/20" />
                            )}
                          </div>
                        </td>
                        {visibleCols.map((c) => {
                          const v = row[c.key];
                          const formatted = c.format ? c.format(v) : (v ?? "—");
                          const isEmpty = v == null || v === "";
                          return (
                            <td
                              key={c.key}
                              className={`px-2 py-1.5 whitespace-nowrap font-mono ${
                                isEmpty ? "text-muted-foreground/40" : "text-foreground"
                              } ${c.key === "smiles" ? "max-w-[160px] truncate" : ""}`}
                              title={c.key === "smiles" ? String(v ?? "") : undefined}
                            >
                              {formatted === "—" ? <span className="text-muted-foreground/30">—</span> : String(formatted)}
                            </td>
                          );
                        })}
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          )}

          {sorted.length > 0 && (
            <p className="text-xs text-muted-foreground">
              Showing {sorted.length} of {rows?.length ?? 0} compounds · {visibleCols.length} columns visible
            </p>
          )}
        </div>
      </DashboardLayout>
    </TooltipProvider>
  );
}
