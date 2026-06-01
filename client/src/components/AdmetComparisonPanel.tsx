import { useState, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Loader2, X, Plus, Download, FileText, Search, GitCompare } from "lucide-react";
import { toast } from "sonner";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";

// ─── Colour palette for up to 4 compounds ────────────────────────────────────
const COLORS = ["#6366f1", "#10b981", "#f59e0b", "#ef4444"];
const COLOR_BG = ["bg-indigo-50", "bg-emerald-50", "bg-amber-50", "bg-red-50"];
const COLOR_BORDER = [
  "border-indigo-300",
  "border-emerald-300",
  "border-amber-300",
  "border-red-300",
];
const COLOR_TEXT = [
  "text-indigo-700",
  "text-emerald-700",
  "text-amber-700",
  "text-red-700",
];

// ─── ADMET property definitions ──────────────────────────────────────────────
interface AdmetProp {
  key: string;
  label: string;
  group: string;
  unit?: string;
  /** lower is better (true) vs higher is better (false) */
  lowerIsBetter?: boolean;
  /** extract numeric value from the parsed result JSON */
  extract: (data: any) => number | string | null;
}

const ADMET_PROPS: AdmetProp[] = [
  // Toxicity profile
  {
    key: "herg_risk",
    label: "hERG Risk Score",
    group: "Toxicity",
    lowerIsBetter: true,
    extract: (d) => d?.toxicity_profile?.hERG?.risk_score ?? null,
  },
  {
    key: "herg_level",
    label: "hERG Risk Level",
    group: "Toxicity",
    extract: (d) => d?.toxicity_profile?.hERG?.risk_level ?? null,
  },
  {
    key: "hepato_risk",
    label: "Hepatotoxicity Score",
    group: "Toxicity",
    lowerIsBetter: true,
    extract: (d) => d?.toxicity_profile?.hepatotoxicity?.risk_score ?? null,
  },
  {
    key: "hepato_level",
    label: "Hepatotoxicity Level",
    group: "Toxicity",
    extract: (d) => d?.toxicity_profile?.hepatotoxicity?.risk_level ?? null,
  },
  {
    key: "mutagenicity",
    label: "Mutagenicity",
    group: "Toxicity",
    extract: (d) => d?.toxicity_profile?.mutagenicity?.prediction ?? null,
  },
  {
    key: "carcinogenicity",
    label: "Carcinogenicity",
    group: "Toxicity",
    extract: (d) => d?.toxicity_profile?.carcinogenicity?.prediction ?? null,
  },
  // Synthetic accessibility
  {
    key: "sa_score",
    label: "SA Score",
    group: "Synthesis",
    lowerIsBetter: true,
    extract: (d) => d?.synthetic_accessibility?.sa_score ?? null,
  },
  {
    key: "sa_difficulty",
    label: "Synthesis Difficulty",
    group: "Synthesis",
    extract: (d) => d?.synthetic_accessibility?.difficulty ?? null,
  },
  // Flat ML properties (from ADMET-AI / Chemprop)
  {
    key: "bbb",
    label: "BBB Permeability",
    group: "ADMET",
    lowerIsBetter: false,
    extract: (d) => d?.bbb ?? d?.BBB ?? null,
  },
  {
    key: "caco2",
    label: "Caco-2 Permeability",
    group: "ADMET",
    unit: "nm/s",
    lowerIsBetter: false,
    extract: (d) => d?.caco2 ?? d?.Caco2 ?? null,
  },
  {
    key: "ames",
    label: "AMES Mutagenicity",
    group: "ADMET",
    lowerIsBetter: true,
    extract: (d) => d?.ames ?? d?.AMES ?? null,
  },
  {
    key: "dili",
    label: "DILI (Liver Injury)",
    group: "ADMET",
    lowerIsBetter: true,
    extract: (d) => d?.dili ?? d?.DILI ?? null,
  },
  {
    key: "half_life",
    label: "Half-life",
    group: "ADMET",
    unit: "h",
    lowerIsBetter: false,
    extract: (d) => d?.half_life ?? d?.HalfLife ?? null,
  },
  {
    key: "clearance",
    label: "Clearance",
    group: "ADMET",
    unit: "mL/min/kg",
    lowerIsBetter: true,
    extract: (d) => d?.clearance ?? d?.Clearance ?? null,
  },
  {
    key: "solubility",
    label: "Solubility",
    group: "ADMET",
    unit: "log mol/L",
    lowerIsBetter: false,
    extract: (d) => d?.solubility ?? d?.Solubility ?? null,
  },
  {
    key: "bioavailability",
    label: "Oral Bioavailability",
    group: "ADMET",
    unit: "%",
    lowerIsBetter: false,
    extract: (d) => d?.bioavailability ?? d?.Bioavailability ?? null,
  },
];

const GROUPS = Array.from(new Set(ADMET_PROPS.map((p) => p.group)));

// ─── Helpers ─────────────────────────────────────────────────────────────────
function parseLatestAdmet(testResults: any[]): any | null {
  if (!testResults) return null;
  const admet = testResults
    .filter((r: any) => r.testType === "admet" && r.testStatus === "completed")
    .sort(
      (a: any, b: any) =>
        new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime()
    );
  if (admet.length === 0) return null;
  try {
    return JSON.parse(admet[0].results || "{}");
  } catch {
    return null;
  }
}

function numericValues(
  prop: AdmetProp,
  compounds: CompoundWithData[]
): (number | null)[] {
  return compounds.map((c) => {
    const v = prop.extract(c.admetData);
    if (v === null || v === undefined) return null;
    const n = typeof v === "number" ? v : parseFloat(String(v));
    return isNaN(n) ? null : n;
  });
}

function cellClass(
  prop: AdmetProp,
  value: number | null,
  allValues: (number | null)[]
): string {
  if (value === null) return "";
  const nums = allValues.filter((v) => v !== null) as number[];
  if (nums.length < 2) return "";
  const best = prop.lowerIsBetter ? Math.min(...nums) : Math.max(...nums);
  const worst = prop.lowerIsBetter ? Math.max(...nums) : Math.min(...nums);
  if (value === best) return "bg-green-50 font-semibold text-green-800";
  if (value === worst) return "bg-red-50 text-red-700";
  return "";
}

// ─── Types ────────────────────────────────────────────────────────────────────
interface CompoundMeta {
  id: number;
  compoundName: string;
  parentCompound: string;
  smiles: string;
  confidenceScore: number;
}

interface CompoundWithData extends CompoundMeta {
  admetData: any | null;
  analyzedAt: string | null;
  colorIdx: number;
}

// ─── Export helpers ───────────────────────────────────────────────────────────
function exportCombinedCsv(compounds: CompoundWithData[]) {
  const header = ["Property", "Group", ...compounds.map((c) => c.compoundName)];
  const rows = ADMET_PROPS.map((prop) => {
    const vals = compounds.map((c) => {
      const v = prop.extract(c.admetData);
      return v === null || v === undefined ? "" : String(v);
    });
    return [prop.label, prop.group, ...vals];
  });
  const csv = [header, ...rows]
    .map((r) => r.map((v) => `"${v.replace(/"/g, '""')}"`).join(","))
    .join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `admet-comparison-${new Date().toISOString().slice(0, 10)}.csv`;
  a.click();
  URL.revokeObjectURL(url);
}

function exportCombinedPdf(compounds: CompoundWithData[]) {
  const colWidth = Math.floor(60 / compounds.length);
  const compoundHeaders = compounds
    .map(
      (c, i) =>
        `<th style="background:${COLORS[i]};color:#fff;padding:8px 12px;min-width:${colWidth}%">${c.compoundName}</th>`
    )
    .join("");

  const groupedProps: Record<string, AdmetProp[]> = {};
  ADMET_PROPS.forEach((p) => {
    if (!groupedProps[p.group]) groupedProps[p.group] = [];
    groupedProps[p.group].push(p);
  });

  const tableRows = GROUPS.flatMap((group) => {
    const props = groupedProps[group] || [];
    const headerRow = `<tr><td colspan="${compounds.length + 1}" style="background:#f3f4f6;font-weight:700;padding:6px 12px;font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:#374151">${group}</td></tr>`;
    const propRows = props.map((prop) => {
      const allNums = numericValues(prop, compounds);
      const cells = compounds.map((c, i) => {
        const v = prop.extract(c.admetData);
        const display =
          v === null || v === undefined
            ? '<span style="color:#9ca3af">—</span>'
            : `${typeof v === "number" ? v.toFixed(2) : v}${prop.unit ? ` <small>${prop.unit}</small>` : ""}`;
        const n = allNums[i];
        const nums = allNums.filter((x) => x !== null) as number[];
        let bg = "";
        if (n !== null && nums.length >= 2) {
          const best = prop.lowerIsBetter ? Math.min(...nums) : Math.max(...nums);
          const worst = prop.lowerIsBetter ? Math.max(...nums) : Math.min(...nums);
          if (n === best) bg = "background:#d1fae5";
          else if (n === worst) bg = "background:#fee2e2";
        }
        return `<td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;${bg}">${display}</td>`;
      });
      return `<tr><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;color:#374151">${prop.label}</td>${cells.join("")}</tr>`;
    });
    return [headerRow, ...propRows];
  }).join("");

  const html = `<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>ADMET Comparison Report</title>
  <style>
    body{font-family:system-ui,sans-serif;margin:32px;color:#111}
    h1{font-size:22px;margin-bottom:4px}
    p{color:#6b7280;font-size:13px;margin-bottom:24px}
    table{border-collapse:collapse;width:100%;font-size:13px}
    th{text-align:left}
    .legend{display:flex;gap:16px;flex-wrap:wrap;margin-bottom:16px;font-size:12px}
    .legend-item{display:flex;align-items:center;gap:6px}
    .dot{width:12px;height:12px;border-radius:50%}
    .note{font-size:11px;color:#9ca3af;margin-top:16px}
  </style>
</head>
<body>
  <h1>ADMET Comparison Report</h1>
  <p>Generated ${new Date().toLocaleString()} &nbsp;|&nbsp; ${compounds.length} compounds</p>
  <div class="legend">
    ${compounds.map((c, i) => `<div class="legend-item"><div class="dot" style="background:${COLORS[i]}"></div><span>${c.compoundName}</span></div>`).join("")}
  </div>
  <table>
    <thead>
      <tr>
        <th style="padding:8px 12px;background:#1e293b;color:#fff;width:40%">Property</th>
        ${compoundHeaders}
      </tr>
    </thead>
    <tbody>${tableRows}</tbody>
  </table>
  <p class="note">Green = best value &nbsp;|&nbsp; Red = worst value &nbsp;|&nbsp; — = not yet analyzed</p>
</body>
</html>`;

  const w = window.open("", "_blank");
  if (w) {
    w.document.write(html);
    w.document.close();
    w.print();
  }
}

// ─── Per-compound data hook ───────────────────────────────────────────────────
function useCompoundAdmet(id: number | null) {
  return trpc.analog.getTestResults.useQuery(
    { analogId: id! },
    { enabled: !!id }
  );
}

// ─── Main component ───────────────────────────────────────────────────────────
export function AdmetComparisonPanel() {
  const [selectedIds, setSelectedIds] = useState<number[]>([]);
  const [searchQuery, setSearchQuery] = useState("");
  const [showPicker, setShowPicker] = useState(true);

  const { data: allAnalogs, isLoading: analogsLoading } = trpc.analog.list.useQuery(
    { limit: 200, offset: 0 },
    { staleTime: 60_000 }
  );

  // Fetch test results for each selected compound (up to 4 queries)
  const r0 = useCompoundAdmet(selectedIds[0] ?? null);
  const r1 = useCompoundAdmet(selectedIds[1] ?? null);
  const r2 = useCompoundAdmet(selectedIds[2] ?? null);
  const r3 = useCompoundAdmet(selectedIds[3] ?? null);
  const resultQueries = [r0, r1, r2, r3];

  const filteredAnalogs = useMemo(() => {
    if (!allAnalogs) return [];
    return (allAnalogs as CompoundMeta[]).filter(
      (a) =>
        a.compoundName &&
        (searchQuery === "" ||
          a.compoundName.toLowerCase().includes(searchQuery.toLowerCase()) ||
          a.parentCompound?.toLowerCase().includes(searchQuery.toLowerCase()))
    );
  }, [allAnalogs, searchQuery]);

  const compoundsWithData: CompoundWithData[] = useMemo(() => {
    return selectedIds
      .map((id, idx) => {
        const meta = (allAnalogs as CompoundMeta[] | undefined)?.find(
          (a) => a.id === id
        );
        if (!meta) return null;
        const results = resultQueries[idx]?.data ?? null;
        const admetData = results ? parseLatestAdmet(results as any[]) : null;
        const latestAdmet = (results as any[] | null)
          ?.filter((r: any) => r.testType === "admet" && r.testStatus === "completed")
          .sort(
            (a: any, b: any) =>
              new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime()
          )[0];
        return {
          ...meta,
          admetData,
          analyzedAt: latestAdmet?.createdAt
            ? new Date(latestAdmet.createdAt).toLocaleString()
            : null,
          colorIdx: idx,
        } as CompoundWithData;
      })
      .filter(Boolean) as CompoundWithData[];
  }, [selectedIds, allAnalogs, r0.data, r1.data, r2.data, r3.data]);

  const isAnyLoading = resultQueries
    .slice(0, selectedIds.length)
    .some((q) => q.isLoading);

  const toggleSelect = (id: number) => {
    setSelectedIds((prev) => {
      if (prev.includes(id)) return prev.filter((x) => x !== id);
      if (prev.length >= 4) {
        toast.warning("Maximum 4 compounds can be compared at once");
        return prev;
      }
      return [...prev, id];
    });
  };

  const removeSelected = (id: number) => {
    setSelectedIds((prev) => prev.filter((x) => x !== id));
  };

  const hasData = compoundsWithData.some((c) => c.admetData !== null);

  return (
    <div className="space-y-4">
      {/* Header bar */}
      <div className="flex items-center justify-between gap-4 flex-wrap">
        <div className="flex items-center gap-2">
          <GitCompare className="h-5 w-5 text-indigo-500" />
          <h3 className="font-semibold text-base">ADMET Side-by-Side Comparison</h3>
          <Badge variant="secondary">{selectedIds.length}/4 selected</Badge>
        </div>
        <div className="flex items-center gap-2">
          <Button
            variant="outline"
            size="sm"
            onClick={() => setShowPicker((v) => !v)}
          >
            {showPicker ? "Hide" : "Show"} Compound Picker
          </Button>
          {hasData && (
            <DropdownMenu>
              <DropdownMenuTrigger asChild>
                <Button size="sm" variant="outline" disabled={!hasData}>
                  <Download className="h-4 w-4 mr-1" />
                  Export Combined
                </Button>
              </DropdownMenuTrigger>
              <DropdownMenuContent align="end">
                <DropdownMenuItem
                  onClick={() => exportCombinedCsv(compoundsWithData)}
                >
                  <Download className="h-4 w-4 mr-2" />
                  Download CSV
                </DropdownMenuItem>
                <DropdownMenuItem
                  onClick={() => exportCombinedPdf(compoundsWithData)}
                >
                  <FileText className="h-4 w-4 mr-2" />
                  Export as PDF
                </DropdownMenuItem>
              </DropdownMenuContent>
            </DropdownMenu>
          )}
        </div>
      </div>

      {/* Compound picker */}
      {showPicker && (
        <Card>
          <CardHeader className="pb-3">
            <CardTitle className="text-sm font-medium">Select Compounds (2–4)</CardTitle>
            <CardDescription className="text-xs">
              Choose compounds that have already had ADMET analysis run. Green = has ADMET data.
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-3">
            {/* Selected chips */}
            {selectedIds.length > 0 && (
              <div className="flex flex-wrap gap-2">
                {compoundsWithData.map((c) => (
                  <div
                    key={c.id}
                    className={`flex items-center gap-1.5 px-3 py-1.5 rounded-full border text-sm font-medium ${COLOR_BG[c.colorIdx]} ${COLOR_BORDER[c.colorIdx]} ${COLOR_TEXT[c.colorIdx]}`}
                  >
                    <span
                      className="w-2 h-2 rounded-full"
                      style={{ background: COLORS[c.colorIdx] }}
                    />
                    {c.compoundName}
                    {c.admetData ? (
                      <span className="text-xs opacity-70">✓</span>
                    ) : (
                      <span className="text-xs opacity-50">no data</span>
                    )}
                    <button
                      onClick={() => removeSelected(c.id)}
                      className="ml-1 hover:opacity-70"
                    >
                      <X className="h-3 w-3" />
                    </button>
                  </div>
                ))}
              </div>
            )}

            {/* Search */}
            <div className="relative">
              <Search className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground" />
              <Input
                placeholder="Search compounds..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                className="pl-8"
              />
            </div>

            {/* Compound list */}
            {analogsLoading ? (
              <div className="flex items-center justify-center py-6">
                <Loader2 className="h-5 w-5 animate-spin text-muted-foreground" />
              </div>
            ) : (
              <div className="max-h-48 overflow-y-auto space-y-1 pr-1">
                {filteredAnalogs.length === 0 && (
                  <p className="text-sm text-muted-foreground text-center py-4">
                    No compounds found
                  </p>
                )}
                {filteredAnalogs.map((analog) => {
                  const isSelected = selectedIds.includes(analog.id);
                  const selIdx = selectedIds.indexOf(analog.id);
                  const isDisabled = !isSelected && selectedIds.length >= 4;
                  return (
                    <button
                      key={analog.id}
                      onClick={() => toggleSelect(analog.id)}
                      disabled={isDisabled}
                      className={`w-full flex items-center justify-between px-3 py-2 rounded-md text-sm transition-colors text-left ${
                        isSelected
                          ? `${COLOR_BG[selIdx]} ${COLOR_BORDER[selIdx]} border`
                          : "hover:bg-muted border border-transparent"
                      } ${isDisabled ? "opacity-40 cursor-not-allowed" : "cursor-pointer"}`}
                    >
                      <div className="flex items-center gap-2 min-w-0">
                        {isSelected && (
                          <span
                            className="w-2.5 h-2.5 rounded-full flex-shrink-0"
                            style={{ background: COLORS[selIdx] }}
                          />
                        )}
                        {!isSelected && (
                          <Plus className="h-3.5 w-3.5 text-muted-foreground flex-shrink-0" />
                        )}
                        <span className="truncate font-medium">{analog.compoundName}</span>
                        <span className="text-muted-foreground text-xs truncate hidden sm:block">
                          {analog.parentCompound}
                        </span>
                      </div>
                      <div className="flex items-center gap-1.5 flex-shrink-0">
                        <Badge variant="secondary" className="text-xs">
                          {analog.confidenceScore}%
                        </Badge>
                      </div>
                    </button>
                  );
                })}
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Loading state */}
      {isAnyLoading && selectedIds.length > 0 && (
        <div className="flex items-center gap-2 text-sm text-muted-foreground py-2">
          <Loader2 className="h-4 w-4 animate-spin" />
          Loading ADMET data for selected compounds…
        </div>
      )}

      {/* Empty state */}
      {selectedIds.length < 2 && !isAnyLoading && (
        <div className="text-center py-10 text-muted-foreground border-2 border-dashed rounded-lg">
          <GitCompare className="h-8 w-8 mx-auto mb-2 opacity-30" />
          <p className="text-sm">Select at least 2 compounds above to compare their ADMET properties.</p>
        </div>
      )}

      {/* Comparison table */}
      {selectedIds.length >= 2 && !isAnyLoading && (
        <div className="overflow-x-auto rounded-lg border">
          <table className="w-full text-sm border-collapse">
            <thead>
              <tr>
                <th className="text-left px-4 py-3 bg-slate-800 text-white font-medium w-[200px] sticky left-0 z-10">
                  Property
                </th>
                {compoundsWithData.map((c) => (
                  <th
                    key={c.id}
                    className="px-4 py-3 text-white font-medium text-left min-w-[160px]"
                    style={{ background: COLORS[c.colorIdx] }}
                  >
                    <div className="font-semibold">{c.compoundName}</div>
                    {c.analyzedAt ? (
                      <div className="text-xs opacity-80 font-normal mt-0.5">
                        {c.analyzedAt}
                      </div>
                    ) : (
                      <div className="text-xs opacity-60 font-normal mt-0.5">
                        No ADMET data yet
                      </div>
                    )}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {GROUPS.map((group) => {
                const props = ADMET_PROPS.filter((p) => p.group === group);
                return [
                  // Group header row
                  <tr key={`group-${group}`}>
                    <td
                      colSpan={compoundsWithData.length + 1}
                      className="px-4 py-2 bg-gray-100 text-xs font-bold uppercase tracking-wide text-gray-500"
                    >
                      {group}
                    </td>
                  </tr>,
                  // Property rows
                  ...props.map((prop) => {
                    const allNums = numericValues(prop, compoundsWithData);
                    return (
                      <tr
                        key={prop.key}
                        className="border-b border-gray-100 hover:bg-gray-50 transition-colors"
                      >
                        <td className="px-4 py-2.5 font-medium text-gray-700 sticky left-0 bg-white border-r border-gray-100">
                          {prop.label}
                          {prop.unit && (
                            <span className="text-xs text-muted-foreground ml-1">
                              ({prop.unit})
                            </span>
                          )}
                        </td>
                        {compoundsWithData.map((c, ci) => {
                          const v = prop.extract(c.admetData);
                          const n = allNums[ci];
                          const cls = cellClass(prop, n, allNums);
                          const display =
                            v === null || v === undefined ? (
                              <span className="text-gray-300">—</span>
                            ) : typeof v === "number" ? (
                              v.toFixed(2)
                            ) : (
                              String(v)
                            );
                          return (
                            <td
                              key={c.id}
                              className={`px-4 py-2.5 ${cls}`}
                            >
                              {display}
                            </td>
                          );
                        })}
                      </tr>
                    );
                  }),
                ];
              })}
            </tbody>
          </table>
          <div className="px-4 py-2 bg-gray-50 text-xs text-muted-foreground border-t">
            <span className="inline-block w-3 h-3 rounded bg-green-100 border border-green-300 mr-1 align-middle" />
            Best value &nbsp;
            <span className="inline-block w-3 h-3 rounded bg-red-100 border border-red-300 mr-1 align-middle" />
            Worst value &nbsp;|&nbsp; Highlighting applies only when ≥2 compounds have numeric data for that property.
          </div>
        </div>
      )}
    </div>
  );
}
