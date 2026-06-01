import { useState, useMemo } from "react";
import { trpc } from "@/lib/trpc";
import DashboardLayout from "@/components/DashboardLayout";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Loader2, Download, ArrowUpDown, ArrowUp, ArrowDown, Thermometer, Search } from "lucide-react";
import { toast } from "sonner";

// Property definitions: key, label, direction (higher = better or lower = better), thresholds
const PROPERTIES: Array<{
  key: string;
  label: string;
  unit?: string;
  higherIsBetter: boolean;
  greenThreshold: number;
  redThreshold: number;
  format?: (v: number) => string;
}> = [
  { key: "bbbPermeability", label: "BBB", unit: "", higherIsBetter: true, greenThreshold: 0.7, redThreshold: 0.3 },
  { key: "oralBioavailability", label: "Oral BA", unit: "", higherIsBetter: true, greenThreshold: 0.7, redThreshold: 0.3 },
  { key: "hia", label: "HIA", unit: "", higherIsBetter: true, greenThreshold: 0.8, redThreshold: 0.5 },
  { key: "caco2", label: "Caco-2", unit: "", higherIsBetter: true, greenThreshold: 0.7, redThreshold: 0.3 },
  { key: "pgp", label: "P-gp", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.7 },
  { key: "ppbr", label: "PPBR", unit: "%", higherIsBetter: false, greenThreshold: 70, redThreshold: 90, format: (v) => v.toFixed(1) },
  { key: "halfLife", label: "t½", unit: "h", higherIsBetter: true, greenThreshold: 4, redThreshold: 1, format: (v) => v.toFixed(1) },
  { key: "clearanceHepatocyte", label: "Clearance", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.7 },
  { key: "herg", label: "hERG", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.7 },
  { key: "ames", label: "AMES", unit: "", higherIsBetter: false, greenThreshold: 0.2, redThreshold: 0.5 },
  { key: "dili", label: "DILI", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "ld50", label: "LD50", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "clintox", label: "ClinTox", unit: "", higherIsBetter: false, greenThreshold: 0.2, redThreshold: 0.5 },
  { key: "cyp1a2", label: "CYP1A2", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "cyp2c9", label: "CYP2C9", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "cyp2c19", label: "CYP2C19", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "cyp2d6", label: "CYP2D6", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "cyp3a4", label: "CYP3A4", unit: "", higherIsBetter: false, greenThreshold: 0.3, redThreshold: 0.6 },
  { key: "solubility", label: "Solubility", unit: "", higherIsBetter: true, greenThreshold: 0.6, redThreshold: 0.3 },
  { key: "lipophilicity", label: "LogP", unit: "", higherIsBetter: false, greenThreshold: 3, redThreshold: 5, format: (v) => v.toFixed(2) },
  { key: "tpsa", label: "TPSA", unit: "Å²", higherIsBetter: false, greenThreshold: 90, redThreshold: 140, format: (v) => v.toFixed(1) },
  { key: "qed", label: "QED", unit: "", higherIsBetter: true, greenThreshold: 0.6, redThreshold: 0.3, format: (v) => v.toFixed(2) },
];

function getCellColor(prop: typeof PROPERTIES[0], value: number): string {
  const { higherIsBetter, greenThreshold, redThreshold } = prop;
  if (higherIsBetter) {
    if (value >= greenThreshold) return "bg-green-100 dark:bg-green-900/30 text-green-800 dark:text-green-200";
    if (value <= redThreshold) return "bg-red-100 dark:bg-red-900/30 text-red-800 dark:text-red-200";
    return "bg-yellow-100 dark:bg-yellow-900/30 text-yellow-800 dark:text-yellow-200";
  } else {
    if (value <= greenThreshold) return "bg-green-100 dark:bg-green-900/30 text-green-800 dark:text-green-200";
    if (value >= redThreshold) return "bg-red-100 dark:bg-red-900/30 text-red-800 dark:text-red-200";
    return "bg-yellow-100 dark:bg-yellow-900/30 text-yellow-800 dark:text-yellow-200";
  }
}

type SortDir = "asc" | "desc";

export default function AdmetHeatmap() {
  const [sortKey, setSortKey] = useState<string>("compoundId");
  const [sortDir, setSortDir] = useState<SortDir>("asc");
  const [search, setSearch] = useState("");
  const [selectedProps, setSelectedProps] = useState<Set<string>>(
    new Set(["bbbPermeability", "herg", "ames", "dili", "oralBioavailability", "qed", "tpsa", "lipophilicity", "cyp3a4", "cyp2d6"])
  );

  const { data: rows, isLoading } = trpc.analog.getAllAdmetHeatmap.useQuery();

  const visibleProps = PROPERTIES.filter((p) => selectedProps.has(p.key));

  const filtered = useMemo(() => {
    if (!rows) return [];
    const q = search.toLowerCase();
    return rows.filter((r: any) =>
      !q ||
      (r.compoundId ?? "").toLowerCase().includes(q) ||
      (r.compoundName ?? "").toLowerCase().includes(q) ||
      (r.smiles ?? "").toLowerCase().includes(q)
    );
  }, [rows, search]);

  const sorted = useMemo(() => {
    return [...filtered].sort((a: any, b: any) => {
      let av: any = a[sortKey];
      let bv: any = b[sortKey];
      if (av == null) av = sortDir === "asc" ? Infinity : -Infinity;
      if (bv == null) bv = sortDir === "asc" ? Infinity : -Infinity;
      const numA = parseFloat(av);
      const numB = parseFloat(bv);
      if (!isNaN(numA) && !isNaN(numB)) {
        return sortDir === "asc" ? numA - numB : numB - numA;
      }
      return sortDir === "asc"
        ? String(av).localeCompare(String(bv))
        : String(bv).localeCompare(String(av));
    });
  }, [filtered, sortKey, sortDir]);

  const handleSort = (key: string) => {
    if (sortKey === key) {
      setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    } else {
      setSortKey(key);
      setSortDir("asc");
    }
  };

  const SortIcon = ({ col }: { col: string }) => {
    if (sortKey !== col) return <ArrowUpDown className="h-3 w-3 opacity-40" />;
    return sortDir === "asc" ? <ArrowUp className="h-3 w-3" /> : <ArrowDown className="h-3 w-3" />;
  };

  const exportCSV = () => {
    if (!sorted.length) return;
    const headers = ["Compound ID", "Name", "SMILES", ...visibleProps.map((p) => `${p.label}${p.unit ? ` (${p.unit})` : ""}`)];
    const csvRows = sorted.map((r: any) => [
      r.compoundId ?? "",
      r.compoundName ?? "",
      r.smiles ?? "",
      ...visibleProps.map((p) => {
        const v = parseFloat(r[p.key]);
        return isNaN(v) ? "" : (p.format ? p.format(v) : v.toFixed(3));
      }),
    ]);
    const csv = [headers, ...csvRows].map((row) => row.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `admet_heatmap_${new Date().toISOString().slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(url);
    toast.success("CSV exported");
  };

  const toggleProp = (key: string) => {
    setSelectedProps((prev) => {
      const next = new Set(prev);
      if (next.has(key)) {
        if (next.size > 1) next.delete(key);
      } else {
        next.add(key);
      }
      return next;
    });
  };

  return (
    <DashboardLayout>
      <div className="space-y-4 p-4">
        {/* Header */}
        <div className="flex items-start justify-between gap-4 flex-wrap">
          <div>
            <h1 className="text-2xl font-bold flex items-center gap-2">
              <Thermometer className="h-6 w-6 text-orange-500" />
              ADMET-AI Heatmap
            </h1>
            <p className="text-sm text-muted-foreground mt-1">
              Chemprop-predicted ADMET properties for all analogs. Color: <span className="text-green-600 font-medium">green = favorable</span>, <span className="text-yellow-600 font-medium">yellow = borderline</span>, <span className="text-red-600 font-medium">red = concern</span>.
            </p>
          </div>
          <Button variant="outline" size="sm" onClick={exportCSV} disabled={!sorted.length}>
            <Download className="h-4 w-4 mr-1" /> Export CSV
          </Button>
        </div>

        {/* Property selector */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm">Visible Properties</CardTitle>
            <CardDescription className="text-xs">Toggle columns to show/hide</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="flex flex-wrap gap-1.5">
              {PROPERTIES.map((p) => (
                <button
                  key={p.key}
                  onClick={() => toggleProp(p.key)}
                  className={`px-2 py-0.5 rounded text-xs font-medium border transition-colors ${
                    selectedProps.has(p.key)
                      ? "bg-primary text-primary-foreground border-primary"
                      : "bg-transparent text-muted-foreground border-border hover:border-primary/50"
                  }`}
                >
                  {p.label}
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

        {/* Heatmap table */}
        {isLoading ? (
          <div className="flex items-center justify-center py-16">
            <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
            <span className="ml-2 text-muted-foreground">Loading ADMET data...</span>
          </div>
        ) : sorted.length === 0 ? (
          <Card>
            <CardContent className="py-12 text-center text-muted-foreground">
              <Thermometer className="h-10 w-10 mx-auto mb-3 opacity-30" />
              <p className="font-medium">No ADMET data yet</p>
              <p className="text-sm mt-1">Run ADMET analysis on analogs from the Compound Testing page to populate this heatmap.</p>
            </CardContent>
          </Card>
        ) : (
          <div className="overflow-x-auto rounded-lg border">
            <table className="w-full text-xs border-collapse">
              <thead>
                <tr className="bg-muted/50 border-b">
                  <th
                    className="sticky left-0 z-10 bg-muted/80 px-3 py-2 text-left font-semibold cursor-pointer whitespace-nowrap"
                    onClick={() => handleSort("compoundId")}
                  >
                    <span className="flex items-center gap-1">Compound <SortIcon col="compoundId" /></span>
                  </th>
                  {visibleProps.map((p) => (
                    <th
                      key={p.key}
                      className="px-2 py-2 text-center font-semibold cursor-pointer whitespace-nowrap"
                      onClick={() => handleSort(p.key)}
                    >
                      <span className="flex items-center justify-center gap-1">
                        {p.label}
                        {p.unit && <span className="text-muted-foreground font-normal">({p.unit})</span>}
                        <SortIcon col={p.key} />
                      </span>
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {sorted.map((row: any, i: number) => (
                  <tr key={row.id ?? i} className="border-b hover:bg-muted/20 transition-colors">
                    <td className="sticky left-0 z-10 bg-background px-3 py-1.5 whitespace-nowrap border-r">
                      <div className="font-mono font-medium">{row.compoundId ?? `#${row.analogId}`}</div>
                      {row.compoundName && (
                        <div className="text-muted-foreground text-xs truncate max-w-[140px]">{row.compoundName}</div>
                      )}
                    </td>
                    {visibleProps.map((p) => {
                      const raw = row[p.key];
                      const v = parseFloat(raw);
                      if (isNaN(v)) {
                        return (
                          <td key={p.key} className="px-2 py-1.5 text-center text-muted-foreground">
                            —
                          </td>
                        );
                      }
                      const colorClass = getCellColor(p, v);
                      const display = p.format ? p.format(v) : v.toFixed(3);
                      return (
                        <td key={p.key} className={`px-2 py-1.5 text-center font-mono font-medium rounded-sm ${colorClass}`}>
                          {display}
                        </td>
                      );
                    })}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}

        {sorted.length > 0 && (
          <p className="text-xs text-muted-foreground">
            Showing {sorted.length} of {rows?.length ?? 0} records
          </p>
        )}
      </div>
    </DashboardLayout>
  );
}
