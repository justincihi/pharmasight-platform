import { useState, useMemo } from "react";
import {
  RadarChart,
  Radar,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  ResponsiveContainer,
  Legend,
  Tooltip,
} from "recharts";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { trpc } from "@/lib/trpc";
import { X, Plus, BarChart3, Info } from "lucide-react";

// Colours for up to 4 compounds
const COLORS = ["#6366f1", "#10b981", "#f59e0b", "#ef4444"];

interface Compound {
  id: number;
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  confidenceScore: number;
  safetyScore: number;
  efficacyScore: number;
  drugLikenessScore: number;
  similarityScore: number;
  patentStatus: string;
}

interface RadarDataPoint {
  metric: string;
  fullMark: number;
  [key: string]: number | string;
}

const METRICS: { key: keyof Compound; label: string; fullMark: number }[] = [
  { key: "confidenceScore", label: "Confidence", fullMark: 100 },
  { key: "safetyScore", label: "Safety", fullMark: 100 },
  { key: "efficacyScore", label: "Efficacy", fullMark: 100 },
  { key: "drugLikenessScore", label: "Drug-likeness", fullMark: 100 },
  { key: "similarityScore", label: "Similarity", fullMark: 100 },
];

function buildRadarData(selected: Compound[]): RadarDataPoint[] {
  return METRICS.map(({ key, label, fullMark }) => {
    const point: RadarDataPoint = { metric: label, fullMark };
    selected.forEach((c) => {
      point[c.compoundName] = (c[key] as number) ?? 0;
    });
    return point;
  });
}

export default function CompoundComparisonChart() {
  const [selectedIds, setSelectedIds] = useState<number[]>([]);
  const [searchQuery, setSearchQuery] = useState("");

  const { data: allAnalogs, isLoading } = trpc.analog.list.useQuery(
    { limit: 200, offset: 0 },
    { staleTime: 60_000 }
  );

  const analogs: Compound[] = useMemo(() => {
    if (!allAnalogs) return [];
    return (allAnalogs as Compound[]).filter(
      (a) =>
        a.compoundName &&
        (searchQuery === "" ||
          a.compoundName.toLowerCase().includes(searchQuery.toLowerCase()) ||
          a.parentCompound?.toLowerCase().includes(searchQuery.toLowerCase()))
    );
  }, [allAnalogs, searchQuery]);

  const selected: Compound[] = useMemo(
    () =>
      selectedIds
        .map((id) => analogs.find((a) => a.id === id))
        .filter(Boolean) as Compound[],
    [selectedIds, analogs]
  );

  const radarData = useMemo(() => buildRadarData(selected), [selected]);

  const toggleSelect = (id: number) => {
    setSelectedIds((prev) => {
      if (prev.includes(id)) return prev.filter((x) => x !== id);
      if (prev.length >= 4) return prev; // max 4
      return [...prev, id];
    });
  };

  const removeSelected = (id: number) => {
    setSelectedIds((prev) => prev.filter((x) => x !== id));
  };

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-2">
        <BarChart3 className="w-5 h-5 text-purple-500" />
        <h2 className="text-xl font-semibold">Compound Comparison</h2>
        <span className="text-sm text-muted-foreground">(select 2–4 analogs)</span>
      </div>

      {/* Selected chips */}
      {selected.length > 0 && (
        <div className="flex flex-wrap gap-2">
          {selected.map((c, i) => (
            <Badge
              key={c.id}
              style={{ backgroundColor: COLORS[i] + "22", borderColor: COLORS[i] }}
              className="border text-sm px-3 py-1 gap-1"
            >
              <span style={{ color: COLORS[i] }} className="font-semibold">
                {c.compoundName}
              </span>
              <button
                onClick={() => removeSelected(c.id)}
                className="ml-1 hover:opacity-70"
                aria-label="Remove"
              >
                <X className="w-3 h-3" />
              </button>
            </Badge>
          ))}
          {selected.length < 4 && (
            <span className="text-xs text-muted-foreground self-center">
              + {4 - selected.length} more can be added
            </span>
          )}
        </div>
      )}

      {/* Radar Chart */}
      {selected.length >= 2 ? (
        <div className="rounded-xl border border-border bg-card p-4">
          <ResponsiveContainer width="100%" height={380}>
            <RadarChart data={radarData} margin={{ top: 20, right: 30, bottom: 20, left: 30 }}>
              <PolarGrid stroke="hsl(var(--border))" />
              <PolarAngleAxis
                dataKey="metric"
                tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 13 }}
              />
              <PolarRadiusAxis
                angle={90}
                domain={[0, 100]}
                tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 11 }}
              />
              {selected.map((c, i) => (
                <Radar
                  key={c.id}
                  name={c.compoundName}
                  dataKey={c.compoundName}
                  stroke={COLORS[i]}
                  fill={COLORS[i]}
                  fillOpacity={0.18}
                  strokeWidth={2}
                />
              ))}
              <Legend
                wrapperStyle={{ fontSize: 13, paddingTop: 12 }}
                formatter={(value) => (
                  <span style={{ color: "hsl(var(--foreground))" }}>{value}</span>
                )}
              />
              <Tooltip
                contentStyle={{
                  background: "hsl(var(--card))",
                  border: "1px solid hsl(var(--border))",
                  borderRadius: 8,
                  fontSize: 13,
                }}
                formatter={(value: number) => [`${value}%`, undefined]}
              />
            </RadarChart>
          </ResponsiveContainer>

          {/* Score table */}
          <div className="mt-4 overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border">
                  <th className="text-left py-2 pr-4 font-medium text-muted-foreground">Metric</th>
                  {selected.map((c, i) => (
                    <th
                      key={c.id}
                      className="text-center py-2 px-3 font-semibold"
                      style={{ color: COLORS[i] }}
                    >
                      {c.compoundName}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {METRICS.map(({ key, label }) => (
                  <tr key={key} className="border-b border-border/50">
                    <td className="py-2 pr-4 text-muted-foreground">{label}</td>
                    {selected.map((c) => {
                      const val = (c[key] as number) ?? 0;
                      const color =
                        val >= 80 ? "text-green-600" : val >= 60 ? "text-yellow-600" : "text-red-500";
                      return (
                        <td key={c.id} className={`text-center py-2 px-3 font-mono font-semibold ${color}`}>
                          {val}%
                        </td>
                      );
                    })}
                  </tr>
                ))}
                <tr>
                  <td className="py-2 pr-4 text-muted-foreground">Patent Status</td>
                  {selected.map((c) => (
                    <td key={c.id} className="text-center py-2 px-3 text-xs">
                      <Badge
                        variant={c.patentStatus === "patent-free" ? "default" : "secondary"}
                        className="text-xs"
                      >
                        {c.patentStatus}
                      </Badge>
                    </td>
                  ))}
                </tr>
              </tbody>
            </table>
          </div>
        </div>
      ) : (
        <div className="rounded-xl border border-dashed border-border bg-muted/20 p-10 text-center">
          <BarChart3 className="w-10 h-10 text-muted-foreground mx-auto mb-3" />
          <p className="text-muted-foreground font-medium">
            Select at least 2 compounds below to compare
          </p>
          <p className="text-sm text-muted-foreground mt-1">
            Up to 4 analogs can be compared at once
          </p>
        </div>
      )}

      {/* Compound Picker */}
      <div>
        <div className="flex items-center justify-between mb-3">
          <h3 className="font-semibold text-sm">Select Compounds</h3>
          <input
            type="text"
            placeholder="Search by name or parent..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            className="px-3 py-1.5 text-sm border border-border rounded-md bg-background w-56"
          />
        </div>

        {isLoading ? (
          <div className="text-muted-foreground text-sm py-4 text-center">Loading analogs…</div>
        ) : analogs.length === 0 ? (
          <div className="text-muted-foreground text-sm py-4 text-center">No analogs found.</div>
        ) : (
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-2 max-h-72 overflow-y-auto pr-1">
            {analogs.slice(0, 60).map((a) => {
              const isSelected = selectedIds.includes(a.id);
              const selIdx = selectedIds.indexOf(a.id);
              const isDisabled = !isSelected && selectedIds.length >= 4;
              return (
                <button
                  key={a.id}
                  onClick={() => toggleSelect(a.id)}
                  disabled={isDisabled}
                  className={`text-left p-3 rounded-lg border transition-all text-sm ${
                    isSelected
                      ? "border-2"
                      : isDisabled
                      ? "border-border opacity-40 cursor-not-allowed"
                      : "border-border hover:border-blue-300 hover:bg-muted/40"
                  }`}
                  style={
                    isSelected
                      ? { borderColor: COLORS[selIdx], backgroundColor: COLORS[selIdx] + "11" }
                      : {}
                  }
                >
                  <div className="flex items-center justify-between gap-2">
                    <span
                      className="font-semibold truncate"
                      style={isSelected ? { color: COLORS[selIdx] } : {}}
                    >
                      {a.compoundName}
                    </span>
                    {isSelected ? (
                      <X className="w-3.5 h-3.5 flex-shrink-0" style={{ color: COLORS[selIdx] }} />
                    ) : (
                      <Plus className="w-3.5 h-3.5 flex-shrink-0 text-muted-foreground" />
                    )}
                  </div>
                  <p className="text-xs text-muted-foreground truncate">
                    Parent: {a.parentCompound}
                  </p>
                  <div className="flex gap-2 mt-1">
                    <span className="text-xs text-blue-600 font-mono">{a.confidenceScore}% conf</span>
                    <span className="text-xs text-green-600 font-mono">{a.safetyScore}% safe</span>
                  </div>
                </button>
              );
            })}
          </div>
        )}
      </div>
    </div>
  );
}
