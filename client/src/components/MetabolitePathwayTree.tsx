import { useMemo, useState } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Download, FlaskConical, ChevronRight } from "lucide-react";

interface Metabolite {
  id: number;
  analogId: number;
  metaboliteSmiles: string;
  metaboliteName: string | null;
  reactionType: string | null;
  enzyme: string | null;
  phase: number | null;
  molecularWeight: number | null;
  logP: number | null;
  predictedActivity: string | null;
  confidenceScore: number | null;
  createdAt: Date | string;
}

interface MetabolitePathwayTreeProps {
  parentName: string;
  parentSmiles: string;
  metabolites: Metabolite[];
}

const PHASE_COLORS: Record<number, { bg: string; text: string; border: string }> = {
  1: { bg: "bg-blue-500/10", text: "text-blue-600 dark:text-blue-400", border: "border-blue-500/30" },
  2: { bg: "bg-purple-500/10", text: "text-purple-600 dark:text-purple-400", border: "border-purple-500/30" },
};

const REACTION_ICONS: Record<string, string> = {
  "Aromatic Hydroxylation": "→OH",
  "Aliphatic Hydroxylation": "→OH",
  "N-Dealkylation": "→NH",
  "O-Dealkylation": "→OH",
  "N-Oxidation": "→N⁺O⁻",
  "S-Oxidation": "→S=O",
  "Glucuronidation": "→GlcUA",
  "Sulfation": "→OSO₃",
  "Acetylation": "→NAc",
  "Methylation": "→Me",
};

function EnzymeTag({ enzyme }: { enzyme: string | null }) {
  if (!enzyme) return null;
  const enzymes = enzyme.split("/").slice(0, 2);
  return (
    <div className="flex flex-wrap gap-1 mt-1">
      {enzymes.map((e) => (
        <span key={e} className="text-[10px] font-mono bg-muted px-1 rounded">{e.trim()}</span>
      ))}
    </div>
  );
}

function MetaboliteNode({ m, index }: { m: Metabolite; index: number }) {
  const phase = m.phase ?? 1;
  const colors = PHASE_COLORS[phase] ?? PHASE_COLORS[1];
  const reactionIcon = m.reactionType ? (REACTION_ICONS[m.reactionType] ?? "→") : "→";
  const confidence = m.confidenceScore != null ? Math.round(m.confidenceScore * 100) : null;

  return (
    <div className={`rounded-lg border ${colors.border} ${colors.bg} p-3 min-w-[200px] max-w-[240px] text-sm`}>
      <div className="flex items-start justify-between gap-2">
        <div className="flex-1 min-w-0">
          <p className="font-semibold truncate text-xs">{m.metaboliteName ?? `M${index + 1}`}</p>
          <p className="font-mono text-[10px] text-muted-foreground truncate mt-0.5">{m.metaboliteSmiles}</p>
        </div>
        {confidence != null && (
          <span className={`text-[10px] font-bold shrink-0 ${colors.text}`}>{confidence}%</span>
        )}
      </div>
      <div className="mt-2 flex items-center gap-1 flex-wrap">
        <Badge variant="outline" className={`text-[10px] px-1 py-0 ${colors.text} ${colors.border}`}>
          Phase {phase}
        </Badge>
        {m.reactionType && (
          <span className={`text-[10px] ${colors.text}`}>{m.reactionType}</span>
        )}
      </div>
      <EnzymeTag enzyme={m.enzyme} />
      {m.molecularWeight != null && (
        <p className="text-[10px] text-muted-foreground mt-1">MW: {m.molecularWeight.toFixed(1)} Da{m.logP != null ? ` · LogP: ${m.logP.toFixed(2)}` : ""}</p>
      )}
    </div>
  );
}

export function MetabolitePathwayTree({ parentName, parentSmiles, metabolites }: MetabolitePathwayTreeProps) {
  const [filter, setFilter] = useState<"all" | 1 | 2>("all");

  const phase1 = useMemo(() => metabolites.filter((m) => (m.phase ?? 1) === 1), [metabolites]);
  const phase2 = useMemo(() => metabolites.filter((m) => (m.phase ?? 1) === 2), [metabolites]);
  const displayed = filter === "all" ? metabolites : metabolites.filter((m) => (m.phase ?? 1) === filter);

  const exportCsv = () => {
    const header = "metaboliteName,smiles,reactionType,enzyme,phase,molecularWeight,logP,confidence";
    const rows = metabolites.map((m) =>
      [
        m.metaboliteName ?? "",
        m.metaboliteSmiles,
        m.reactionType ?? "",
        m.enzyme ?? "",
        m.phase ?? 1,
        m.molecularWeight?.toFixed(2) ?? "",
        m.logP?.toFixed(2) ?? "",
        m.confidenceScore != null ? (m.confidenceScore * 100).toFixed(0) + "%" : "",
      ].map((v) => `"${v}"`).join(",")
    );
    const blob = new Blob([[header, ...rows].join("\n")], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `metabolites-${parentName.replace(/\s+/g, "-")}-${Date.now()}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  };

  if (metabolites.length === 0) {
    return (
      <Card>
        <CardContent className="flex flex-col items-center justify-center py-12 text-muted-foreground gap-2">
          <FlaskConical className="h-8 w-8 opacity-30" />
          <p className="text-sm">No metabolites predicted yet. Run metabolite prediction first.</p>
        </CardContent>
      </Card>
    );
  }

  return (
    <TooltipProvider>
      <Card>
        <CardHeader className="pb-3">
          <div className="flex items-start justify-between flex-wrap gap-2">
            <div>
              <CardTitle className="text-base">Metabolite Pathway</CardTitle>
              <CardDescription className="text-xs mt-1">
                {phase1.length} Phase I · {phase2.length} Phase II metabolites predicted via SMARTS engine
              </CardDescription>
            </div>
            <div className="flex items-center gap-2">
              <div className="flex rounded-md border overflow-hidden text-xs">
                {(["all", 1, 2] as const).map((f) => (
                  <button
                    key={f}
                    onClick={() => setFilter(f)}
                    className={`px-2 py-1 transition-colors ${filter === f ? "bg-primary text-primary-foreground" : "hover:bg-muted"}`}
                  >
                    {f === "all" ? "All" : `Phase ${f}`}
                  </button>
                ))}
              </div>
              <Button variant="outline" size="sm" onClick={exportCsv} className="gap-1 text-xs h-7">
                <Download className="h-3 w-3" />
                CSV
              </Button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          {/* Tree layout: parent → metabolites */}
          <div className="overflow-x-auto pb-2">
            <div className="flex items-start gap-0 min-w-max">
              {/* Parent node */}
              <div className="flex flex-col items-center">
                <div className="rounded-lg border-2 border-primary/40 bg-primary/5 p-3 min-w-[200px] max-w-[240px]">
                  <p className="font-bold text-sm truncate">{parentName}</p>
                  <p className="font-mono text-[10px] text-muted-foreground truncate mt-0.5">{parentSmiles}</p>
                  <Badge variant="outline" className="mt-2 text-[10px] px-1 py-0">Parent Compound</Badge>
                </div>
              </div>

              {/* Arrow + metabolites column */}
              <div className="flex items-center self-center mx-2">
                <div className="h-px w-8 bg-border" />
                <ChevronRight className="h-4 w-4 text-muted-foreground -ml-1" />
              </div>

              {/* Metabolite nodes stacked vertically */}
              <div className="flex flex-col gap-3">
                {displayed.map((m, i) => (
                  <div key={m.id} className="flex items-center gap-2">
                    <div className="flex flex-col items-center">
                      <div className="h-px w-6 bg-border" />
                    </div>
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <div className="cursor-default">
                          <MetaboliteNode m={m} index={i} />
                        </div>
                      </TooltipTrigger>
                      <TooltipContent side="right" className="max-w-xs text-xs space-y-1">
                        <p><strong>SMILES:</strong> <span className="font-mono">{m.metaboliteSmiles}</span></p>
                        {m.reactionType && <p><strong>Reaction:</strong> {m.reactionType}</p>}
                        {m.enzyme && <p><strong>Enzyme(s):</strong> {m.enzyme}</p>}
                        {m.predictedActivity && <p><strong>Activity:</strong> {m.predictedActivity}</p>}
                        {m.molecularWeight != null && <p><strong>MW:</strong> {m.molecularWeight.toFixed(2)} Da</p>}
                        {m.logP != null && <p><strong>LogP:</strong> {m.logP.toFixed(2)}</p>}
                      </TooltipContent>
                    </Tooltip>
                  </div>
                ))}
              </div>
            </div>
          </div>

          {/* Summary table */}
          <div className="mt-4 border rounded-md overflow-hidden">
            <table className="w-full text-xs">
              <thead className="bg-muted/50">
                <tr>
                  <th className="text-left px-3 py-2 font-medium">Metabolite</th>
                  <th className="text-left px-3 py-2 font-medium">Reaction</th>
                  <th className="text-left px-3 py-2 font-medium">Enzyme</th>
                  <th className="text-left px-3 py-2 font-medium">Phase</th>
                  <th className="text-right px-3 py-2 font-medium">MW</th>
                  <th className="text-right px-3 py-2 font-medium">LogP</th>
                  <th className="text-right px-3 py-2 font-medium">Confidence</th>
                </tr>
              </thead>
              <tbody className="divide-y">
                {displayed.map((m, i) => {
                  const phase = m.phase ?? 1;
                  const colors = PHASE_COLORS[phase] ?? PHASE_COLORS[1];
                  return (
                    <tr key={m.id} className="hover:bg-muted/30">
                      <td className="px-3 py-2 font-medium">{m.metaboliteName ?? `M${i + 1}`}</td>
                      <td className={`px-3 py-2 ${colors.text}`}>{m.reactionType ?? "—"}</td>
                      <td className="px-3 py-2 font-mono text-muted-foreground">{m.enzyme?.split("/")[0] ?? "—"}</td>
                      <td className="px-3 py-2">
                        <Badge variant="outline" className={`text-[10px] px-1 py-0 ${colors.text} ${colors.border}`}>
                          Phase {phase}
                        </Badge>
                      </td>
                      <td className="px-3 py-2 text-right font-mono">{m.molecularWeight?.toFixed(1) ?? "—"}</td>
                      <td className="px-3 py-2 text-right font-mono">{m.logP?.toFixed(2) ?? "—"}</td>
                      <td className="px-3 py-2 text-right">
                        {m.confidenceScore != null
                          ? <span className={colors.text + " font-semibold"}>{(m.confidenceScore * 100).toFixed(0)}%</span>
                          : "—"}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>
    </TooltipProvider>
  );
}
