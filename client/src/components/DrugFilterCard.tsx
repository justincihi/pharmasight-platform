import { useState } from "react";
import { GlassCard } from "./GlassCard";
import { Button } from "./ui/button";
import { Badge } from "./ui/badge";
import { Loader2, ChevronDown, ChevronUp, AlertTriangle, CheckCircle, Brain, Activity } from "lucide-react";
import { trpc } from "@/lib/trpc";
import { motion, AnimatePresence } from "framer-motion";

interface DrugFilterCardProps {
  smiles: string;
  compoundName: string;
}

interface PainsResult {
  pains_violations: string[];
  brenk_violations: string[];
  passes_filters: boolean;
  message: string;
}

interface CnsResult {
  cns_mpo_score: number;
  mpo_score: number;
  interpretation: string;
  components: {
    logp: number;
    mw: number;
    tpsa: number;
    hbd: number;
    hba: number;
  };
}

interface BbbResult {
  bbb_permeable: boolean;
  probability: number;
  interpretation: string;
}

export function DrugFilterCard({ smiles, compoundName }: DrugFilterCardProps) {
  const [isExpanded, setIsExpanded] = useState(false);
  const [activeFilter, setActiveFilter] = useState<"pains" | "cns" | "bbb" | null>(null);

  const painsMutation = trpc.drugFilters.painsBrenk.useMutation();
  const cnsMutation = trpc.drugFilters.cnsMpo.useMutation();
  const bbbMutation = trpc.drugFilters.bbbPermeability.useMutation();

  const handleRunFilter = async (filterType: "pains" | "cns" | "bbb") => {
    setActiveFilter(filterType);
    setIsExpanded(true);

    if (filterType === "pains") {
      await painsMutation.mutateAsync({ smiles });
    } else if (filterType === "cns") {
      await cnsMutation.mutateAsync({ smiles });
    } else if (filterType === "bbb") {
      await bbbMutation.mutateAsync({ smiles });
    }
  };

  const isLoading = painsMutation.isPending || cnsMutation.isPending || bbbMutation.isPending;

  const isPainsResult = (data: unknown): data is PainsResult => {
    return typeof data === 'object' && data !== null && 'passes_filters' in data;
  };

  const isCnsResult = (data: unknown): data is CnsResult => {
    return typeof data === 'object' && data !== null && 'cns_mpo_score' in data;
  };

  const isBbbResult = (data: unknown): data is BbbResult => {
    return typeof data === 'object' && data !== null && 'bbb_permeable' in data;
  };

  return (
    <GlassCard className="p-6">
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-purple-500 to-pink-500 flex items-center justify-center">
            <Activity className="w-5 h-5 text-white" />
          </div>
          <div>
            <h3 className="text-lg font-semibold text-foreground">Drug Filters & Scores</h3>
            <p className="text-sm text-muted-foreground">Structural alerts, CNS MPO, BBB permeability</p>
          </div>
        </div>
        <Button
          variant="ghost"
          size="sm"
          onClick={() => setIsExpanded(!isExpanded)}
        >
          {isExpanded ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
        </Button>
      </div>

      {/* Filter Buttons */}
      <div className="grid grid-cols-3 gap-3 mb-4">
        <Button
          variant="outline"
          size="sm"
          onClick={() => handleRunFilter("pains")}
          disabled={isLoading}
          className="glass border-glass-border"
        >
          {painsMutation.isPending ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <AlertTriangle className="w-4 h-4 mr-2" />}
          PAINS/Brenk
        </Button>
        <Button
          variant="outline"
          size="sm"
          onClick={() => handleRunFilter("cns")}
          disabled={isLoading}
          className="glass border-glass-border"
        >
          {cnsMutation.isPending ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Brain className="w-4 h-4 mr-2" />}
          CNS MPO
        </Button>
        <Button
          variant="outline"
          size="sm"
          onClick={() => handleRunFilter("bbb")}
          disabled={isLoading}
          className="glass border-glass-border"
        >
          {bbbMutation.isPending ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Activity className="w-4 h-4 mr-2" />}
          BBB
        </Button>
      </div>

      {/* Results */}
      <AnimatePresence>
        {isExpanded && activeFilter && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: "auto" }}
            exit={{ opacity: 0, height: 0 }}
            className="space-y-4 pt-4 border-t border-glass-border"
          >
            {activeFilter === "pains" && painsMutation.data && isPainsResult(painsMutation.data) ? (
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <span className="font-medium text-foreground">PAINS/Brenk Filter</span>
                  <Badge variant={painsMutation.data.passes_filters ? "default" : "destructive"}>
                    {painsMutation.data.passes_filters ? "PASS" : "FAIL"}
                  </Badge>
                </div>
                {(painsMutation.data.pains_violations.length > 0 || painsMutation.data.brenk_violations.length > 0) && (
                  <div className="space-y-2">
                    <p className="text-sm text-muted-foreground">Structural Alerts:</p>
                    {[...painsMutation.data.pains_violations, ...painsMutation.data.brenk_violations].map((alert: string, idx: number) => (
                      <div key={idx} className="flex items-start gap-2 text-sm">
                        <AlertTriangle className="w-4 h-4 text-orange-500 mt-0.5" />
                        <span className="text-foreground">{alert}</span>
                      </div>
                    ))}
                  </div>
                )}
                {painsMutation.data.passes_filters && (
                  <div className="flex items-center gap-2 text-sm text-green-600">
                    <CheckCircle className="w-4 h-4" />
                    <span>No structural alerts detected</span>
                  </div>
                )}
              </div>
            ) : null}

            {activeFilter === "cns" && cnsMutation.data && isCnsResult(cnsMutation.data) ? (
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <span className="font-medium text-foreground">CNS MPO Score</span>
                  <Badge variant={cnsMutation.data.cns_mpo_score >= 4 ? "default" : "secondary"}>
                    {cnsMutation.data.cns_mpo_score.toFixed(2)} / 6.0
                  </Badge>
                </div>
                <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-2">
                  <div
                    className="bg-gradient-to-r from-blue-500 to-indigo-500 h-2 rounded-full transition-all"
                    style={{ width: `${(cnsMutation.data.cns_mpo_score / 6) * 100}%` }}
                  />
                </div>
                <p className="text-sm text-muted-foreground">
                  {cnsMutation.data.cns_mpo_score >= 5
                    ? "Excellent CNS drug-likeness"
                    : cnsMutation.data.cns_mpo_score >= 4
                    ? "Good CNS drug-likeness"
                    : "Poor CNS drug-likeness"}
                </p>
                {cnsMutation.data.components && (
                  <div className="grid grid-cols-2 gap-2 text-sm">
                    <div>
                      <span className="text-muted-foreground">LogP:</span>{" "}
                      <span className="text-foreground font-medium">{cnsMutation.data.components.logp?.toFixed(2)}</span>
                    </div>
                    <div>
                      <span className="text-muted-foreground">TPSA:</span>{" "}
                      <span className="text-foreground font-medium">{cnsMutation.data.components.tpsa?.toFixed(1)}</span>
                    </div>
                  </div>
                )}
              </div>
            ) : null}

            {activeFilter === "bbb" && bbbMutation.data && isBbbResult(bbbMutation.data) ? (
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <span className="font-medium text-foreground">BBB Permeability</span>
                  <Badge variant={bbbMutation.data.bbb_permeable ? "default" : "secondary"}>
                    {bbbMutation.data.bbb_permeable ? "HIGH" : "LOW"}
                  </Badge>
                </div>
                <div className="flex items-center gap-2">
                  <span className="text-sm text-muted-foreground">Probability:</span>
                  <span className="text-foreground font-medium">{(bbbMutation.data.probability * 100).toFixed(1)}%</span>
                </div>
                <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-2">
                  <div
                    className="bg-gradient-to-r from-green-500 to-emerald-500 h-2 rounded-full transition-all"
                    style={{ width: `${bbbMutation.data.probability * 100}%` }}
                  />
                </div>
                <p className="text-sm text-muted-foreground">
                  {bbbMutation.data.bbb_permeable
                    ? "Likely to cross the blood-brain barrier"
                    : "Unlikely to cross the blood-brain barrier"}
                </p>
              </div>
            ) : null}
          </motion.div>
        )}
      </AnimatePresence>
    </GlassCard>
  );
}
