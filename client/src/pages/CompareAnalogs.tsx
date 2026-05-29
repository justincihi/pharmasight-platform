import DashboardLayout from "@/components/DashboardLayout";
import CompoundComparisonChart from "@/components/CompoundComparisonChart";
import { GlassCard } from "@/components/GlassCard";
import { BarChart3, Info } from "lucide-react";

export default function CompareAnalogs() {
  return (
    <DashboardLayout>
      <div className="p-6 max-w-6xl mx-auto space-y-6">
        {/* Page Header */}
        <div>
          <h1 className="text-3xl font-bold bg-gradient-to-r from-purple-600 to-blue-600 bg-clip-text text-transparent">
            Compound Comparison
          </h1>
          <p className="text-muted-foreground mt-2">
            Side-by-side radar chart comparison of up to 4 analogs across key pharmaceutical
            metrics
          </p>
        </div>

        {/* Info Banner */}
        <div className="flex items-start gap-3 p-4 rounded-lg bg-blue-50 dark:bg-blue-950/20 border border-blue-200 dark:border-blue-800">
          <Info className="w-5 h-5 text-blue-600 flex-shrink-0 mt-0.5" />
          <div className="text-sm text-blue-800 dark:text-blue-200">
            <span className="font-semibold">How to use:</span> Browse the compound list below and
            click any card to add it to the comparison. Select 2–4 analogs to generate the radar
            chart. Metrics include Confidence, Safety, Efficacy, Drug-likeness, and Structural
            Similarity scores.
          </div>
        </div>

        {/* Main Chart Component */}
        <GlassCard className="p-6">
          <CompoundComparisonChart />
        </GlassCard>
      </div>
    </DashboardLayout>
  );
}
