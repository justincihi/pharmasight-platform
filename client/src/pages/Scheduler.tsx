import { useState } from "react";
import { Button } from "@/components/ui/button";
import { GlassCard } from "@/components/GlassCard";
import { trpc } from "@/lib/trpc";
import { toast } from "sonner";
import { Play, Clock, CheckCircle, Calendar } from "lucide-react";

export default function Scheduler() {
  const [isRunning, setIsRunning] = useState(false);

  const statusQuery = trpc.scheduler.status.useQuery(undefined, {
    refetchInterval: 5000, // Refresh every 5 seconds
  });

  const runNowMutation = trpc.scheduler.runNow.useMutation({
    onSuccess: () => {
      toast.success("Autonomous research engine started!");
      setIsRunning(false);
      statusQuery.refetch();
    },
    onError: (error) => {
      toast.error(`Failed to run scheduler: ${error.message}`);
      setIsRunning(false);
    },
  });

  const handleRunNow = () => {
    setIsRunning(true);
    runNowMutation.mutate();
  };

  const status = statusQuery.data;

  return (
    <div className="p-6 space-y-6">
      <div>
        <h1 className="text-3xl font-bold bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent">
          Autonomous Research Scheduler
        </h1>
        <p className="text-muted-foreground mt-2">
          Configure and manage automated compound discovery
        </p>
      </div>

      {/* Status Card */}
      <GlassCard className="p-6">
        <div className="flex items-center justify-between mb-6">
          <div className="flex items-center gap-3">
            <div className={`w-3 h-3 rounded-full ${status?.running ? 'bg-green-500 animate-pulse' : 'bg-gray-400'}`} />
            <h2 className="text-xl font-semibold">
              Scheduler Status: {status?.running ? "Running" : "Stopped"}
            </h2>
          </div>
          <Button
            onClick={handleRunNow}
            disabled={isRunning || runNowMutation.isPending}
            size="lg"
            className="bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700"
          >
            {isRunning || runNowMutation.isPending ? (
              <>
                <Clock className="w-4 h-4 mr-2 animate-spin" />
                Running...
              </>
            ) : (
              <>
                <Play className="w-4 h-4 mr-2" />
                Run Now
              </>
            )}
          </Button>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <div className="flex items-center gap-3 p-4 bg-blue-50 dark:bg-blue-950/20 rounded-lg">
            <Calendar className="w-8 h-8 text-blue-600" />
            <div>
              <p className="text-sm text-muted-foreground">Schedule</p>
              <p className="font-semibold">{status?.config?.cronSchedule || "0 9 * * *"}</p>
              <p className="text-xs text-muted-foreground">Daily at 9 AM</p>
            </div>
          </div>

          <div className="flex items-center gap-3 p-4 bg-green-50 dark:bg-green-950/20 rounded-lg">
            <Clock className="w-8 h-8 text-green-600" />
            <div>
              <p className="text-sm text-muted-foreground">Next Run</p>
              <p className="font-semibold text-sm">
                {status?.nextRun ? new Date(status.nextRun).toLocaleString() : "Not scheduled"}
              </p>
            </div>
          </div>

          <div className="flex items-center gap-3 p-4 bg-purple-50 dark:bg-purple-950/20 rounded-lg">
            <CheckCircle className="w-8 h-8 text-purple-600" />
            <div>
              <p className="text-sm text-muted-foreground">Min Confidence</p>
              <p className="font-semibold">{status?.config?.minConfidenceForNotification || 85}%</p>
              <p className="text-xs text-muted-foreground">For notifications</p>
            </div>
          </div>
        </div>
      </GlassCard>

      {/* How It Works */}
      <GlassCard className="p-6">
        <h2 className="text-xl font-semibold mb-4">How It Works</h2>
        <div className="space-y-3 text-muted-foreground">
          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-full bg-blue-100 dark:bg-blue-900/30 flex items-center justify-center flex-shrink-0 mt-1">
              <span className="text-blue-600 font-semibold">1</span>
            </div>
            <div>
              <p className="font-semibold text-foreground">Research Goals Analysis</p>
              <p className="text-sm">The engine analyzes your research goals (psychedelics, nootropics, anxiolytics) and AI-discovered medical trends.</p>
            </div>
          </div>

          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-full bg-purple-100 dark:bg-purple-900/30 flex items-center justify-center flex-shrink-0 mt-1">
              <span className="text-purple-600 font-semibold">2</span>
            </div>
            <div>
              <p className="font-semibold text-foreground">Compound Discovery</p>
              <p className="text-sm">Uses Perplexity AI and Gemini to search literature, generate novel analogs with RDKit, and predict ADMET properties.</p>
            </div>
          </div>

          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-full bg-green-100 dark:bg-green-900/30 flex items-center justify-center flex-shrink-0 mt-1">
              <span className="text-green-600 font-semibold">3</span>
            </div>
            <div>
              <p className="font-semibold text-foreground">Import & Notify</p>
              <p className="text-sm">Imports discoveries into the database and sends notifications for high-confidence compounds (≥85%).</p>
            </div>
          </div>
        </div>
      </GlassCard>
    </div>
  );
}
