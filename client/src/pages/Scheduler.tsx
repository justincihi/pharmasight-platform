import { useState, useEffect, useCallback } from "react";
import { Button } from "@/components/ui/button";
import { GlassCard } from "@/components/GlassCard";
import { trpc } from "@/lib/trpc";
import { toast } from "sonner";
import {
  Play,
  Clock,
  CheckCircle,
  Calendar,
  Settings,
  RefreshCw,
  Timer,
  Zap,
  ChevronDown,
} from "lucide-react";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import { Badge } from "@/components/ui/badge";

// Preset schedule options
const SCHEDULE_PRESETS = [
  { label: "Every 6 Hours", cron: "0 */6 * * *", description: "Runs 4× per day" },
  { label: "Every 12 Hours", cron: "0 */12 * * *", description: "Runs twice daily" },
  { label: "Daily at 9 AM", cron: "0 9 * * *", description: "Once per day (recommended)" },
  { label: "Daily at Midnight", cron: "0 0 * * *", description: "Once per day at midnight" },
  { label: "Every 2 Days", cron: "0 9 */2 * *", description: "Every other day at 9 AM" },
  { label: "Weekly (Monday)", cron: "0 9 * * 1", description: "Once per week on Monday" },
];

function getPresetLabel(cron: string): string {
  const preset = SCHEDULE_PRESETS.find((p) => p.cron === cron);
  return preset?.label ?? cron;
}

/** Parse a cron string and compute the next fire time */
function getNextFireTime(cron: string): Date | null {
  try {
    const parts = cron.trim().split(/\s+/);
    if (parts.length !== 5) return null;
    const [minute, hour, dom, , dow] = parts;

    const now = new Date();
    const next = new Date(now);
    next.setSeconds(0);
    next.setMilliseconds(0);

    // Handle simple interval patterns like */6
    if (hour.startsWith("*/")) {
      const interval = parseInt(hour.slice(2), 10);
      const currentHour = now.getHours();
      const minuteVal = parseInt(minute, 10) || 0;
      let nextHour = Math.ceil((currentHour + (now.getMinutes() > minuteVal ? 1 : 0)) / interval) * interval;
      if (nextHour <= currentHour && now.getMinutes() >= minuteVal) nextHour += interval;
      if (nextHour >= 24) {
        next.setDate(next.getDate() + 1);
        nextHour = nextHour % 24;
      }
      next.setHours(nextHour, minuteVal, 0, 0);
      return next;
    }

    // Simple daily at fixed hour
    const hourVal = parseInt(hour, 10);
    const minuteVal = parseInt(minute, 10);
    if (!isNaN(hourVal) && !isNaN(minuteVal)) {
      next.setHours(hourVal, minuteVal, 0, 0);
      if (next <= now) {
        // Handle weekly
        if (dow !== "*") {
          const targetDow = parseInt(dow, 10);
          let daysAhead = (targetDow - now.getDay() + 7) % 7 || 7;
          next.setDate(next.getDate() + daysAhead);
        } else if (dom.startsWith("*/")) {
          const interval = parseInt(dom.slice(2), 10);
          next.setDate(next.getDate() + interval);
        } else {
          next.setDate(next.getDate() + 1);
        }
      }
      return next;
    }
    return null;
  } catch {
    return null;
  }
}

/** Format a duration in seconds into a human-readable countdown */
function formatCountdown(seconds: number): string {
  if (seconds <= 0) return "Imminent";
  const d = Math.floor(seconds / 86400);
  const h = Math.floor((seconds % 86400) / 3600);
  const m = Math.floor((seconds % 3600) / 60);
  const s = seconds % 60;
  if (d > 0) return `${d}d ${h}h ${m}m`;
  if (h > 0) return `${h}h ${m}m ${s}s`;
  if (m > 0) return `${m}m ${s}s`;
  return `${s}s`;
}

export default function Scheduler() {
  const [isRunning, setIsRunning] = useState(false);
  const [countdown, setCountdown] = useState<string>("—");
  const [nextFireDate, setNextFireDate] = useState<Date | null>(null);

  const statusQuery = trpc.scheduler.status.useQuery(undefined, {
    refetchInterval: 10000,
  });

  const cronQuery = trpc.scheduler.getCronInterval.useQuery(undefined, {
    refetchInterval: 30000,
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

  const setCronMutation = trpc.scheduler.setCronInterval.useMutation({
    onSuccess: (data) => {
      toast.success(`Schedule updated to: ${getPresetLabel(data.cronSchedule)}`);
      cronQuery.refetch();
      statusQuery.refetch();
    },
    onError: (error) => {
      toast.error(`Failed to update schedule: ${error.message}`);
    },
  });

  const startMutation = trpc.scheduler.start.useMutation({
    onSuccess: () => {
      toast.success("Scheduler started");
      statusQuery.refetch();
    },
    onError: (err) => toast.error(`Failed to start: ${err.message}`),
  });

  const stopMutation = trpc.scheduler.stop.useMutation({
    onSuccess: () => {
      toast.success("Scheduler stopped");
      statusQuery.refetch();
    },
    onError: (err) => toast.error(`Failed to stop: ${err.message}`),
  });

  const handleRunNow = () => {
    setIsRunning(true);
    runNowMutation.mutate();
  };

  const handleSetSchedule = (cron: string) => {
    setCronMutation.mutate({ cronSchedule: cron });
  };

  // Compute next fire time from cron string
  const updateNextFire = useCallback(() => {
    const cron = cronQuery.data?.cronSchedule ?? statusQuery.data?.config?.cronSchedule ?? "0 9 * * *";
    // Prefer server-provided nextRun if available
    if (statusQuery.data?.nextRun) {
      const serverNext = new Date(statusQuery.data.nextRun);
      if (serverNext > new Date()) {
        setNextFireDate(serverNext);
        return;
      }
    }
    const computed = getNextFireTime(cron);
    setNextFireDate(computed);
  }, [cronQuery.data, statusQuery.data]);

  useEffect(() => {
    updateNextFire();
  }, [updateNextFire]);

  // Live countdown ticker
  useEffect(() => {
    const tick = () => {
      if (!nextFireDate) {
        setCountdown("—");
        return;
      }
      const diff = Math.max(0, Math.floor((nextFireDate.getTime() - Date.now()) / 1000));
      setCountdown(formatCountdown(diff));
    };
    tick();
    const interval = setInterval(tick, 1000);
    return () => clearInterval(interval);
  }, [nextFireDate]);

  const status = statusQuery.data;
  const currentCron = cronQuery.data?.cronSchedule ?? status?.config?.cronSchedule ?? "0 9 * * *";
  const currentPreset = SCHEDULE_PRESETS.find((p) => p.cron === currentCron);

  return (
    <div className="p-6 space-y-6 max-w-5xl">
      <div>
        <h1 className="text-3xl font-bold bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent">
          Autonomous Research Scheduler
        </h1>
        <p className="text-muted-foreground mt-2">
          Configure and manage automated compound discovery
        </p>
      </div>

      {/* Status + Run Now */}
      <GlassCard className="p-6">
        <div className="flex items-center justify-between mb-6 flex-wrap gap-3">
          <div className="flex items-center gap-3">
            <div
              className={`w-3 h-3 rounded-full ${
                status?.running ? "bg-green-500 animate-pulse" : "bg-gray-400"
              }`}
            />
            <h2 className="text-xl font-semibold">
              Scheduler:{" "}
              <span className={status?.running ? "text-green-600" : "text-muted-foreground"}>
                {status?.running ? "Active" : "Stopped"}
              </span>
            </h2>
            {status?.running && (
              <Badge variant="outline" className="text-green-600 border-green-500">
                Running
              </Badge>
            )}
          </div>
          <div className="flex items-center gap-2">
            {status?.running ? (
              <Button
                variant="outline"
                size="sm"
                onClick={() => stopMutation.mutate()}
                disabled={stopMutation.isPending}
              >
                Stop Scheduler
              </Button>
            ) : (
              <Button
                variant="outline"
                size="sm"
                onClick={() => startMutation.mutate()}
                disabled={startMutation.isPending}
              >
                <Zap className="w-4 h-4 mr-1" />
                Start Scheduler
              </Button>
            )}
            <Button
              onClick={handleRunNow}
              disabled={isRunning || runNowMutation.isPending}
              size="lg"
              className="bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700"
            >
              {isRunning || runNowMutation.isPending ? (
                <>
                  <RefreshCw className="w-4 h-4 mr-2 animate-spin" />
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
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          {/* Countdown Timer */}
          <div className="flex items-center gap-3 p-4 bg-blue-50 dark:bg-blue-950/20 rounded-lg">
            <Timer className="w-8 h-8 text-blue-600 flex-shrink-0" />
            <div className="min-w-0">
              <p className="text-sm text-muted-foreground">Next Run In</p>
              <p className="font-bold text-lg text-blue-600 font-mono tabular-nums">
                {countdown}
              </p>
              {nextFireDate && (
                <p className="text-xs text-muted-foreground truncate">
                  {nextFireDate.toLocaleString()}
                </p>
              )}
            </div>
          </div>

          {/* Current Schedule */}
          <div className="flex items-center gap-3 p-4 bg-purple-50 dark:bg-purple-950/20 rounded-lg">
            <Calendar className="w-8 h-8 text-purple-600 flex-shrink-0" />
            <div className="min-w-0 flex-1">
              <p className="text-sm text-muted-foreground">Schedule</p>
              <p className="font-semibold truncate">{currentPreset?.label ?? currentCron}</p>
              <p className="text-xs text-muted-foreground">
                {currentPreset?.description ?? "Custom cron"}
              </p>
            </div>
          </div>

          {/* Min Confidence */}
          <div className="flex items-center gap-3 p-4 bg-green-50 dark:bg-green-950/20 rounded-lg">
            <CheckCircle className="w-8 h-8 text-green-600 flex-shrink-0" />
            <div>
              <p className="text-sm text-muted-foreground">Min Confidence</p>
              <p className="font-semibold">{status?.config?.minConfidenceForNotification ?? 85}%</p>
              <p className="text-xs text-muted-foreground">For notifications</p>
            </div>
          </div>
        </div>
      </GlassCard>

      {/* Interval Configurator */}
      <GlassCard className="p-6">
        <div className="flex items-center gap-2 mb-5">
          <Settings className="w-5 h-5 text-muted-foreground" />
          <h2 className="text-xl font-semibold">Schedule Configuration</h2>
        </div>

        <p className="text-sm text-muted-foreground mb-4">
          Choose how often the autonomous research engine runs. More frequent runs discover
          compounds faster but consume more API credits. Daily is recommended for cost-efficiency.
        </p>

        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
          {SCHEDULE_PRESETS.map((preset) => {
            const isActive = currentCron === preset.cron;
            return (
              <button
                key={preset.cron}
                onClick={() => handleSetSchedule(preset.cron)}
                disabled={setCronMutation.isPending}
                className={`text-left p-4 rounded-lg border-2 transition-all ${
                  isActive
                    ? "border-blue-500 bg-blue-50 dark:bg-blue-950/30"
                    : "border-border hover:border-blue-300 hover:bg-muted/50"
                }`}
              >
                <div className="flex items-center justify-between mb-1">
                  <span className="font-semibold text-sm">{preset.label}</span>
                  {isActive && (
                    <Badge className="bg-blue-500 text-white text-xs">Active</Badge>
                  )}
                </div>
                <p className="text-xs text-muted-foreground">{preset.description}</p>
                <p className="text-xs font-mono text-muted-foreground mt-1">{preset.cron}</p>
              </button>
            );
          })}
        </div>

        {/* Custom cron input */}
        <div className="mt-4 pt-4 border-t border-border">
          <p className="text-sm font-medium mb-2">Custom Cron Expression</p>
          <div className="flex gap-2">
            <input
              type="text"
              defaultValue={currentCron}
              id="custom-cron-input"
              placeholder="e.g. 0 9 * * *"
              className="flex-1 px-3 py-2 text-sm border border-border rounded-md bg-background font-mono"
            />
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                const input = document.getElementById("custom-cron-input") as HTMLInputElement;
                if (input?.value) handleSetSchedule(input.value.trim());
              }}
              disabled={setCronMutation.isPending}
            >
              Apply
            </Button>
          </div>
          <p className="text-xs text-muted-foreground mt-1">
            Format: minute hour day-of-month month day-of-week (UTC)
          </p>
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
              <p className="text-sm">
                The engine analyzes your research goals (psychedelics, nootropics, anxiolytics)
                and AI-discovered medical trends.
              </p>
            </div>
          </div>

          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-full bg-purple-100 dark:bg-purple-900/30 flex items-center justify-center flex-shrink-0 mt-1">
              <span className="text-purple-600 font-semibold">2</span>
            </div>
            <div>
              <p className="font-semibold text-foreground">Compound Discovery</p>
              <p className="text-sm">
                Uses Perplexity AI and Gemini to search literature, generate novel analogs, and
                predict ADMET properties.
              </p>
            </div>
          </div>

          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-full bg-green-100 dark:bg-green-900/30 flex items-center justify-center flex-shrink-0 mt-1">
              <span className="text-green-600 font-semibold">3</span>
            </div>
            <div>
              <p className="font-semibold text-foreground">Import &amp; Notify</p>
              <p className="text-sm">
                Imports discoveries into the database and sends notifications for high-confidence
                compounds (≥{status?.config?.minConfidenceForNotification ?? 85}%).
              </p>
            </div>
          </div>
        </div>
      </GlassCard>
    </div>
  );
}
