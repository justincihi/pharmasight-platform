import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { trpc } from "@/lib/trpc";
import { Activity, Calendar, CheckCircle2, Clock, PlayCircle, RefreshCw, AlertCircle } from "lucide-react";
import { toast } from "sonner";

export default function SchedulerDashboard() {
  const { data: schedulerStatus, isLoading, refetch } = trpc.scheduler.status.useQuery();
  const runNow = trpc.scheduler.runNow.useMutation({
    onSuccess: () => {
      toast.success("Scheduler triggered successfully");
      refetch();
    },
    onError: (error) => {
      toast.error(`Failed to run scheduler: ${error.message}`);
    },
  });

  if (isLoading) {
    return (
      <div className="container mx-auto py-8">
        <div className="flex items-center justify-center h-64">
          <RefreshCw className="h-8 w-8 animate-spin text-muted-foreground" />
        </div>
      </div>
    );
  }

  const getStatusColor = (running: boolean) => {
    return running
      ? "bg-green-100 text-green-800 border-green-200"
      : "bg-blue-100 text-blue-800 border-blue-200";
  };

  const formatDate = (dateString: string | null) => {
    if (!dateString) return "Never";
    return new Date(dateString).toLocaleString();
  };

  return (
    <div className="container mx-auto py-8 space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold">Scheduler Dashboard</h1>
          <p className="text-muted-foreground mt-1">
            Monitor autonomous research engine and scheduled tasks
          </p>
        </div>
        <Button
          onClick={() => runNow.mutate()}
          disabled={runNow.isPending || schedulerStatus?.running}
        >
          <PlayCircle className="h-4 w-4 mr-2" />
          Run Now
        </Button>
      </div>

      {/* Status Overview */}
      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Status</CardTitle>
            <Activity className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <Badge className={getStatusColor(schedulerStatus?.running || false)}>
              {schedulerStatus?.running ? "RUNNING" : "IDLE"}
            </Badge>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Schedule</CardTitle>
            <Clock className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-sm font-medium">
              {schedulerStatus?.config?.cronSchedule || "Not configured"}
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Next Run</CardTitle>
            <Calendar className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-sm font-medium">
              {formatDate(schedulerStatus?.nextRun || null)}
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Enabled</CardTitle>
            <CheckCircle2 className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <Badge variant={schedulerStatus?.config?.enabled ? "default" : "outline"}>
              {schedulerStatus?.config?.enabled ? "Yes" : "No"}
            </Badge>
          </CardContent>
        </Card>
      </div>

      {/* Scheduler Configuration */}
      <Card>
        <CardHeader>
          <CardTitle>Scheduler Configuration</CardTitle>
          <CardDescription>
            Autonomous research engine configuration and status
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            <div className="flex items-center justify-between p-4 border rounded-lg">
              <div className="flex items-center gap-4">
                <div className="flex h-10 w-10 items-center justify-center rounded-full bg-primary/10">
                  <Activity className="h-5 w-5 text-primary" />
                </div>
                <div>
                  <h4 className="font-semibold">Autonomous Discovery Engine</h4>
                  <p className="text-sm text-muted-foreground">
                    Imports analog discoveries from research engine
                  </p>
                </div>
              </div>
              <div className="text-right">
                <Badge className={getStatusColor(schedulerStatus?.running || false)}>
                  {schedulerStatus?.running ? "RUNNING" : "IDLE"}
                </Badge>
                <p className="text-xs text-muted-foreground mt-1">
                  Next: {formatDate(schedulerStatus?.nextRun || null)}
                </p>
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Configuration Details */}
      <Card>
        <CardHeader>
          <CardTitle>Configuration Details</CardTitle>
          <CardDescription>Scheduler settings and parameters</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid gap-4 md:grid-cols-2">
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Cron Schedule</p>
              <p className="text-sm font-mono">{schedulerStatus?.config?.cronSchedule || "N/A"}</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Timezone</p>
              <p className="text-sm">America/New_York (EST)</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Enabled</p>
              <p className="text-sm">{schedulerStatus?.config?.enabled ? "Yes" : "No"}</p>
            </div>
            <div className="space-y-1">
              <p className="text-sm font-medium text-muted-foreground">Next Execution</p>
              <p className="text-sm">{formatDate(schedulerStatus?.nextRun || null)}</p>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
