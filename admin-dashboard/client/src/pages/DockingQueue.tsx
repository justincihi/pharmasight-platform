import { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { trpc } from '@/lib/trpc';
import { Loader2, Play, Pause, Trash2, RefreshCw, TrendingUp, TrendingDown } from 'lucide-react';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import { Link } from 'wouter';

export default function DockingQueue() {
  const [selectedTarget, setSelectedTarget] = useState<string>('all');
  const [selectedStatus, setSelectedStatus] = useState<string>('all');

  const { data: queueStatus, isLoading: statusLoading, refetch: refetchStatus } = 
    trpc.docking.getQueueStatus.useQuery(undefined, {
      refetchInterval: 5000, // Auto-refresh every 5 seconds
    });

  const { data: queueJobs, isLoading: jobsLoading, refetch: refetchJobs } = 
    trpc.docking.getQueueJobs.useQuery(
      { target: selectedTarget === 'all' ? undefined : selectedTarget, status: selectedStatus === 'all' ? undefined : selectedStatus },
      { refetchInterval: 5000 }
    );

  const clearCompletedMutation = trpc.docking.clearCompleted.useMutation({
    onSuccess: () => {
      refetchStatus();
      refetchJobs();
    },
  });

  const retryFailedMutation = trpc.docking.retryFailed.useMutation({
    onSuccess: () => {
      refetchStatus();
      refetchJobs();
    },
  });

  const handleClearCompleted = () => {
    if (confirm('Clear all completed jobs from the queue?')) {
      clearCompletedMutation.mutate();
    }
  };

  const handleRetryFailed = () => {
    if (confirm('Retry all failed docking jobs?')) {
      retryFailedMutation.mutate();
    }
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'completed':
        return 'default';
      case 'running':
        return 'secondary';
      case 'pending':
        return 'outline';
      case 'failed':
        return 'destructive';
      default:
        return 'outline';
    }
  };

  const getAffinityTrend = (affinity: string | null) => {
    if (!affinity) return null;
    const value = parseFloat(affinity);
    if (value < -8) return <TrendingDown className="w-4 h-4 text-green-600" />;
    if (value < -6) return <TrendingDown className="w-4 h-4 text-blue-600" />;
    return <TrendingUp className="w-4 h-4 text-gray-400" />;
  };

  if (statusLoading || jobsLoading) {
    return (
      <div className="flex items-center justify-center h-screen">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  return (
    <div className="container mx-auto py-8 space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold">Docking Queue</h1>
          <p className="text-muted-foreground mt-1">
            Monitor and manage molecular docking jobs across multiple receptor targets
          </p>
        </div>
        <div className="flex gap-2">
          <Button variant="outline" size="sm" onClick={() => { refetchStatus(); refetchJobs(); }}>
            <RefreshCw className="w-4 h-4 mr-2" />
            Refresh
          </Button>
          <Button 
            variant="outline" 
            size="sm" 
            onClick={handleRetryFailed}
            disabled={!queueStatus || queueStatus.failed === 0}
          >
            <Play className="w-4 h-4 mr-2" />
            Retry Failed
          </Button>
          <Button 
            variant="outline" 
            size="sm" 
            onClick={handleClearCompleted}
            disabled={!queueStatus || queueStatus.completed === 0}
          >
            <Trash2 className="w-4 h-4 mr-2" />
            Clear Completed
          </Button>
        </div>
      </div>

      {/* Queue Status Cards */}
      <div className="grid grid-cols-1 md:grid-cols-5 gap-4">
        <Card>
          <CardHeader className="pb-3">
            <CardDescription>Total Jobs</CardDescription>
            <CardTitle className="text-3xl">{queueStatus?.total || 0}</CardTitle>
          </CardHeader>
        </Card>
        <Card>
          <CardHeader className="pb-3">
            <CardDescription>Pending</CardDescription>
            <CardTitle className="text-3xl text-gray-600">{queueStatus?.pending || 0}</CardTitle>
          </CardHeader>
        </Card>
        <Card>
          <CardHeader className="pb-3">
            <CardDescription>Running</CardDescription>
            <CardTitle className="text-3xl text-blue-600">{queueStatus?.running || 0}</CardTitle>
          </CardHeader>
        </Card>
        <Card>
          <CardHeader className="pb-3">
            <CardDescription>Completed</CardDescription>
            <CardTitle className="text-3xl text-green-600">{queueStatus?.completed || 0}</CardTitle>
          </CardHeader>
        </Card>
        <Card>
          <CardHeader className="pb-3">
            <CardDescription>Failed</CardDescription>
            <CardTitle className="text-3xl text-red-600">{queueStatus?.failed || 0}</CardTitle>
          </CardHeader>
        </Card>
      </div>

      {/* Filters */}
      <Card>
        <CardHeader>
          <CardTitle>Filter Jobs</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex gap-4">
            <div className="flex gap-2">
              <span className="text-sm font-medium">Target:</span>
              {['all', 'NMDA', '5-HT2A', 'D2'].map((target) => (
                <Button
                  key={target}
                  variant={selectedTarget === target ? 'default' : 'outline'}
                  size="sm"
                  onClick={() => setSelectedTarget(target)}
                >
                  {target === 'all' ? 'All' : target}
                </Button>
              ))}
            </div>
            <div className="flex gap-2">
              <span className="text-sm font-medium">Status:</span>
              {['all', 'pending', 'running', 'completed', 'failed'].map((status) => (
                <Button
                  key={status}
                  variant={selectedStatus === status ? 'default' : 'outline'}
                  size="sm"
                  onClick={() => setSelectedStatus(status)}
                >
                  {status.charAt(0).toUpperCase() + status.slice(1)}
                </Button>
              ))}
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Jobs Table */}
      <Card>
        <CardHeader>
          <CardTitle>Docking Jobs</CardTitle>
          <CardDescription>
            {queueJobs?.length || 0} jobs matching filters
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Analog</TableHead>
                <TableHead>Target</TableHead>
                <TableHead>Status</TableHead>
                <TableHead>Priority</TableHead>
                <TableHead>Binding Affinity</TableHead>
                <TableHead>Score</TableHead>
                <TableHead>Created</TableHead>
                <TableHead>Completed</TableHead>
                <TableHead>Actions</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {queueJobs && queueJobs.length > 0 ? (
                queueJobs.map((job: any) => (
                  <TableRow key={job.id}>
                    <TableCell>
                      <Link href={`/analogs/${job.analogId}`}>
                        <span className="text-blue-600 hover:underline cursor-pointer">
                          Analog #{job.analogId}
                        </span>
                      </Link>
                    </TableCell>
                    <TableCell>
                      <Badge variant="outline">{job.target}</Badge>
                    </TableCell>
                    <TableCell>
                      <Badge variant={getStatusColor(job.status)}>
                        {job.status}
                      </Badge>
                    </TableCell>
                    <TableCell>{job.priority}</TableCell>
                    <TableCell>
                      <div className="flex items-center gap-2">
                        {job.bindingAffinity || '—'}
                        {getAffinityTrend(job.bindingAffinity)}
                      </div>
                    </TableCell>
                    <TableCell>{job.dockingScore || '—'}</TableCell>
                    <TableCell className="text-sm text-muted-foreground">
                      {new Date(job.createdAt).toLocaleDateString()}
                    </TableCell>
                    <TableCell className="text-sm text-muted-foreground">
                      {job.completedAt ? new Date(job.completedAt).toLocaleDateString() : '—'}
                    </TableCell>
                    <TableCell>
                      <Link href={`/analogs/${job.analogId}`}>
                        <Button variant="ghost" size="sm">View</Button>
                      </Link>
                    </TableCell>
                  </TableRow>
                ))
              ) : (
                <TableRow>
                  <TableCell colSpan={9} className="text-center text-muted-foreground py-8">
                    No jobs found matching filters
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          </Table>
        </CardContent>
      </Card>
    </div>
  );
}
