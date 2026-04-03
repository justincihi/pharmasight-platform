import React, { useState, useEffect } from 'react';
import { trpc } from '@/lib/trpc';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { AlertCircle, CheckCircle2, Clock, Download, Trash2, Play, Pause } from 'lucide-react';
import { formatDistanceToNow } from 'date-fns';

interface BatchJob {
  jobId: string;
  status: 'pending' | 'running' | 'completed' | 'failed' | 'cancelled';
  totalCompounds: number;
  processedCompounds: number;
  successCount: number;
  failureCount: number;
  progress: number;
  createdAt: Date;
  completedAt?: Date;
  statistics?: {
    successRate: number;
    averageAffinity: number;
    bestAffinity: number;
    worstAffinity: number;
  };
}

export function BatchDockingDashboard() {
  const [selectedJob, setSelectedJob] = useState<string | null>(null);
  const [autoRefresh, setAutoRefresh] = useState(true);
  const [refreshInterval, setRefreshInterval] = useState(2000); // 2 seconds

  // Fetch batch jobs
  const { data: jobs = [], isLoading: jobsLoading, refetch: refetchJobs } = trpc.batchDocking.listJobs.useQuery({
    limit: 50,
    offset: 0,
  });

  // Fetch selected job details
  const { data: jobDetailsResponse, refetch: refetchJobDetails } = trpc.batchDocking.getJobDetails.useQuery(
    { jobId: selectedJob || '' },
    { enabled: !!selectedJob }
  );
  
  const jobDetails = jobDetailsResponse?.job;

  // Auto-refresh for running jobs
  useEffect(() => {
    if (!autoRefresh) return;

    const hasRunningJobs = jobs.some((job: any) => job.status === 'running' || job.status === 'pending');
    if (!hasRunningJobs) return;

    const interval = setInterval(() => {
      refetchJobs();
      if (selectedJob) {
        refetchJobDetails();
      }
    }, refreshInterval);

    return () => clearInterval(interval);
  }, [autoRefresh, refreshInterval, jobs, selectedJob, refetchJobs, refetchJobDetails]);

  // Cancel job mutation
  const { mutate: cancelJob, isPending: isCancelling } = trpc.batchDocking.cancelJob.useMutation({
    onSuccess: () => {
      refetchJobs();
      if (selectedJob) {
        refetchJobDetails();
      }
    },
  });

  // Export job handler
  const handleExport = (jobId: string, format: 'csv' | 'json') => {
    // Construct the download URL directly
    const baseUrl = window.location.origin;
    const downloadUrl = `${baseUrl}/api/trpc/batchDocking.exportResults?input=${encodeURIComponent(JSON.stringify({ jobId, format }))}`;
    
    const element = document.createElement('a');
    element.setAttribute('href', downloadUrl);
    element.setAttribute('download', `batch-docking-${jobId}.${format}`);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'completed':
        return 'bg-green-100 text-green-800';
      case 'running':
        return 'bg-blue-100 text-blue-800';
      case 'pending':
        return 'bg-yellow-100 text-yellow-800';
      case 'failed':
        return 'bg-red-100 text-red-800';
      case 'cancelled':
        return 'bg-gray-100 text-gray-800';
      default:
        return 'bg-gray-100 text-gray-800';
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckCircle2 className="w-4 h-4" />;
      case 'running':
        return <Clock className="w-4 h-4 animate-spin" />;
      case 'pending':
        return <Clock className="w-4 h-4" />;
      case 'failed':
        return <AlertCircle className="w-4 h-4" />;
      default:
        return null;
    }
  };

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex justify-between items-center">
        <div>
          <h1 className="text-3xl font-bold">Batch Docking Dashboard</h1>
          <p className="text-gray-600 mt-2">Monitor and manage molecular docking jobs</p>
        </div>
        <div className="flex gap-2">
          <Button
            variant={autoRefresh ? 'default' : 'outline'}
            onClick={() => setAutoRefresh(!autoRefresh)}
          >
            {autoRefresh ? <Pause className="w-4 h-4 mr-2" /> : <Play className="w-4 h-4 mr-2" />}
            {autoRefresh ? 'Auto-refresh ON' : 'Auto-refresh OFF'}
          </Button>
          <Button variant="outline" onClick={() => refetchJobs()}>
            Refresh
          </Button>
        </div>
      </div>

      {/* Main Content */}
      <Tabs defaultValue="jobs" className="w-full">
        <TabsList>
          <TabsTrigger value="jobs">Active Jobs</TabsTrigger>
          <TabsTrigger value="details">Job Details</TabsTrigger>
          <TabsTrigger value="results">Results</TabsTrigger>
        </TabsList>

        {/* Jobs List Tab */}
        <TabsContent value="jobs" className="space-y-4">
          {jobsLoading ? (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">Loading jobs...</p>
              </CardContent>
            </Card>
          ) : jobs.length === 0 ? (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">No batch docking jobs yet. Start a new batch to begin.</p>
              </CardContent>
            </Card>
          ) : (
            <div className="space-y-3">
              {jobs.map((job: any) => (
                <Card
                  key={job.jobId}
                  className={`cursor-pointer transition-all ${
                    selectedJob === job.jobId ? 'ring-2 ring-blue-500' : ''
                  }`}
                  onClick={() => setSelectedJob(job.jobId)}
                >
                  <CardContent className="pt-6">
                    <div className="flex items-start justify-between">
                      <div className="flex-1">
                        <div className="flex items-center gap-2 mb-2">
                          <h3 className="font-semibold">{job.jobId}</h3>
                          <Badge className={getStatusColor(job.status)}>
                            {getStatusIcon(job.status)}
                            <span className="ml-1">{job.status.toUpperCase()}</span>
                          </Badge>
                        </div>

                        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-3">
                          <div>
                            <p className="text-sm text-gray-600">Total Compounds</p>
                            <p className="text-lg font-semibold">{job.totalCompounds}</p>
                          </div>
                          <div>
                            <p className="text-sm text-gray-600">Processed</p>
                            <p className="text-lg font-semibold">{job.processedCompounds}</p>
                          </div>
                          <div>
                            <p className="text-sm text-gray-600">Success</p>
                            <p className="text-lg font-semibold text-green-600">{job.completedCompounds}</p>
                          </div>
                          <div>
                            <p className="text-sm text-gray-600">Failed</p>
                            <p className="text-lg font-semibold text-red-600">{job.failedCompounds}</p>
                          </div>
                        </div>

                        <div className="space-y-2">
                          <div className="flex justify-between text-sm">
                            <span className="text-gray-600">Progress</span>
                            <span className="font-semibold">{Math.round(((job.completedCompounds + job.failedCompounds) / job.totalCompounds) * 100)}%</span>
                          </div>
                          <Progress value={Math.round(((job.completedCompounds + job.failedCompounds) / job.totalCompounds) * 100)} className="h-2" />
                        </div>

                        <p className="text-xs text-gray-500 mt-2">
                          {formatDistanceToNow(new Date(job.createdAt), { addSuffix: true })}
                        </p>
                      </div>

                      <div className="flex gap-2">
                        {job.status === 'running' || job.status === 'pending' ? (
                          <Button
                            variant="destructive"
                            size="sm"
                            onClick={(e) => {
                              e.stopPropagation();
                              cancelJob({ jobId: job.jobId });
                            }}
                            disabled={isCancelling}
                          >
                            Cancel
                          </Button>
                        ) : null}
                        {job.status === 'completed' ? (
                          <Button
                            variant="outline"
                            size="sm"
                            onClick={(e: React.MouseEvent) => {
                              e.stopPropagation();
                              handleExport(job.jobId, 'csv');
                            }}
                            disabled={job.status !== 'completed'}
                          >
                            <Download className="w-4 h-4 mr-1" />
                            Export
                          </Button>
                        ) : null}
                      </div>
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>
          )}
        </TabsContent>

        {/* Job Details Tab */}
        <TabsContent value="details" className="space-y-4">
          {!selectedJob ? (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">Select a job to view details</p>
              </CardContent>
            </Card>
          ) : jobDetails ? (
            <div className="space-y-4">
              {/* Job Info */}
              <Card>
                <CardHeader>
                  <CardTitle>Job Information</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div>
                      <p className="text-sm text-gray-600">Job ID</p>
                      <p className="font-mono text-sm">{jobDetails.jobId}</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Status</p>
                      <Badge className={getStatusColor(jobDetails.status)}>
                        {jobDetails.status.toUpperCase()}
                      </Badge>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Created</p>
                      <p className="text-sm">{new Date(jobDetails.createdAt).toLocaleString()}</p>
                    </div>
                    {jobDetails.completedAt && (
                      <div>
                        <p className="text-sm text-gray-600">Completed</p>
                        <p className="text-sm">{new Date(jobDetails.completedAt).toLocaleString()}</p>
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>

              {/* Statistics */}
              {jobDetails.statistics && (
                <Card>
                  <CardHeader>
                    <CardTitle>Statistics</CardTitle>
                  </CardHeader>
                  <CardContent className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div>
                      <p className="text-sm text-gray-600">Success Rate</p>
                      <p className="text-2xl font-bold">{jobDetails.statistics.successRate}%</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Avg Affinity</p>
                      <p className="text-2xl font-bold">{jobDetails.statistics.averageAffinity.toFixed(2)}</p>
                      <p className="text-xs text-gray-500">kcal/mol</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Best Affinity</p>
                      <p className="text-2xl font-bold text-green-600">
                        {jobDetails.statistics.bestAffinity.toFixed(2)}
                      </p>
                      <p className="text-xs text-gray-500">kcal/mol</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Worst Affinity</p>
                      <p className="text-2xl font-bold text-red-600">
                        {jobDetails.statistics.worstAffinity.toFixed(2)}
                      </p>
                      <p className="text-xs text-gray-500">kcal/mol</p>
                    </div>
                  </CardContent>
                </Card>
              )}

              {/* Progress */}
              <Card>
                <CardHeader>
                  <CardTitle>Progress</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div>
                      <p className="text-sm text-gray-600">Total</p>
                      <p className="text-2xl font-bold">{jobDetails.totalCompounds}</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Processed</p>
                      <p className="text-2xl font-bold">{jobDetails.processedCompounds}</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Success</p>
                      <p className="text-2xl font-bold text-green-600">{jobDetails.completedCompounds}</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">Failed</p>
                      <p className="text-2xl font-bold text-red-600">{jobDetails.failedCompounds}</p>
                    </div>
                  </div>

                  <div className="space-y-2">
                    <div className="flex justify-between text-sm">
                      <span className="text-gray-600">Overall Progress</span>
                      <span className="font-semibold">{jobDetailsResponse?.progressPercentage}%</span>
                    </div>
                    <Progress value={jobDetailsResponse?.progressPercentage || 0} className="h-3" />
                  </div>
                </CardContent>
              </Card>
            </div>
          ) : (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">Loading job details...</p>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        {/* Results Tab */}
        <TabsContent value="results" className="space-y-4">
          {!selectedJob ? (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">Select a job to view results</p>
              </CardContent>
            </Card>
          ) : jobDetailsResponse?.results && jobDetailsResponse.results.length > 0 ? (
            <div className="space-y-3">
              {jobDetailsResponse.results.map((result: any, idx: number) => (
                <Card key={idx}>
                  <CardContent className="pt-6">
                    <div className="flex justify-between items-start">
                      <div className="flex-1">
                        <h4 className="font-semibold mb-2">{result.compoundName}</h4>
                        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                          <div>
                            <p className="text-gray-600">Status</p>
                            <Badge className={result.status === 'completed' ? 'bg-green-100 text-green-800' : 'bg-red-100 text-red-800'}>
                              {result.status}
                            </Badge>
                          </div>
                          {result.bindingAffinity && (
                            <div>
                              <p className="text-gray-600">Binding Affinity</p>
                              <p className="font-semibold">{result.bindingAffinity} kcal/mol</p>
                            </div>
                          )}
                          {result.numPoses && (
                            <div>
                              <p className="text-gray-600">Poses</p>
                              <p className="font-semibold">{result.numPoses}</p>
                            </div>
                          )}
                          {result.errorMessage && (
                            <div>
                              <p className="text-gray-600">Error</p>
                              <p className="text-red-600 text-xs">{result.errorMessage}</p>
                            </div>
                          )}
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>
          ) : (
            <Card>
              <CardContent className="pt-6">
                <p className="text-gray-600">No results available for this job</p>
              </CardContent>
            </Card>
          )}
        </TabsContent>
      </Tabs>
    </div>
  );
}
