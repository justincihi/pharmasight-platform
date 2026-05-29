import { useEffect, useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { RefreshCw, Trash2, TrendingUp } from 'lucide-react';
import { trpc } from '@/lib/trpc';

interface CacheStats {
  totalEntries: number;
  totalHits: number;
  hitRate: number;
  averageAge: number;
  cacheSize: number;
}

interface CacheStatsDashboardProps {
  autoRefresh?: boolean;
  refreshInterval?: number;
}

/**
 * CacheStatsDashboard
 * Displays cache performance metrics for all analysis types
 * Shows hit rates, memory usage, and cache management options
 */
export default function CacheStatsDashboard({
  autoRefresh = true,
  refreshInterval = 5000,
}: CacheStatsDashboardProps) {
  const [stats, setStats] = useState<Record<string, CacheStats> | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [lastUpdated, setLastUpdated] = useState<Date | null>(null);

  const statsQuery = trpc.cache.stats.useQuery();
  const clearMutation = trpc.cache.clear.useMutation();

  // Fetch cache statistics
  const fetchStats = async () => {
    setIsLoading(true);
    try {
      const result = await statsQuery.refetch();
      if (result.data) {
        setStats(result.data);
        setLastUpdated(new Date());
      }
    } catch (error) {
      console.error('Failed to fetch cache stats:', error);
    } finally {
      setIsLoading(false);
    }
  };

  // Initial fetch
  useEffect(() => {
    fetchStats();
  }, []);

  // Auto-refresh
  useEffect(() => {
    if (!autoRefresh) return;

    const interval = setInterval(() => {
      fetchStats();
    }, refreshInterval);

    return () => clearInterval(interval);
  }, [autoRefresh, refreshInterval]);

  // Handle clear cache
  const handleClearCache = async () => {
    if (!confirm('Are you sure you want to clear all caches? This will remove all cached analysis results.')) {
      return;
    }

    try {
      await clearMutation.mutateAsync();
      setStats(null);
      setLastUpdated(null);
      await fetchStats();
    } catch (error) {
      console.error('Failed to clear cache:', error);
    }
  };

  // Format bytes to human-readable size
  const formatBytes = (bytes: number): string => {
    if (bytes === 0) return '0 B';
    const k = 1024;
    const sizes = ['B', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round((bytes / Math.pow(k, i)) * 100) / 100 + ' ' + sizes[i];
  };

  // Format milliseconds to human-readable time
  const formatTime = (ms: number): string => {
    if (ms < 1000) return Math.round(ms) + ' ms';
    if (ms < 60000) return Math.round(ms / 1000) + ' s';
    return Math.round(ms / 60000) + ' m';
  };

  // Get hit rate color
  const getHitRateColor = (hitRate: number): string => {
    if (hitRate >= 0.7) return 'text-green-600';
    if (hitRate >= 0.4) return 'text-yellow-600';
    return 'text-red-600';
  };

  if (!stats) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Cache Statistics</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-gray-500">Loading cache statistics...</p>
        </CardContent>
      </Card>
    );
  }

  const analysisTypes = ['docking', 'toxicity', 'admet', 'pkpd'] as const;
  const totalCacheSize = Object.values(stats).reduce((sum, s) => sum + s.cacheSize, 0);
  const totalEntries = Object.values(stats).reduce((sum, s) => sum + s.totalEntries, 0);
  const totalHits = Object.values(stats).reduce((sum, s) => sum + s.totalHits, 0);

  return (
    <div className="space-y-6">
      {/* Summary Cards */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium text-gray-600">Total Entries</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{totalEntries}</div>
            <p className="text-xs text-gray-500 mt-1">Cached compounds</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium text-gray-600">Total Hits</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{totalHits}</div>
            <p className="text-xs text-gray-500 mt-1">Cache accesses</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium text-gray-600">Cache Size</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{formatBytes(totalCacheSize)}</div>
            <p className="text-xs text-gray-500 mt-1">Memory usage</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium text-gray-600">Last Updated</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-sm font-mono">
              {lastUpdated ? lastUpdated.toLocaleTimeString() : 'Never'}
            </div>
            <p className="text-xs text-gray-500 mt-1">Auto-refresh: {autoRefresh ? 'On' : 'Off'}</p>
          </CardContent>
        </Card>
      </div>

      {/* Detailed Statistics by Analysis Type */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <CardTitle>Cache Performance by Analysis Type</CardTitle>
            <div className="flex gap-2">
              <Button
                onClick={fetchStats}
                disabled={isLoading}
                variant="outline"
                size="sm"
                className="flex items-center gap-2"
              >
                <RefreshCw className="w-4 h-4" />
                Refresh
              </Button>
              <Button
                onClick={handleClearCache}
                disabled={clearMutation.isPending || totalEntries === 0}
                variant="destructive"
                size="sm"
                className="flex items-center gap-2"
              >
                <Trash2 className="w-4 h-4" />
                Clear All
              </Button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            {analysisTypes.map((type) => {
              const typeStats = stats[type];
              if (!typeStats) return null;

              const hitRatePercent = Math.round(typeStats.hitRate * 100);
              const hitRateColor = getHitRateColor(typeStats.hitRate);

              return (
                <div key={type} className="border rounded-lg p-4">
                  <div className="flex items-center justify-between mb-3">
                    <div className="flex items-center gap-2">
                      <h3 className="font-semibold capitalize">{type} Analysis</h3>
                      <Badge variant="outline">{typeStats.totalEntries} entries</Badge>
                    </div>
                    <div className="flex items-center gap-2">
                      <TrendingUp className="w-4 h-4 text-blue-600" />
                      <span className={`font-semibold ${hitRateColor}`}>{hitRatePercent}%</span>
                    </div>
                  </div>

                  <div className="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm">
                    <div>
                      <p className="text-gray-600">Hit Rate</p>
                      <p className="font-mono font-semibold">{hitRatePercent}%</p>
                    </div>
                    <div>
                      <p className="text-gray-600">Total Hits</p>
                      <p className="font-mono font-semibold">{typeStats.totalHits}</p>
                    </div>
                    <div>
                      <p className="text-gray-600">Avg Age</p>
                      <p className="font-mono font-semibold">{formatTime(typeStats.averageAge)}</p>
                    </div>
                    <div>
                      <p className="text-gray-600">Size</p>
                      <p className="font-mono font-semibold">{formatBytes(typeStats.cacheSize)}</p>
                    </div>
                  </div>

                  {/* Hit rate progress bar */}
                  <div className="mt-3 w-full bg-gray-200 rounded-full h-2">
                    <div
                      className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                      style={{ width: `${hitRatePercent}%` }}
                    />
                  </div>
                </div>
              );
            })}
          </div>
        </CardContent>
      </Card>

      {/* Cache Management Info */}
      <Card className="bg-blue-50 border-blue-200">
        <CardHeader>
          <CardTitle className="text-base">Cache Management</CardTitle>
        </CardHeader>
        <CardContent className="text-sm space-y-2">
          <p>
            <strong>Hit Rate:</strong> Percentage of cache accesses that found cached results. Higher is better.
          </p>
          <p>
            <strong>Average Age:</strong> Average time since cache entries were created. Older entries may be expired.
          </p>
          <p>
            <strong>Cache Size:</strong> Total memory used by cached results. Larger caches may impact performance.
          </p>
          <p>
            <strong>Clear Cache:</strong> Removes all cached results. Use this if you suspect stale data or to free memory.
          </p>
        </CardContent>
      </Card>
    </div>
  );
}
