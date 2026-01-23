import DashboardLayout from "@/components/DashboardLayout";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Loader2 } from "lucide-react";
import { LineChart, Line, PieChart, Pie, Cell, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, BarChart, Bar } from "recharts";
import { trpc } from "@/lib/trpc";

export default function Analytics() {
  const { data: stats, isLoading: statsLoading } = trpc.analytics.getStats.useQuery();
  const { data: timeline, isLoading: timelineLoading } = trpc.analytics.getTimeline.useQuery({ days: 30 });

  // Prepare chart data
  const patentData = stats ? [
    { name: "Patent-Free", value: stats.patentFreeCount, color: "#10b981" },
    { name: "Patented", value: stats.patentedCount, color: "#ef4444" },
  ] : [];

  const confidenceData = stats ? [
    { name: "85%+ Confidence", value: stats.highConfidenceCount, color: "#3b82f6" },
    { name: "Below 85%", value: stats.totalDiscovered - stats.highConfidenceCount, color: "#d1d5db" },
  ] : [];

  const timelineData = timeline
    ? timeline.map((item: any, idx: number) => ({
        date: new Date(item.discoveredAt).toLocaleDateString(),
        count: idx + 1,
      }))
    : [];

  if (statsLoading || timelineLoading) {
    return (
      <DashboardLayout>
        <div className="flex items-center justify-center py-12">
          <Loader2 className="w-8 h-8 animate-spin text-blue-600" />
        </div>
      </DashboardLayout>
    );
  }

  return (
    <DashboardLayout>
      <div className="space-y-6">
        {/* Header */}
        <div>
          <h1 className="text-3xl font-bold text-gray-900">Analytics Dashboard</h1>
          <p className="text-gray-600 mt-1">
            Discovery statistics and performance metrics
          </p>
        </div>

        {/* Key Metrics */}
        {stats && (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium text-gray-600">
                  Total Discovered
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-gray-900">
                  {stats.totalDiscovered}
                </p>
                <p className="text-xs text-gray-500 mt-1">analogs</p>
              </CardContent>
            </Card>

            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium text-gray-600">
                  High Confidence (85%+)
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-green-600">
                  {stats.highConfidencePercentage}%
                </p>
                <p className="text-xs text-gray-500 mt-1">
                  {stats.highConfidenceCount} compounds
                </p>
              </CardContent>
            </Card>

            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium text-gray-600">
                  Patent-Free
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-blue-600">
                  {stats.patentFreePercentage}%
                </p>
                <p className="text-xs text-gray-500 mt-1">
                  {stats.patentFreeCount} compounds
                </p>
              </CardContent>
            </Card>

            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium text-gray-600">
                  Avg Confidence
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-3xl font-bold text-purple-600">
                  {stats.averageConfidence}%
                </p>
                <p className="text-xs text-gray-500 mt-1">across all analogs</p>
              </CardContent>
            </Card>
          </div>
        )}

        {/* Charts */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Patent Status Pie Chart */}
          <Card>
            <CardHeader>
              <CardTitle>Patent Status Distribution</CardTitle>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={300}>
                <PieChart>
                  <Pie
                    data={patentData}
                    cx="50%"
                    cy="50%"
                    labelLine={false}
                    label={({ name, value }) => `${name}: ${value}`}
                    outerRadius={80}
                    fill="#8884d8"
                    dataKey="value"
                  >
                    {patentData.map((entry: any, index: number) => (
                      <Cell key={`cell-${index}`} fill={entry.color} />
                    ))}
                  </Pie>
                  <Tooltip />
                </PieChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>

          {/* Confidence Level Pie Chart */}
          <Card>
            <CardHeader>
              <CardTitle>Confidence Level Distribution</CardTitle>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={300}>
                <PieChart>
                  <Pie
                    data={confidenceData}
                    cx="50%"
                    cy="50%"
                    labelLine={false}
                    label={({ name, value }) => `${name}: ${value}`}
                    outerRadius={80}
                    fill="#8884d8"
                    dataKey="value"
                  >
                    {confidenceData.map((entry: any, index: number) => (
                      <Cell key={`cell-${index}`} fill={entry.color} />
                    ))}
                  </Pie>
                  <Tooltip />
                </PieChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>
        </div>

        {/* Discovery Timeline */}
        {timelineData.length > 0 && (
          <Card>
            <CardHeader>
              <CardTitle>Discovery Timeline (Last 30 Days)</CardTitle>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={300}>
                <LineChart data={timelineData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="date" />
                  <YAxis />
                  <Tooltip />
                  <Legend />
                  <Line
                    type="monotone"
                    dataKey="count"
                    stroke="#3b82f6"
                    dot={{ fill: "#3b82f6", r: 4 }}
                    activeDot={{ r: 6 }}
                    name="Cumulative Discoveries"
                  />
                </LineChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>
        )}
      </div>
    </DashboardLayout>
  );
}
