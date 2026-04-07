import { useState } from "react";
import DashboardLayout from "@/components/DashboardLayout";
import { AnalogCard } from "@/components/AnalogCard";
import SDFUploader from "@/components/SDFUploader";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Loader2, Search, Filter, FlaskConical } from "lucide-react";
import { trpc } from "@/lib/trpc";
import type { AnalogDiscovery } from "@/types";
import { BatchAnalysisModal } from "@/components/BatchAnalysisModal";
import { Checkbox } from "@/components/ui/checkbox";

export default function AdminDashboard() {
  const [searchQuery, setSearchQuery] = useState("");
  const [patentFilter, setPatentFilter] = useState<string | undefined>();
  const [confidenceFilter, setConfidenceFilter] = useState<number | undefined>();
  const [page, setPage] = useState(0);
  const [selectedAnalogs, setSelectedAnalogs] = useState<Set<number>>(new Set());
  const [showBatchAnalysis, setShowBatchAnalysis] = useState(false);

  // Fetch analogs
  const { data: analogs, isLoading, error } = trpc.analog.list.useQuery({
    limit: 12,
    offset: page * 12,
    patentStatus: patentFilter,
    minConfidence: confidenceFilter,
  });

  // Search functionality
  const { data: searchResults } = trpc.analog.search.useQuery(
    { query: searchQuery },
    { enabled: searchQuery.length > 2 }
  );

  const displayAnalogs = searchQuery.length > 2 ? searchResults : analogs;

  const toggleAnalogSelection = (analogId: number) => {
    setSelectedAnalogs(prev => {
      const newSet = new Set(prev);
      if (newSet.has(analogId)) {
        newSet.delete(analogId);
      } else {
        newSet.add(analogId);
      }
      return newSet;
    });
  };

  const toggleSelectAll = () => {
    if (selectedAnalogs.size === displayAnalogs?.length) {
      setSelectedAnalogs(new Set());
    } else {
      setSelectedAnalogs(new Set(displayAnalogs?.map((a: AnalogDiscovery) => a.id) || []));
    }
  };

  const getSelectedAnalogsData = () => {
    if (!displayAnalogs) return [];
    return displayAnalogs
      .filter((a: AnalogDiscovery) => selectedAnalogs.has(a.id))
      .map((a: AnalogDiscovery) => ({
        id: a.id,
        compoundName: a.compoundName,
        smiles: a.smiles
      }));
  };

  const handleRunTest = (analogId: number, testType: string) => {
    console.log(`Running ${testType} test on analog ${analogId}`);
    // TODO: Implement test running logic
  };

  return (
    <DashboardLayout>
      <div className="space-y-6">
        {/* Header */}
        <div>
          <h1 className="text-3xl font-bold text-gray-900">Analog Discoveries</h1>
          <p className="text-gray-600 mt-1">
            Manage and analyze pharmaceutical analog compounds
          </p>
        </div>
        
        {/* SDF Uploader */}
        <SDFUploader />

        {/* Search and Filters */}
        <div className="bg-white p-4 rounded-lg border border-gray-200 space-y-4">
          <div className="flex gap-2">
            <div className="flex-1 relative">
              <Search className="absolute left-3 top-3 w-4 h-4 text-gray-400" />
              <Input
                placeholder="Search by compound name or SMILES..."
                value={searchQuery}
                onChange={(e) => {
                  setSearchQuery(e.target.value);
                  setPage(0);
                }}
                className="pl-10"
              />
            </div>
            {selectedAnalogs.size > 0 && (
              <Button 
                variant="default" 
                size="sm"
                onClick={() => setShowBatchAnalysis(true)}
              >
                <FlaskConical className="w-4 h-4 mr-2" />
                Analyze Selected ({selectedAnalogs.size})
              </Button>
            )}
            <Button variant="outline" size="sm" onClick={toggleSelectAll}>
              {selectedAnalogs.size === displayAnalogs?.length ? 'Deselect All' : 'Select All'}
            </Button>
            <Button variant="outline" size="sm">
              <Filter className="w-4 h-4 mr-2" />
              Advanced
            </Button>
          </div>

          <div className="grid grid-cols-2 gap-4">
            <div>
              <label className="text-sm font-medium text-gray-700 block mb-2">
                Patent Status
              </label>
              <Select value={patentFilter || "all"} onValueChange={(val) => {
                setPatentFilter(val === "all" ? undefined : val);
                setPage(0);
              }}>
                <SelectTrigger>
                  <SelectValue placeholder="All statuses" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">All Statuses</SelectItem>
                  <SelectItem value="patent-free">Patent-Free</SelectItem>
                  <SelectItem value="patent-opportunity">Patent Opportunity</SelectItem>
                  <SelectItem value="patented">Patented</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div>
              <label className="text-sm font-medium text-gray-700 block mb-2">
                Minimum Confidence
              </label>
              <Select value={confidenceFilter?.toString() || "all"} onValueChange={(val) => {
                setConfidenceFilter(val === "all" ? undefined : parseInt(val));
                setPage(0);
              }}>
                <SelectTrigger>
                  <SelectValue placeholder="Any confidence" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">Any Confidence</SelectItem>
                  <SelectItem value="85">85%+</SelectItem>
                  <SelectItem value="90">90%+</SelectItem>
                  <SelectItem value="95">95%+</SelectItem>
                </SelectContent>
              </Select>
            </div>
          </div>
        </div>

        {/* Loading State */}
        {isLoading && (
          <div className="flex items-center justify-center py-12">
            <Loader2 className="w-8 h-8 animate-spin text-blue-600" />
          </div>
        )}

        {/* Error State */}
        {error && (
          <div className="bg-red-50 border border-red-200 rounded-lg p-4">
            <p className="text-red-800 font-medium">Error loading analogs</p>
            <p className="text-red-700 text-sm mt-1">{error.message}</p>
          </div>
        )}

        {/* Analogs Grid */}
        {!isLoading && displayAnalogs && displayAnalogs.length > 0 && (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            {displayAnalogs.map((analog: AnalogDiscovery) => (
              <div key={analog.id} className="relative">
                <div className="absolute top-2 left-2 z-10">
                  <Checkbox
                    checked={selectedAnalogs.has(analog.id)}
                    onCheckedChange={() => toggleAnalogSelection(analog.id)}
                    className="bg-white border-2"
                  />
                </div>
                <AnalogCard
                  analog={analog}
                  onRunTest={handleRunTest}
                />
              </div>
            ))}
          </div>
        )}

        {/* Batch Analysis Modal */}
        <BatchAnalysisModal
          open={showBatchAnalysis}
          onOpenChange={setShowBatchAnalysis}
          selectedAnalogs={getSelectedAnalogsData()}
        />

        {/* Empty State */}
        {!isLoading && displayAnalogs && displayAnalogs.length === 0 && (
          <div className="text-center py-12">
            <p className="text-gray-600 mb-4">No analogs found matching your criteria</p>
            <Button
              variant="outline"
              onClick={() => {
                setSearchQuery("");
                setPatentFilter(undefined);
                setConfidenceFilter(undefined);
                setPage(0);
              }}
            >
              Clear Filters
            </Button>
          </div>
        )}

        {/* Pagination */}
        {!isLoading && displayAnalogs && displayAnalogs.length > 0 && (
          <div className="flex justify-center gap-2 pt-4">
            <Button
              variant="outline"
              onClick={() => setPage(Math.max(0, page - 1))}
              disabled={page === 0}
            >
              Previous
            </Button>
            <span className="px-4 py-2 text-sm text-gray-600">
              Page {page + 1}
            </span>
            <Button
              variant="outline"
              onClick={() => setPage(page + 1)}
              disabled={!displayAnalogs || displayAnalogs.length < 12}
            >
              Next
            </Button>
          </div>
        )}
      </div>
    </DashboardLayout>
  );
}
