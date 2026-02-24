import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Checkbox } from "@/components/ui/checkbox";
import { Badge } from "@/components/ui/badge";
import { toast } from "sonner";
import { Download, Filter, CheckCircle, XCircle } from "lucide-react";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Input } from "@/components/ui/input";

export default function BatchOperations() {
  const [selectedAnalogs, setSelectedAnalogs] = useState<number[]>([]);
  const [selectAll, setSelectAll] = useState(false);
  const [patentFilter, setPatentFilter] = useState<string>("all");
  const [minConfidence, setMinConfidence] = useState<number>(0);
  const [parentFilter, setParentFilter] = useState<string>("all");

  // Fetch all analogs with filters
  const { data: analogs, isLoading } = trpc.analog.list.useQuery({
    limit: 1000,
    offset: 0,
    patentStatus: patentFilter !== "all" ? patentFilter : undefined,
    minConfidence: minConfidence > 0 ? minConfidence : undefined,
  });

  // Get unique parent compounds for filter
  const parentCompounds: string[] = Array.from(
    new Set(analogs?.analogs?.map((a: any) => a.parentCompound) || [])
  ) as string[];

  // Filter by parent compound on client side
  const filteredAnalogs = analogs?.analogs?.filter((a: any) => 
    parentFilter === "all" || a.parentCompound === parentFilter
  ) || [];

  const handleSelectAll = () => {
    if (selectAll) {
      setSelectedAnalogs([]);
    } else {
      setSelectedAnalogs(filteredAnalogs?.map((a: any) => a.id) || []);
    }
    setSelectAll(!selectAll);
  };

  const handleToggleAnalog = (id: number) => {
    if (selectedAnalogs.includes(id)) {
      setSelectedAnalogs(selectedAnalogs.filter(aid => aid !== id));
    } else {
      setSelectedAnalogs([...selectedAnalogs, id]);
    }
  };

  const bulkApproveMutation = trpc.analog.bulkApprove.useMutation({
    onSuccess: (data) => {
      toast.success(`Successfully approved ${data.count} analog(s)`);
      setSelectedAnalogs([]);
      setSelectAll(false);
    },
    onError: (error) => {
      toast.error(`Failed to approve analogs: ${error.message}`);
    },
  });

  const bulkRejectMutation = trpc.analog.bulkReject.useMutation({
    onSuccess: (data) => {
      toast.success(`Successfully rejected ${data.count} analog(s)`);
      setSelectedAnalogs([]);
      setSelectAll(false);
    },
    onError: (error) => {
      toast.error(`Failed to reject analogs: ${error.message}`);
    },
  });

  const handleBulkApprove = () => {
    if (selectedAnalogs.length === 0) {
      toast.error("Please select analogs to approve");
      return;
    }
    bulkApproveMutation.mutate({ analogIds: selectedAnalogs });
  };

  const handleBulkReject = () => {
    if (selectedAnalogs.length === 0) {
      toast.error("Please select analogs to reject");
      return;
    }
    bulkRejectMutation.mutate({ analogIds: selectedAnalogs });
  };

  const handleExportSelected = () => {
    if (selectedAnalogs.length === 0) {
    toast.error("Please select analogs to export");
      return;
    }

    const selectedData = analogs?.analogs.filter((a: any) => selectedAnalogs.includes(a.id));
    const csv = [
      ["Compound ID", "Name", "Parent", "SMILES", "Confidence", "Similarity", "Safety", "Efficacy", "Patent Status"].join(","),
      ...(selectedData || []).map((a: any) => [
        a.compoundId,
        a.compoundName,
        a.parentCompound,
        a.smiles,
        a.confidenceScore,
        a.similarityScore,
        a.safetyScore,
        a.efficacyScore,
        a.patentStatus
      ].join(","))
    ].join("\n");

    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `pharmasight-selected-analogs-${new Date().toISOString().split('T')[0]}.csv`;
    a.click();
    URL.revokeObjectURL(url);

      toast.success(`Exported ${selectedAnalogs.length} analogs to CSV`);
  };

  const handleExportAll = () => {
    if (!analogs?.analogs || analogs.analogs.length === 0) {
    toast.error("No analogs available to export");
      return;
    }

    const csv = [
      ["Compound ID", "Name", "Parent", "SMILES", "Confidence", "Similarity", "Safety", "Efficacy", "Patent Status"].join(","),
      ...analogs.analogs.map((a: any) => [
        a.compoundId,
        a.compoundName,
        a.parentCompound,
        a.smiles,
        a.confidenceScore,
        a.similarityScore,
        a.safetyScore,
        a.efficacyScore,
        a.patentStatus
      ].join(","))
    ].join("\n");

    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `pharmasight-all-analogs-${new Date().toISOString().split('T')[0]}.csv`;
    a.click();
    URL.revokeObjectURL(url);

    toast.success(`Exported all ${analogs.analogs.length} analogs to CSV`);
  };

  if (isLoading) {
    return (
      <div className="container py-8">
        <div className="flex items-center justify-center h-64">
          <p className="text-muted-foreground">Loading analogs...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="container py-8">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Batch Operations</h1>
        <p className="text-muted-foreground">
          Select and export multiple analogs at once for further analysis
        </p>
      </div>

      {/* Filters */}
      <Card className="mb-6">
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Filter className="w-5 h-5" />
            Filters
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div>
              <label className="text-sm font-medium mb-2 block">Patent Status</label>
              <Select value={patentFilter} onValueChange={setPatentFilter}>
                <SelectTrigger>
                  <SelectValue placeholder="All statuses" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">All Statuses</SelectItem>
                  <SelectItem value="patent_free">Patent-Free</SelectItem>
                  <SelectItem value="patent_opportunity">Patent Opportunity</SelectItem>
                  <SelectItem value="patented">Patented</SelectItem>
                </SelectContent>
              </Select>
            </div>
            <div>
              <label className="text-sm font-medium mb-2 block">Parent Compound</label>
              <Select value={parentFilter} onValueChange={setParentFilter}>
                <SelectTrigger>
                  <SelectValue placeholder="All compounds" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">All Compounds</SelectItem>
                  {parentCompounds.map((parent: string) => (
                    <SelectItem key={parent} value={parent}>{parent}</SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>
            <div>
              <label className="text-sm font-medium mb-2 block">Min Confidence (%)</label>
              <Input
                type="number"
                min="0"
                max="100"
                value={minConfidence}
                onChange={(e) => setMinConfidence(parseInt(e.target.value) || 0)}
                placeholder="0"
              />
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Action Bar */}
      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Bulk Actions</CardTitle>
          <CardDescription>
            {selectedAnalogs.length} of {filteredAnalogs?.length || 0} analogs selected
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-2">
            <Button
              onClick={() => handleBulkApprove()}
              disabled={selectedAnalogs.length === 0}
              variant="default"
              className="bg-green-600 hover:bg-green-700"
            >
              <CheckCircle className="mr-2 h-4 w-4" />
              Approve Selected ({selectedAnalogs.length})
            </Button>
            <Button
              onClick={() => handleBulkReject()}
              disabled={selectedAnalogs.length === 0}
              variant="destructive"
            >
              <XCircle className="mr-2 h-4 w-4" />
              Reject Selected ({selectedAnalogs.length})
            </Button>
            <Button
              onClick={handleExportSelected}
              disabled={selectedAnalogs.length === 0}
              variant="outline"
            >
              <Download className="mr-2 h-4 w-4" />
              Export Selected ({selectedAnalogs.length})
            </Button>
            <Button
              onClick={handleExportAll}
              variant="outline"
            >
              <Download className="mr-2 h-4 w-4" />
              Export All ({analogs?.analogs?.length || 0})
            </Button>
          </div>
        </CardContent>
      </Card>

      {/* Analog List */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <CardTitle>All Analogs ({filteredAnalogs?.length || 0})</CardTitle>
            <div className="flex items-center gap-2">
              <Checkbox
                id="select-all"
                checked={selectAll}
                onCheckedChange={handleSelectAll}
              />
              <label htmlFor="select-all" className="text-sm font-medium cursor-pointer">
                Select All
              </label>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <div className="space-y-2 max-h-[600px] overflow-y-auto">
            {filteredAnalogs?.map((analog: any) => (
              <div
                key={analog.id}
                className={`flex items-center gap-4 p-4 rounded-lg border ${
                  selectedAnalogs.includes(analog.id) ? "bg-accent" : "bg-background"
                }`}
              >
                <Checkbox
                  checked={selectedAnalogs.includes(analog.id)}
                  onCheckedChange={() => handleToggleAnalog(analog.id)}
                />
                <div className="flex-1 grid grid-cols-1 md:grid-cols-4 gap-4">
                  <div>
                    <p className="font-semibold">{analog.compoundId}</p>
                    <p className="text-sm text-muted-foreground">{analog.compoundName}</p>
                  </div>
                  <div>
                    <p className="text-sm text-muted-foreground">Parent</p>
                    <p className="font-medium">{analog.parentCompound}</p>
                  </div>
                  <div className="flex gap-2">
                    <Badge variant={analog.patentStatus === "patent-free" ? "default" : "secondary"}>
                      {analog.patentStatus}
                    </Badge>
                    <Badge variant="outline">{analog.confidenceScore}%</Badge>
                  </div>
                  <div className="text-right">
                    <p className="text-sm text-muted-foreground">Safety: {analog.safetyScore}/100</p>
                    <p className="text-sm text-muted-foreground">Efficacy: {analog.efficacyScore}/100</p>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
