import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Checkbox } from "@/components/ui/checkbox";
import { Badge } from "@/components/ui/badge";
import { toast } from "sonner";
import { Download } from "lucide-react";

export default function BatchOperations() {
  const [selectedAnalogs, setSelectedAnalogs] = useState<number[]>([]);
  const [selectAll, setSelectAll] = useState(false);

  // Fetch all analogs
  const { data: analogs, isLoading } = trpc.analog.list.useQuery({
    limit: 1000, // Get all analogs
    offset: 0
  });

  const handleSelectAll = () => {
    if (selectAll) {
      setSelectedAnalogs([]);
    } else {
      setSelectedAnalogs(analogs?.analogs.map((a: any) => a.id) || []);
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

      {/* Action Bar */}
      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Bulk Actions</CardTitle>
          <CardDescription>
            {selectedAnalogs.length} of {analogs?.analogs?.length || 0} analogs selected
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-2">
            <Button
              onClick={handleExportSelected}
              disabled={selectedAnalogs.length === 0}
              variant="default"
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
            <CardTitle>All Analogs ({analogs?.analogs?.length || 0})</CardTitle>
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
            {analogs?.analogs?.map((analog: any) => (
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
