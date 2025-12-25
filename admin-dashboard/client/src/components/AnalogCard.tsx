import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Link } from "wouter";
import { Download, Beaker, FileText, Box } from "lucide-react";
import { useState } from "react";
import { MoleculeViewer3D } from "./MoleculeViewer3D";
import { trpc } from "@/lib/trpc";
import { toast } from "sonner";
import type { AnalogDiscovery } from "../types";

interface AnalogCardProps {
  analog: AnalogDiscovery;
  onRunTest?: (analogId: number, testType: string) => void;
}

export function AnalogCard({ analog, onRunTest }: AnalogCardProps) {
  const [show3D, setShow3D] = useState(false);

  const handleExport = async (format: 'smiles' | 'sdf' | 'pdf') => {
    try {
      let data;
      if (format === 'smiles') {
        data = await trpc.export.smiles.useQuery({ analogId: analog.id }).data;
      } else if (format === 'sdf') {
        data = await trpc.export.sdf.useQuery({ analogId: analog.id }).data;
      } else {
        data = await trpc.export.pdf.useQuery({ analogId: analog.id }).data;
      }

      if (data) {
        const blob = new Blob([data.content], { type: data.mimeType });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = data.filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        toast.success(`Exported ${analog.compoundName} as ${format.toUpperCase()}`);
      }
    } catch (error) {
      toast.error(`Export failed: ${error}`);
    }
  };

  const getConfidenceBadgeColor = (score: number) => {
    if (score >= 90) return "bg-green-100 text-green-800";
    if (score >= 80) return "bg-blue-100 text-blue-800";
    if (score >= 70) return "bg-yellow-100 text-yellow-800";
    return "bg-orange-100 text-orange-800";
  };

  const getPatentStatusColor = (status: string) => {
    if (status === "patent-free") return "bg-green-50 text-green-700 border-green-200";
    if (status === "patent-opportunity") return "bg-yellow-50 text-yellow-700 border-yellow-200";
    return "bg-red-50 text-red-700 border-red-200";
  };

  return (
    <Card className="w-full hover:shadow-lg transition-shadow">
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between">
          <div className="flex-1">
            <Link href={`/admin/analog/${analog.id}`}>
              <CardTitle className="text-xl mb-1 hover:text-primary cursor-pointer transition-colors">{analog.compoundName}</CardTitle>
            </Link>
            <p className="text-sm text-gray-600">Parent: {analog.parentCompound}</p>
          </div>
          <Badge className={getConfidenceBadgeColor(analog.confidenceScore)}>
            {analog.confidenceScore}%
          </Badge>
        </div>
      </CardHeader>

      <CardContent className="space-y-4">
        {/* SMILES Notation */}
        <div className="bg-slate-900 text-slate-100 p-3 rounded font-mono text-xs overflow-x-auto">
          <p className="text-gray-400 text-xs mb-1">SMILES</p>
          {analog.smiles}
        </div>

        {/* Patent Status */}
        <div className={`p-3 rounded border ${getPatentStatusColor(analog.patentStatus)}`}>
          <p className="font-semibold text-sm capitalize">
            {analog.patentStatus === "patent-free"
              ? "Patent-Free"
              : analog.patentStatus === "patent-opportunity"
                ? "Patent Opportunity"
                : "Patented"}
          </p>
        </div>

        {/* Scores Grid */}
        <div className="grid grid-cols-2 gap-3">
          <div className="bg-gray-50 p-3 rounded">
            <p className="text-xs text-gray-600 font-medium">Similarity</p>
            <p className="text-lg font-bold text-gray-900">{analog.similarityScore}%</p>
          </div>
          <div className="bg-gray-50 p-3 rounded">
            <p className="text-xs text-gray-600 font-medium">Safety</p>
            <p className="text-lg font-bold text-gray-900">{analog.safetyScore}/100</p>
          </div>
          <div className="bg-gray-50 p-3 rounded">
            <p className="text-xs text-gray-600 font-medium">Efficacy</p>
            <p className="text-lg font-bold text-gray-900">{analog.efficacyScore}/100</p>
          </div>
          <div className="bg-gray-50 p-3 rounded">
            <p className="text-xs text-gray-600 font-medium">Drug-Likeness</p>
            <p className="text-lg font-bold text-gray-900">{analog.drugLikenessScore}/100</p>
          </div>
        </div>

        {/* Market Value */}
        {analog.marketValue && (
          <div className="border-t pt-3">
            <p className="text-xs text-gray-600 font-medium mb-1">Market Value</p>
            <p className="text-lg font-bold text-blue-600">{analog.marketValue}</p>
          </div>
        )}

        {/* Therapeutic Potential */}
        {analog.therapeuticPotential && (
          <div className="border-t pt-3">
            <p className="text-xs text-gray-600 font-medium mb-1">Therapeutic Potential</p>
            <p className="text-sm text-gray-700">{analog.therapeuticPotential}</p>
          </div>
        )}

        {/* Key Differences */}
        {analog.keyDifferences && (
          <div className="border-t pt-3">
            <p className="text-xs text-gray-600 font-medium mb-1">Key Differences</p>
            <p className="text-sm text-gray-700">{analog.keyDifferences}</p>
          </div>
        )}

        {/* 3D Structure Viewer */}
        {show3D && (
          <div className="border-t pt-3">
            <MoleculeViewer3D
              smiles={analog.smiles}
              compoundName={analog.compoundName}
              width={350}
              height={250}
            />
          </div>
        )}

        {/* Action Buttons */}
        <div className="grid grid-cols-2 gap-2 pt-3 border-t">
          <Button
            size="sm"
            variant={show3D ? "default" : "outline"}
            onClick={() => setShow3D(!show3D)}
            className="col-span-2"
          >
            <Box className="w-4 h-4 mr-1" />
            {show3D ? "Hide" : "View"} 3D Structure
          </Button>
          <Button
            size="sm"
            variant="outline"
            onClick={() => handleExport('smiles')}
          >
            <Download className="w-4 h-4 mr-1" />
            SMILES
          </Button>
          <Button
            size="sm"
            variant="outline"
            onClick={() => handleExport('sdf')}
          >
            <Beaker className="w-4 h-4 mr-1" />
            SDF
          </Button>
        </div>
      </CardContent>
    </Card>
  );
}
