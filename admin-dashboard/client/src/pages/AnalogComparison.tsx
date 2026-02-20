import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Badge } from "@/components/ui/badge";
import { ArrowLeft, TrendingUp, TrendingDown, Minus } from "lucide-react";
import { useLocation } from "wouter";

export default function AnalogComparison() {
  const [, setLocation] = useLocation();
  const [analog1Id, setAnalog1Id] = useState<number | null>(null);
  const [analog2Id, setAnalog2Id] = useState<number | null>(null);

  const { data: analogsList } = trpc.analog.list.useQuery({ limit: 1000, offset: 0 });
  const { data: analog1 } = trpc.analog.getById.useQuery(
    { id: analog1Id || 0 },
    { enabled: !!analog1Id }
  );
  const { data: analog2 } = trpc.analog.getById.useQuery(
    { id: analog2Id || 0 },
    { enabled: !!analog2Id }
  );

  const analogs = analogsList?.analogs || [];

  const compareValue = (val1: number, val2: number) => {
    if (val1 > val2) return <TrendingUp className="w-4 h-4 text-green-600" />;
    if (val1 < val2) return <TrendingDown className="w-4 h-4 text-red-600" />;
    return <Minus className="w-4 h-4 text-gray-400" />;
  };

  const renderPropertyBar = (value: number, color: string) => (
    <div className="w-full bg-muted rounded-full h-3">
      <div
        className={`h-3 rounded-full ${color}`}
        style={{ width: `${value}%` }}
      />
    </div>
  );

  return (
    <div className="container py-8">
      <Button
        variant="ghost"
        onClick={() => setLocation("/admin/dashboard")}
        className="mb-4"
      >
        <ArrowLeft className="mr-2 h-4 w-4" />
        Back to Dashboard
      </Button>

      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Analog Comparison</h1>
        <p className="text-muted-foreground">
          Compare two analogs side-by-side with visual property charts
        </p>
      </div>

      {/* Analog Selectors */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-8">
        <Card>
          <CardHeader>
            <CardTitle>Select Analog 1</CardTitle>
          </CardHeader>
          <CardContent>
            <Select
              value={analog1Id?.toString() || ""}
              onValueChange={(val) => setAnalog1Id(parseInt(val))}
            >
              <SelectTrigger>
                <SelectValue placeholder="Choose first analog" />
              </SelectTrigger>
              <SelectContent>
                {analogs.map((a: any) => (
                  <SelectItem key={a.id} value={a.id.toString()}>
                    {a.compoundName}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Select Analog 2</CardTitle>
          </CardHeader>
          <CardContent>
            <Select
              value={analog2Id?.toString() || ""}
              onValueChange={(val) => setAnalog2Id(parseInt(val))}
            >
              <SelectTrigger>
                <SelectValue placeholder="Choose second analog" />
              </SelectTrigger>
              <SelectContent>
                {analogs.map((a: any) => (
                  <SelectItem key={a.id} value={a.id.toString()}>
                    {a.compoundName}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </CardContent>
        </Card>
      </div>

      {/* Comparison View */}
      {analog1 && analog2 ? (
        <div className="space-y-6">
          {/* Basic Info Comparison */}
          <Card>
            <CardHeader>
              <CardTitle>Basic Information</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div>
                  <h3 className="font-bold text-lg mb-2">{analog1.compoundName}</h3>
                  <p className="text-sm text-muted-foreground mb-1">
                    <span className="font-medium">Compound ID:</span> {analog1.compoundId}
                  </p>
                  <p className="text-sm text-muted-foreground mb-1">
                    <span className="font-medium">Parent:</span> {analog1.parentCompound}
                  </p>
                  <Badge
                    variant={analog1.patentStatus === 'patent_free' ? 'default' : 'secondary'}
                    className="mt-2"
                  >
                    {analog1.patentStatus}
                  </Badge>
                </div>
                <div>
                  <h3 className="font-bold text-lg mb-2">{analog2.compoundName}</h3>
                  <p className="text-sm text-muted-foreground mb-1">
                    <span className="font-medium">Compound ID:</span> {analog2.compoundId}
                  </p>
                  <p className="text-sm text-muted-foreground mb-1">
                    <span className="font-medium">Parent:</span> {analog2.parentCompound}
                  </p>
                  <Badge
                    variant={analog2.patentStatus === 'patent_free' ? 'default' : 'secondary'}
                    className="mt-2"
                  >
                    {analog2.patentStatus}
                  </Badge>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Property Scores Comparison */}
          <Card>
            <CardHeader>
              <CardTitle>Property Scores Comparison</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-6">
                {/* Confidence Score */}
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">Confidence Score</span>
                    {compareValue(analog1.confidenceScore, analog2.confidenceScore)}
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog1.compoundName}</span>
                        <span className="text-sm font-bold">{analog1.confidenceScore}%</span>
                      </div>
                      {renderPropertyBar(analog1.confidenceScore, "bg-primary")}
                    </div>
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog2.compoundName}</span>
                        <span className="text-sm font-bold">{analog2.confidenceScore}%</span>
                      </div>
                      {renderPropertyBar(analog2.confidenceScore, "bg-primary")}
                    </div>
                  </div>
                </div>

                {/* Similarity Score */}
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">Similarity Score</span>
                    {compareValue(analog1.similarityScore, analog2.similarityScore)}
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog1.compoundName}</span>
                        <span className="text-sm font-bold">{analog1.similarityScore}%</span>
                      </div>
                      {renderPropertyBar(analog1.similarityScore, "bg-blue-500")}
                    </div>
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog2.compoundName}</span>
                        <span className="text-sm font-bold">{analog2.similarityScore}%</span>
                      </div>
                      {renderPropertyBar(analog2.similarityScore, "bg-blue-500")}
                    </div>
                  </div>
                </div>

                {/* Safety Score */}
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">Safety Score</span>
                    {compareValue(analog1.safetyScore, analog2.safetyScore)}
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog1.compoundName}</span>
                        <span className="text-sm font-bold">{analog1.safetyScore}/100</span>
                      </div>
                      {renderPropertyBar(analog1.safetyScore, "bg-green-500")}
                    </div>
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog2.compoundName}</span>
                        <span className="text-sm font-bold">{analog2.safetyScore}/100</span>
                      </div>
                      {renderPropertyBar(analog2.safetyScore, "bg-green-500")}
                    </div>
                  </div>
                </div>

                {/* Efficacy Score */}
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">Efficacy Score</span>
                    {compareValue(analog1.efficacyScore, analog2.efficacyScore)}
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog1.compoundName}</span>
                        <span className="text-sm font-bold">{analog1.efficacyScore}/100</span>
                      </div>
                      {renderPropertyBar(analog1.efficacyScore, "bg-purple-500")}
                    </div>
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog2.compoundName}</span>
                        <span className="text-sm font-bold">{analog2.efficacyScore}/100</span>
                      </div>
                      {renderPropertyBar(analog2.efficacyScore, "bg-purple-500")}
                    </div>
                  </div>
                </div>

                {/* Drug-Likeness Score */}
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">Drug-Likeness Score</span>
                    {compareValue(analog1.drugLikenessScore, analog2.drugLikenessScore)}
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog1.compoundName}</span>
                        <span className="text-sm font-bold">{analog1.drugLikenessScore}/100</span>
                      </div>
                      {renderPropertyBar(analog1.drugLikenessScore, "bg-orange-500")}
                    </div>
                    <div>
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-sm">{analog2.compoundName}</span>
                        <span className="text-sm font-bold">{analog2.drugLikenessScore}/100</span>
                      </div>
                      {renderPropertyBar(analog2.drugLikenessScore, "bg-orange-500")}
                    </div>
                  </div>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Chemical Structure Comparison */}
          <Card>
            <CardHeader>
              <CardTitle>Chemical Structures (SMILES)</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div>
                  <p className="text-sm font-medium mb-2">{analog1.compoundName}</p>
                  <div className="bg-slate-900 text-slate-100 p-4 rounded font-mono text-xs overflow-x-auto">
                    {analog1.smiles}
                  </div>
                </div>
                <div>
                  <p className="text-sm font-medium mb-2">{analog2.compoundName}</p>
                  <div className="bg-slate-900 text-slate-100 p-4 rounded font-mono text-xs overflow-x-auto">
                    {analog2.smiles}
                  </div>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Market Value Comparison */}
          {(analog1.marketValue || analog2.marketValue) && (
            <Card>
              <CardHeader>
                <CardTitle>Market Value</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  <div>
                    <p className="text-sm text-muted-foreground mb-1">{analog1.compoundName}</p>
                    <p className="text-2xl font-bold text-primary">
                      {analog1.marketValue || "N/A"}
                    </p>
                  </div>
                  <div>
                    <p className="text-sm text-muted-foreground mb-1">{analog2.compoundName}</p>
                    <p className="text-2xl font-bold text-primary">
                      {analog2.marketValue || "N/A"}
                    </p>
                  </div>
                </div>
              </CardContent>
            </Card>
          )}
        </div>
      ) : (
        <Card>
          <CardContent className="py-12 text-center text-muted-foreground">
            <p>Select two analogs above to compare their properties</p>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
