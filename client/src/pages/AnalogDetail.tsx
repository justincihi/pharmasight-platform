import { useParams, useLocation } from "wouter";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { ArrowLeft, Download, Beaker, FileText, ExternalLink } from "lucide-react";
import { Loader2 } from "lucide-react";
import { toast } from "sonner";

export default function AnalogDetail() {
  const params = useParams();
  const [, setLocation] = useLocation();
  const analogId = parseInt(params.id || "0");

  const { data: analog, isLoading } = trpc.analog.getById.useQuery({ id: analogId });
  // Test results will be added later
  const testResults: any[] = [];

  if (isLoading) {
    return (
      <div className="container mx-auto py-8 flex items-center justify-center min-h-[400px]">
        <Loader2 className="h-8 w-8 animate-spin text-primary" />
      </div>
    );
  }

  if (!analog) {
    return (
      <div className="container mx-auto py-8">
        <div className="text-center">
          <h2 className="text-2xl font-bold mb-2">Analog Not Found</h2>
          <p className="text-muted-foreground mb-4">The requested analog could not be found.</p>
          <Button onClick={() => setLocation("/admin/dashboard")}>
            <ArrowLeft className="mr-2 h-4 w-4" />
            Back to Dashboard
          </Button>
        </div>
      </div>
    );
  }

  const handleExport = async (format: 'smiles' | 'sdf' | 'pdf') => {
    try {
      let data;
      if (format === 'smiles') {
        data = await trpc.export.smiles.useQuery({ analogId }).data;
      } else if (format === 'sdf') {
        data = await trpc.export.sdf.useQuery({ analogId }).data;
      } else {
        data = await trpc.export.pdf.useQuery({ analogId }).data;
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
        toast.success(`Exported as ${format.toUpperCase()}`);
      }
    } catch (error) {
      toast.error(`Export failed: ${error}`);
    }
  };

  return (
    <div className="container mx-auto py-8">
      {/* Header */}
      <div className="mb-6">
        <Button
          variant="ghost"
          onClick={() => setLocation("/admin/dashboard")}
          className="mb-4"
        >
          <ArrowLeft className="mr-2 h-4 w-4" />
          Back to Dashboard
        </Button>

        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-3xl font-bold mb-2">{analog.compoundName}</h1>
            <p className="text-muted-foreground">
              Parent Compound: <span className="font-medium">{analog.parentCompound}</span>
            </p>
          </div>

          <div className="flex gap-2">
            <Button variant="outline" size="sm" onClick={() => handleExport('smiles')}>
              <Download className="mr-2 h-4 w-4" />
              SMILES
            </Button>
            <Button variant="outline" size="sm" onClick={() => handleExport('sdf')}>
              <Beaker className="mr-2 h-4 w-4" />
              SDF
            </Button>
            <Button variant="outline" size="sm" onClick={() => handleExport('pdf')}>
              <FileText className="mr-2 h-4 w-4" />
              PDF Report
            </Button>
          </div>
        </div>
      </div>

      <div className="grid gap-6 md:grid-cols-[2fr_1fr]">
        {/* Main Content */}
        <div className="space-y-6">
          {/* Chemical Structure */}
          <Card>
            <CardHeader>
              <CardTitle>Chemical Structure</CardTitle>
              <CardDescription>SMILES Notation</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="bg-slate-900 text-slate-100 p-4 rounded font-mono text-sm overflow-x-auto">
                {analog.smiles}
              </div>
            </CardContent>
          </Card>

          {/* Therapeutic Information */}
          <Card>
            <CardHeader>
              <CardTitle>Therapeutic Potential</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              {analog.therapeuticPotential && (
                <div>
                  <h4 className="font-semibold mb-2">Mechanism of Action</h4>
                  <p className="text-sm text-muted-foreground">{analog.therapeuticPotential}</p>
                </div>
              )}

              {analog.keyDifferences && (
                <div>
                  <h4 className="font-semibold mb-2">Key Differences from Parent</h4>
                  <p className="text-sm text-muted-foreground">{analog.keyDifferences}</p>
                </div>
              )}
            </CardContent>
          </Card>

          {/* Test Results */}
          <Card>
            <CardHeader>
              <CardTitle>Cheminformatics Test Results</CardTitle>
              <CardDescription>
                {testResults?.length || 0} test(s) completed
              </CardDescription>
            </CardHeader>
            <CardContent>
              {testResults && testResults.length > 0 ? (
                <Tabs defaultValue={testResults[0]?.testType || "admet"}>
                  <TabsList className="grid w-full grid-cols-4">
                    <TabsTrigger value="admet">ADMET</TabsTrigger>
                    <TabsTrigger value="docking">Docking</TabsTrigger>
                    <TabsTrigger value="toxicity">Toxicity</TabsTrigger>
                    <TabsTrigger value="pkpd">PK/PD</TabsTrigger>
                  </TabsList>

                  {testResults.map((result: any) => (
                    <TabsContent key={result.id} value={result.testType} className="space-y-4">
                      <div className="grid gap-4">
                        <div className="flex items-center justify-between">
                          <Badge variant={result.testStatus === 'completed' ? 'default' : 'secondary'}>
                            {result.testStatus}
                          </Badge>
                          <span className="text-sm text-muted-foreground">
                            {new Date(result.createdAt).toLocaleDateString()}
                          </span>
                        </div>

                        <div className="bg-muted p-4 rounded">
                          <pre className="text-sm overflow-x-auto">
                            {JSON.stringify(JSON.parse(result.results), null, 2)}
                          </pre>
                        </div>
                      </div>
                    </TabsContent>
                  ))}
                </Tabs>
              ) : (
                <div className="text-center py-8 text-muted-foreground">
                  <p>No test results available yet.</p>
                  <Button
                    variant="outline"
                    className="mt-4"
                    onClick={() => setLocation("/testing")}
                  >
                    Run Tests
                  </Button>
                </div>
              )}
            </CardContent>
          </Card>
        </div>

        {/* Sidebar */}
        <div className="space-y-6">
          {/* Scores Card */}
          <Card>
            <CardHeader>
              <CardTitle>Scores & Metrics</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <div className="flex items-center justify-between mb-1">
                  <span className="text-sm font-medium">Confidence</span>
                  <Badge
                    variant={analog.confidenceScore >= 85 ? 'default' : 'secondary'}
                  >
                    {analog.confidenceScore}%
                  </Badge>
                </div>
                <div className="w-full bg-muted rounded-full h-2">
                  <div
                    className="bg-primary h-2 rounded-full"
                    style={{ width: `${analog.confidenceScore}%` }}
                  />
                </div>
              </div>

              <div>
                <div className="flex items-center justify-between mb-1">
                  <span className="text-sm font-medium">Similarity</span>
                  <span className="text-sm">{analog.similarityScore}%</span>
                </div>
                <div className="w-full bg-muted rounded-full h-2">
                  <div
                    className="bg-blue-500 h-2 rounded-full"
                    style={{ width: `${analog.similarityScore}%` }}
                  />
                </div>
              </div>

              <div>
                <div className="flex items-center justify-between mb-1">
                  <span className="text-sm font-medium">Safety</span>
                  <span className="text-sm">{analog.safetyScore}/100</span>
                </div>
                <div className="w-full bg-muted rounded-full h-2">
                  <div
                    className="bg-green-500 h-2 rounded-full"
                    style={{ width: `${analog.safetyScore}%` }}
                  />
                </div>
              </div>

              <div>
                <div className="flex items-center justify-between mb-1">
                  <span className="text-sm font-medium">Efficacy</span>
                  <span className="text-sm">{analog.efficacyScore}/100</span>
                </div>
                <div className="w-full bg-muted rounded-full h-2">
                  <div
                    className="bg-purple-500 h-2 rounded-full"
                    style={{ width: `${analog.efficacyScore}%` }}
                  />
                </div>
              </div>

              <div>
                <div className="flex items-center justify-between mb-1">
                  <span className="text-sm font-medium">Drug-Likeness</span>
                  <span className="text-sm">{analog.drugLikenessScore}/100</span>
                </div>
                <div className="w-full bg-muted rounded-full h-2">
                  <div
                    className="bg-orange-500 h-2 rounded-full"
                    style={{ width: `${analog.drugLikenessScore}%` }}
                  />
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Patent Status */}
          <Card>
            <CardHeader>
              <CardTitle>Patent Status</CardTitle>
            </CardHeader>
            <CardContent>
              <Badge
                variant={analog.patentStatus === 'patent_free' ? 'default' : 'secondary'}
                className="w-full justify-center py-2"
              >
                {analog.patentStatus === 'patent_free'
                  ? 'Patent-Free'
                  : analog.patentStatus === 'patent_opportunity'
                    ? 'Patent Opportunity'
                    : 'Patented'}
              </Badge>

              {analog.patentStatus === 'patent_free' && (
                <p className="text-sm text-muted-foreground mt-3">
                  This compound appears to be free from existing patents and may be eligible
                  for new patent filing.
                </p>
              )}
            </CardContent>
          </Card>

          {/* Market Value */}
          {analog.marketValue && (
            <Card>
              <CardHeader>
                <CardTitle>Market Value</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="text-3xl font-bold text-primary">
                  {analog.marketValue}
                </div>
                <p className="text-sm text-muted-foreground mt-2">
                  Estimated market potential based on therapeutic area and novelty
                </p>
              </CardContent>
            </Card>
          )}

          {/* External Resources */}
          <Card>
            <CardHeader>
              <CardTitle>External Resources</CardTitle>
            </CardHeader>
            <CardContent className="space-y-2">
              <Button variant="outline" className="w-full justify-start" asChild>
                <a
                  href={`https://pubchem.ncbi.nlm.nih.gov/#query=${encodeURIComponent(analog.smiles)}`}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <ExternalLink className="mr-2 h-4 w-4" />
                  Search PubChem
                </a>
              </Button>
              <Button variant="outline" className="w-full justify-start" asChild>
                <a
                  href={`https://www.ebi.ac.uk/chembl/`}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <ExternalLink className="mr-2 h-4 w-4" />
                  Search ChEMBL
                </a>
              </Button>
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}
