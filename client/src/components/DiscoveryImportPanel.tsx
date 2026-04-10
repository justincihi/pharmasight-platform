import React, { useState, useEffect } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { AlertCircle, CheckCircle2, Download, Plus, Trash2, Loader2 } from 'lucide-react';
import { toast } from 'sonner';

interface Discovery {
  id: number;
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  smiles: string;
  confidenceScore: number;
  patentStatus: string;
  discoveredAt: Date;
  discoveredBy: string;
  isPartial: boolean;
  dataCompleteness: number;
}

interface ValidationResult {
  isValid: boolean;
  isPartial: boolean;
  issues: string[];
  missingFields: string[];
  dataCompleteness: number;
  pharmacophoreMatch: number;
  recommendations: string[];
}

export function DiscoveryImportPanel() {
  const [discoveries, setDiscoveries] = useState<Discovery[]>([]);
  const [selectedDiscoveries, setSelectedDiscoveries] = useState<number[]>([]);
  const [validationResults, setValidationResults] = useState<Record<number, ValidationResult>>({});
  const [loading, setLoading] = useState(false);
  const [activeTab, setActiveTab] = useState('all');

  useEffect(() => {
    loadDiscoveries();
  }, []);

  const loadDiscoveries = async () => {
    setLoading(true);
    try {
      // In a real implementation, this would call a tRPC mutation
      // For now, we'll use mock data
      const mockDiscoveries: Discovery[] = [
        {
          id: 1,
          compoundId: 'DISC-001',
          compoundName: 'Ketamine Analog A',
          parentCompound: 'Ketamine',
          smiles: 'CN1C[C@H]2C[C@H]1[C@@H](C(=O)c3ccccc3)C2',
          confidenceScore: 85,
          patentStatus: 'patent-free',
          discoveredAt: new Date(Date.now() - 30 * 24 * 60 * 60 * 1000),
          discoveredBy: 'autonomous-engine',
          isPartial: false,
          dataCompleteness: 95,
        },
        {
          id: 2,
          compoundId: 'DISC-002',
          compoundName: 'Partial Scaffold B',
          parentCompound: 'Ketamine',
          smiles: 'C1=CC=C(C=C1)*',
          confidenceScore: 45,
          patentStatus: 'unknown',
          discoveredAt: new Date(Date.now() - 15 * 24 * 60 * 60 * 1000),
          discoveredBy: 'automated-research',
          isPartial: true,
          dataCompleteness: 35,
        },
      ];
      setDiscoveries(mockDiscoveries);

      // Validate each discovery
      const results: Record<number, ValidationResult> = {};
      mockDiscoveries.forEach((d) => {
        results[d.id] = {
          isValid: d.dataCompleteness >= 50 && !d.isPartial,
          isPartial: d.isPartial,
          issues: d.isPartial ? ['Fragment detected', 'Incomplete SMILES'] : [],
          missingFields: d.dataCompleteness < 100 ? ['mechanism_of_action'] : [],
          dataCompleteness: d.dataCompleteness,
          pharmacophoreMatch: 78,
          recommendations: d.isPartial
            ? ['Auto-query via chatbot for more data', 'Save to Scaffold Library']
            : ['Ready to import to Master List'],
        };
      });
      setValidationResults(results);
    } catch (error) {
      toast.error('Failed to load discoveries');
    } finally {
      setLoading(false);
    }
  };

  const handleSelectDiscovery = (id: number) => {
    setSelectedDiscoveries((prev) =>
      prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]
    );
  };

  const handleImportToMasterList = async () => {
    const validDiscoveries = selectedDiscoveries.filter(
      (id) => validationResults[id]?.isValid
    );

    if (validDiscoveries.length === 0) {
      toast.error('No valid discoveries selected');
      return;
    }

    setLoading(true);
    try {
      // In a real implementation, this would call a tRPC mutation
      toast.success(`${validDiscoveries.length} discoveries imported to Master List`);
      setSelectedDiscoveries([]);
    } catch (error) {
      toast.error('Failed to import discoveries');
    } finally {
      setLoading(false);
    }
  };

  const handleSaveToScaffoldLibrary = async () => {
    const partialDiscoveries = selectedDiscoveries.filter(
      (id) => validationResults[id]?.isPartial
    );

    if (partialDiscoveries.length === 0) {
      toast.error('No partial discoveries selected');
      return;
    }

    setLoading(true);
    try {
      // In a real implementation, this would call a tRPC mutation
      toast.success(`${partialDiscoveries.length} scaffolds saved to Scaffold Library`);
      setSelectedDiscoveries([]);
    } catch (error) {
      toast.error('Failed to save scaffolds');
    } finally {
      setLoading(false);
    }
  };

  const filteredDiscoveries = discoveries.filter((d) => {
    if (activeTab === 'complete') return !d.isPartial;
    if (activeTab === 'partial') return d.isPartial;
    return true;
  });

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>Discovery Audit & Import</CardTitle>
          <CardDescription>
            Review and import discoveries from the last 60 days
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <TabsList className="grid w-full grid-cols-3">
              <TabsTrigger value="all">All ({discoveries.length})</TabsTrigger>
              <TabsTrigger value="complete">
                Complete ({discoveries.filter((d) => !d.isPartial).length})
              </TabsTrigger>
              <TabsTrigger value="partial">
                Partial ({discoveries.filter((d) => d.isPartial).length})
              </TabsTrigger>
            </TabsList>

            <TabsContent value={activeTab} className="space-y-4 mt-4">
              {loading ? (
                <div className="flex items-center justify-center py-8">
                  <Loader2 className="h-6 w-6 animate-spin" />
                </div>
              ) : filteredDiscoveries.length === 0 ? (
                <div className="text-center py-8 text-muted-foreground">
                  No discoveries found
                </div>
              ) : (
                <>
                  {filteredDiscoveries.map((discovery) => {
                    const validation = validationResults[discovery.id];
                    const isSelected = selectedDiscoveries.includes(discovery.id);

                    return (
                      <Card
                        key={discovery.id}
                        className={`cursor-pointer transition-all ${
                          isSelected ? 'ring-2 ring-blue-500' : ''
                        }`}
                        onClick={() => handleSelectDiscovery(discovery.id)}
                      >
                        <CardContent className="pt-6">
                          <div className="space-y-4">
                            {/* Header */}
                            <div className="flex items-start justify-between">
                              <div className="flex-1">
                                <h4 className="font-semibold">{discovery.compoundName}</h4>
                                <p className="text-sm text-muted-foreground">
                                  Parent: {discovery.parentCompound}
                                </p>
                              </div>
                              <div className="flex gap-2">
                                <Badge
                                  variant={
                                    discovery.isPartial ? 'destructive' : 'default'
                                  }
                                >
                                  {discovery.isPartial ? 'Partial' : 'Complete'}
                                </Badge>
                                <Badge variant="outline">
                                  {discovery.patentStatus}
                                </Badge>
                              </div>
                            </div>

                            {/* SMILES */}
                            <div className="bg-muted p-3 rounded-lg">
                              <div className="text-xs text-muted-foreground mb-1">
                                SMILES
                              </div>
                              <code className="text-xs break-all">
                                {discovery.smiles}
                              </code>
                            </div>

                            {/* Metrics */}
                            <div className="grid grid-cols-3 gap-3">
                              <div className="p-2 border rounded">
                                <div className="text-xs text-muted-foreground">
                                  Confidence
                                </div>
                                <div className="font-semibold">
                                  {discovery.confidenceScore}%
                                </div>
                              </div>
                              <div className="p-2 border rounded">
                                <div className="text-xs text-muted-foreground">
                                  Data Completeness
                                </div>
                                <div className="font-semibold">
                                  {validation?.dataCompleteness || 0}%
                                </div>
                              </div>
                              <div className="p-2 border rounded">
                                <div className="text-xs text-muted-foreground">
                                  Pharmacophore Match
                                </div>
                                <div className="font-semibold">
                                  {validation?.pharmacophoreMatch || 0}%
                                </div>
                              </div>
                            </div>

                            {/* Validation Issues */}
                            {validation?.issues && validation.issues.length > 0 && (
                              <div className="bg-red-50 border border-red-200 rounded-lg p-3">
                                <div className="flex gap-2">
                                  <AlertCircle className="h-4 w-4 text-red-600 flex-shrink-0 mt-0.5" />
                                  <div className="text-sm text-red-800">
                                    {validation.issues.join(', ')}
                                  </div>
                                </div>
                              </div>
                            )}

                            {/* Recommendations */}
                            {validation?.recommendations && (
                              <div className="bg-blue-50 border border-blue-200 rounded-lg p-3">
                                <div className="text-sm font-medium text-blue-900 mb-2">
                                  Recommendations:
                                </div>
                                <ul className="text-sm text-blue-800 space-y-1">
                                  {validation.recommendations.map((rec, idx) => (
                                    <li key={idx}>• {rec}</li>
                                  ))}
                                </ul>
                              </div>
                            )}

                            {/* Validation Status */}
                            <div className="flex items-center gap-2">
                              {validation?.isValid ? (
                                <>
                                  <CheckCircle2 className="h-4 w-4 text-green-600" />
                                  <span className="text-sm text-green-600">
                                    Ready for Master List
                                  </span>
                                </>
                              ) : validation?.isPartial ? (
                                <>
                                  <AlertCircle className="h-4 w-4 text-yellow-600" />
                                  <span className="text-sm text-yellow-600">
                                    Suitable for Scaffold Library
                                  </span>
                                </>
                              ) : (
                                <>
                                  <AlertCircle className="h-4 w-4 text-red-600" />
                                  <span className="text-sm text-red-600">
                                    Needs more data
                                  </span>
                                </>
                              )}
                            </div>
                          </div>
                        </CardContent>
                      </Card>
                    );
                  })}

                  {/* Action Buttons */}
                  {selectedDiscoveries.length > 0 && (
                    <div className="flex gap-2 sticky bottom-0 bg-background pt-4 border-t">
                      <Button
                        onClick={handleImportToMasterList}
                        disabled={
                          selectedDiscoveries.filter(
                            (id) => validationResults[id]?.isValid
                          ).length === 0
                        }
                      >
                        <Plus className="mr-2 h-4 w-4" />
                        Import to Master List
                      </Button>
                      <Button
                        variant="outline"
                        onClick={handleSaveToScaffoldLibrary}
                        disabled={
                          selectedDiscoveries.filter(
                            (id) => validationResults[id]?.isPartial
                          ).length === 0
                        }
                      >
                        <Download className="mr-2 h-4 w-4" />
                        Save to Scaffold Library
                      </Button>
                      <Button
                        variant="ghost"
                        onClick={() => setSelectedDiscoveries([])}
                      >
                        <Trash2 className="mr-2 h-4 w-4" />
                        Clear Selection
                      </Button>
                    </div>
                  )}
                </>
              )}
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>
    </div>
  );
}
