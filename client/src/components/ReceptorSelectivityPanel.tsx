import React, { useMemo } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Target, TrendingUp } from 'lucide-react';

interface DockingResult {
  compoundName: string;
  receptor: string;
  bindingAffinity: number;
}

interface ReceptorSelectivityPanelProps {
  results: DockingResult[];
  primaryReceptor?: string;
}

export const ReceptorSelectivityPanel: React.FC<ReceptorSelectivityPanelProps> = ({
  results,
  primaryReceptor = 'NMDA_GluN2A',
}) => {
  const selectivityData = useMemo(() => {
    if (results.length === 0) return [];

    // Group by compound
    const grouped = results.reduce(
      (acc, result) => {
        if (!acc[result.compoundName]) {
          acc[result.compoundName] = {};
        }
        acc[result.compoundName][result.receptor] = result.bindingAffinity;
        return acc;
      },
      {} as Record<string, Record<string, number>>
    );

    // Calculate selectivity scores
    return Object.entries(grouped)
      .map(([compound, affinities]) => {
        const primaryAffinity = affinities[primaryReceptor] ?? 0;
        const otherAffinities = Object.entries(affinities)
          .filter(([receptor]) => receptor !== primaryReceptor)
          .map(([, affinity]) => affinity);

        const avgOffTarget = otherAffinities.length > 0
          ? otherAffinities.reduce((a, b) => a + b, 0) / otherAffinities.length
          : 0;

        const selectivityScore = primaryAffinity - avgOffTarget;

        return {
          compound,
          primaryAffinity,
          avgOffTarget,
          selectivityScore,
          affinities,
        };
      })
      .sort((a, b) => b.selectivityScore - a.selectivityScore);
  }, [results, primaryReceptor]);

  const receptors = useMemo(() => {
    return Array.from(new Set(results.map((r) => r.receptor)));
  }, [results]);

  const getAffinityColor = (affinity: number) => {
    if (affinity < -8) return 'bg-green-100 text-green-900';
    if (affinity < -6) return 'bg-blue-100 text-blue-900';
    if (affinity < -4) return 'bg-yellow-100 text-yellow-900';
    return 'bg-red-100 text-red-900';
  };

  const getSelectivityColor = (score: number) => {
    if (score > 2) return 'text-green-600';
    if (score > 0) return 'text-blue-600';
    return 'text-red-600';
  };

  return (
    <Card className="w-full">
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Target className="w-5 h-5 text-blue-500" />
          Receptor Selectivity Analysis
        </CardTitle>
        <CardDescription>Analyze binding selectivity across multiple receptors</CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <Tabs defaultValue="selectivity" className="w-full">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="selectivity">Selectivity</TabsTrigger>
            <TabsTrigger value="heatmap">Heatmap</TabsTrigger>
            <TabsTrigger value="analysis">Analysis</TabsTrigger>
          </TabsList>

          <TabsContent value="selectivity" className="space-y-4 mt-4">
            {selectivityData.length === 0 ? (
              <div className="text-center py-8 text-gray-500">
                <p>No data available. Run batch docking first.</p>
              </div>
            ) : (
              <div className="space-y-3">
                {selectivityData.map((item) => (
                  <div key={item.compound} className="p-4 border rounded-lg">
                    <div className="flex justify-between items-start mb-3">
                      <div>
                        <p className="font-semibold">{item.compound}</p>
                        <p className="text-sm text-gray-600">Primary: {primaryReceptor}</p>
                      </div>
                      <div className={`text-lg font-bold ${getSelectivityColor(item.selectivityScore)}`}>
                        {item.selectivityScore > 0 ? '+' : ''}{item.selectivityScore.toFixed(2)}
                      </div>
                    </div>

                    <div className="grid grid-cols-3 gap-2 text-sm">
                      <div className="p-2 bg-gray-50 rounded">
                        <p className="text-gray-600">Primary Affinity</p>
                        <p className="font-mono font-bold">{item.primaryAffinity.toFixed(2)}</p>
                      </div>
                      <div className="p-2 bg-gray-50 rounded">
                        <p className="text-gray-600">Avg Off-Target</p>
                        <p className="font-mono font-bold">{item.avgOffTarget.toFixed(2)}</p>
                      </div>
                      <div className="p-2 bg-blue-50 rounded">
                        <p className="text-gray-600">Selectivity</p>
                        <p className="font-mono font-bold text-blue-600">{item.selectivityScore.toFixed(2)}</p>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            )}
          </TabsContent>

          <TabsContent value="heatmap" className="space-y-4 mt-4">
            {results.length === 0 ? (
              <div className="text-center py-8 text-gray-500">
                <p>No data available.</p>
              </div>
            ) : (
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b">
                      <th className="text-left p-2 font-semibold">Compound</th>
                      {receptors.map((receptor) => (
                        <th key={receptor} className="text-center p-2 font-semibold text-xs">
                          {receptor.split('_')[0]}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {selectivityData.map((item) => (
                      <tr key={item.compound} className="border-b hover:bg-gray-50">
                        <td className="p-2 font-semibold">{item.compound}</td>
                        {receptors.map((receptor) => {
                          const affinity = item.affinities[receptor];
                          return (
                            <td key={`${item.compound}-${receptor}`} className="text-center p-2">
                              <Badge className={getAffinityColor(affinity)}>
                                {affinity?.toFixed(1) ?? 'N/A'}
                              </Badge>
                            </td>
                          );
                        })}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </TabsContent>

          <TabsContent value="analysis" className="space-y-4 mt-4">
            <div className="space-y-3">
              <div className="p-4 bg-blue-50 rounded-lg border border-blue-200">
                <h4 className="font-semibold text-blue-900 mb-2 flex items-center gap-2">
                  <TrendingUp className="w-4 h-4" />
                  Selectivity Interpretation
                </h4>
                <ul className="text-sm text-blue-800 space-y-1">
                  <li>• <strong>Score &gt; 2:</strong> Excellent selectivity (strong primary, weak off-target)</li>
                  <li>• <strong>Score 0-2:</strong> Good selectivity (moderate primary, some off-target)</li>
                  <li>• <strong>Score &lt; 0:</strong> Poor selectivity (weak primary or strong off-target)</li>
                </ul>
              </div>

              <div className="p-4 bg-green-50 rounded-lg border border-green-200">
                <h4 className="font-semibold text-green-900 mb-2">Binding Affinity Guide</h4>
                <div className="grid grid-cols-2 gap-2 text-sm text-green-800">
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 bg-green-400 rounded"></div>
                    <span>&lt; -8 kcal/mol: Excellent</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 bg-blue-400 rounded"></div>
                    <span>-8 to -6: Good</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 bg-yellow-400 rounded"></div>
                    <span>-6 to -4: Moderate</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 bg-red-400 rounded"></div>
                    <span>&gt; -4: Weak</span>
                  </div>
                </div>
              </div>

              {selectivityData.length > 0 && (
                <div className="p-4 bg-purple-50 rounded-lg border border-purple-200">
                  <h4 className="font-semibold text-purple-900 mb-2">Top Compound</h4>
                  <p className="text-sm text-purple-800">
                    <strong>{selectivityData[0].compound}</strong> shows the best selectivity profile with a score of{' '}
                    <strong>{selectivityData[0].selectivityScore.toFixed(2)}</strong>, indicating strong primary receptor
                    binding with minimal off-target effects.
                  </p>
                </div>
              )}
            </div>
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  );
};

export default ReceptorSelectivityPanel;
