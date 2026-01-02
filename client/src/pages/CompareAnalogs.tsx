import { useState } from 'react';
import DashboardLayout from '@/components/DashboardLayout';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Loader2, Search, X, BarChart3 } from 'lucide-react';
import { trpc } from '@/lib/trpc';
import { RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar, Legend, ResponsiveContainer } from 'recharts';

export default function CompareAnalogs() {
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedAnalogs, setSelectedAnalogs] = useState<any[]>([]);

  const { data: searchResults, isLoading } = trpc.analog.search.useQuery(
    { query: searchQuery },
    { enabled: searchQuery.length > 2 }
  );

  const handleAddAnalog = (analog: any) => {
    if (selectedAnalogs.length >= 5) {
      alert('⚠️ Maximum 5 analogs can be compared at once');
      return;
    }
    if (selectedAnalogs.find(a => a.id === analog.id)) {
      alert('⚠️ This analog is already selected');
      return;
    }
    setSelectedAnalogs([...selectedAnalogs, analog]);
    setSearchQuery('');
  };

  const handleRemoveAnalog = (analogId: number) => {
    setSelectedAnalogs(selectedAnalogs.filter(a => a.id !== analogId));
  };

  // Prepare radar chart data
  const radarData = [
    {
      property: 'Confidence',
      ...Object.fromEntries(
        selectedAnalogs.map(a => [a.compoundName, a.confidenceScore])
      ),
    },
    {
      property: 'Safety',
      ...Object.fromEntries(
        selectedAnalogs.map(a => [a.compoundName, a.safetyScore])
      ),
    },
    {
      property: 'Efficacy',
      ...Object.fromEntries(
        selectedAnalogs.map(a => [a.compoundName, a.efficacyScore])
      ),
    },
    {
      property: 'Drug Likeness',
      ...Object.fromEntries(
        selectedAnalogs.map(a => [a.compoundName, a.drugLikenessScore])
      ),
    },
    {
      property: 'Similarity',
      ...Object.fromEntries(
        selectedAnalogs.map(a => [a.compoundName, a.similarityScore])
      ),
    },
  ];

  const colors = ['#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6'];

  return (
    <DashboardLayout>
      <div className="space-y-6">
        {/* Header */}
        <div>
          <h1 className="text-3xl font-bold text-gray-900">Compare Analogs</h1>
          <p className="text-gray-600 mt-1">
            Side-by-side comparison of molecular properties and ADMET profiles
          </p>
        </div>

        {/* Search and Selection */}
        <Card>
          <CardHeader>
            <CardTitle>Select Analogs to Compare (Max 5)</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="relative">
              <Search className="absolute left-3 top-3 w-4 h-4 text-gray-400" />
              <Input
                placeholder="Search by compound name or SMILES..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                className="pl-10"
              />
            </div>

            {/* Search Results */}
            {searchQuery.length > 2 && (
              <div className="border border-gray-200 rounded-lg max-h-60 overflow-y-auto">
                {isLoading ? (
                  <div className="p-4 text-center">
                    <Loader2 className="w-6 h-6 animate-spin mx-auto text-blue-600" />
                  </div>
                ) : searchResults && searchResults.length > 0 ? (
                  <div className="divide-y">
                    {searchResults.map((analog: any) => (
                      <div
                        key={analog.id}
                        className="p-3 hover:bg-gray-50 cursor-pointer flex justify-between items-center"
                        onClick={() => handleAddAnalog(analog)}
                      >
                        <div>
                          <p className="font-semibold text-gray-900">{analog.compoundName}</p>
                          <p className="text-sm text-gray-600">{analog.parentCompound}</p>
                        </div>
                        <Button size="sm" variant="outline">
                          Add
                        </Button>
                      </div>
                    ))}
                  </div>
                ) : (
                  <p className="p-4 text-center text-gray-500">No results found</p>
                )}
              </div>
            )}

            {/* Selected Analogs */}
            {selectedAnalogs.length > 0 && (
              <div className="space-y-2">
                <p className="text-sm font-medium text-gray-700">
                  Selected ({selectedAnalogs.length}/5):
                </p>
                <div className="flex flex-wrap gap-2">
                  {selectedAnalogs.map((analog, idx) => (
                    <div
                      key={analog.id}
                      className="flex items-center gap-2 px-3 py-1.5 rounded-full text-sm"
                      style={{ backgroundColor: colors[idx] + '20', color: colors[idx] }}
                    >
                      <span className="font-medium">{analog.compoundName}</span>
                      <button
                        onClick={() => handleRemoveAnalog(analog.id)}
                        className="hover:opacity-70"
                      >
                        <X className="w-4 h-4" />
                      </button>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Comparison Results */}
        {selectedAnalogs.length >= 2 && (
          <>
            {/* Radar Chart */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <BarChart3 className="w-5 h-5" />
                  Property Comparison
                </CardTitle>
              </CardHeader>
              <CardContent>
                <ResponsiveContainer width="100%" height={400}>
                  <RadarChart data={radarData}>
                    <PolarGrid />
                    <PolarAngleAxis dataKey="property" />
                    <PolarRadiusAxis angle={90} domain={[0, 100]} />
                    {selectedAnalogs.map((analog, idx) => (
                      <Radar
                        key={analog.id}
                        name={analog.compoundName}
                        dataKey={analog.compoundName}
                        stroke={colors[idx]}
                        fill={colors[idx]}
                        fillOpacity={0.3}
                      />
                    ))}
                    <Legend />
                  </RadarChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>

            {/* Side-by-Side Table */}
            <Card>
              <CardHeader>
                <CardTitle>Detailed Comparison</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b">
                        <th className="text-left p-3 font-semibold text-gray-700">Property</th>
                        {selectedAnalogs.map((analog, idx) => (
                          <th
                            key={analog.id}
                            className="text-left p-3 font-semibold"
                            style={{ color: colors[idx] }}
                          >
                            {analog.compoundName}
                          </th>
                        ))}
                      </tr>
                    </thead>
                    <tbody className="divide-y">
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Parent Compound</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.parentCompound}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">SMILES</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3 text-xs font-mono">{a.smiles}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Confidence Score</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.confidenceScore}%</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Safety Score</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.safetyScore}%</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Efficacy Score</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.efficacyScore}%</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Drug Likeness</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.drugLikenessScore}%</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Similarity Score</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.similarityScore}%</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Patent Status</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3 capitalize">{a.patentStatus}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Molecular Weight</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.molecularWeight || 'N/A'}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">LogP</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.logP || 'N/A'}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">H-Bond Donors</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.hBondDonors ?? 'N/A'}</td>
                        ))}
                      </tr>
                      <tr>
                        <td className="p-3 font-medium text-gray-700">H-Bond Acceptors</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">{a.hBondAcceptors ?? 'N/A'}</td>
                        ))}
                      </tr>
                      {selectedAnalogs.some(a => a.bindingAffinity) && (
                        <tr>
                          <td className="p-3 font-medium text-gray-700">Binding Affinity</td>
                          {selectedAnalogs.map(a => (
                            <td key={a.id} className="p-3">
                              {a.bindingAffinity ? `${a.bindingAffinity} kcal/mol` : 'Not tested'}
                            </td>
                          ))}
                        </tr>
                      )}
                      <tr>
                        <td className="p-3 font-medium text-gray-700">Discovered</td>
                        {selectedAnalogs.map(a => (
                          <td key={a.id} className="p-3">
                            {new Date(a.discoveredAt).toLocaleDateString()}
                          </td>
                        ))}
                      </tr>
                    </tbody>
                  </table>
                </div>
              </CardContent>
            </Card>

            {/* Export Button */}
            <div className="flex justify-end">
              <Button onClick={() => alert('Export functionality coming soon!')}>
                Export Comparison Report
              </Button>
            </div>
          </>
        )}

        {/* Empty State */}
        {selectedAnalogs.length < 2 && (
          <Card>
            <CardContent className="py-12 text-center">
              <BarChart3 className="w-16 h-16 mx-auto text-gray-300 mb-4" />
              <p className="text-gray-600 text-lg">
                Select at least 2 analogs to compare
              </p>
              <p className="text-gray-500 text-sm mt-2">
                Use the search box above to find and add compounds
              </p>
            </CardContent>
          </Card>
        )}
      </div>
    </DashboardLayout>
  );
}
