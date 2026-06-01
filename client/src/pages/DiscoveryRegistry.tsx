import { useState, useMemo } from 'react';
import { trpc } from '@/lib/trpc';
import DashboardLayout from '@/components/DashboardLayout';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Loader2, Search, Award, Download, Filter, ExternalLink, Copy, CheckCircle2, TrendingUp, Shield, AlertTriangle } from 'lucide-react';
import { toast } from 'sonner';

const CONFIDENCE_TIER = (score: number) => {
  if (score >= 90) return { label: 'Tier 1 — High Confidence', color: 'text-green-700 bg-green-50 border-green-300', dot: 'bg-green-500' };
  if (score >= 70) return { label: 'Tier 2 — Moderate Confidence', color: 'text-yellow-700 bg-yellow-50 border-yellow-300', dot: 'bg-yellow-500' };
  if (score >= 50) return { label: 'Tier 3 — Low Confidence', color: 'text-orange-700 bg-orange-50 border-orange-300', dot: 'bg-orange-500' };
  return { label: 'Tier 4 — Speculative', color: 'text-gray-600 bg-gray-50 border-gray-300', dot: 'bg-gray-400' };
};

const PATENT_CONFIG: Record<string, { color: string; icon: typeof Shield }> = {
  'patent-free': { color: 'text-green-700 bg-green-50 border-green-200', icon: CheckCircle2 },
  'patent-opportunity': { color: 'text-blue-700 bg-blue-50 border-blue-200', icon: TrendingUp },
  patented: { color: 'text-red-700 bg-red-50 border-red-200', icon: Shield },
  unknown: { color: 'text-gray-600 bg-gray-50 border-gray-200', icon: AlertTriangle },
};

export default function DiscoveryRegistry() {
  const [search, setSearch] = useState('');
  const [minConfidence, setMinConfidence] = useState<number>(0);
  const [patentFilter, setPatentFilter] = useState('all');
  const [tierFilter, setTierFilter] = useState('all');
  const [sortBy, setSortBy] = useState<'confidence' | 'date' | 'name'>('confidence');
  const [copied, setCopied] = useState<number | null>(null);

  const registryQuery = trpc.encyclopedia.search.useQuery({
    query: search || undefined,
    patentStatus: patentFilter !== 'all' ? patentFilter : undefined,
    minConfidence: minConfidence > 0 ? minConfidence : undefined,
    limit: 500,
  });

  const compounds = useMemo(() => {
    const data = registryQuery.data ?? [];
    let filtered = data;

    if (tierFilter !== 'all') {
      const tierMin = tierFilter === 'tier1' ? 90 : tierFilter === 'tier2' ? 70 : tierFilter === 'tier3' ? 50 : 0;
      const tierMax = tierFilter === 'tier1' ? 100 : tierFilter === 'tier2' ? 89 : tierFilter === 'tier3' ? 69 : 49;
      filtered = filtered.filter(c => c.confidenceScore >= tierMin && c.confidenceScore <= tierMax);
    }

    return [...filtered].sort((a, b) => {
      if (sortBy === 'confidence') return b.confidenceScore - a.confidenceScore;
      if (sortBy === 'name') return a.compoundName.localeCompare(b.compoundName);
      return new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime();
    });
  }, [registryQuery.data, tierFilter, sortBy]);

  // Stats
  const stats = useMemo(() => {
    const all = registryQuery.data ?? [];
    return {
      total: all.length,
      tier1: all.filter(c => c.confidenceScore >= 90).length,
      tier2: all.filter(c => c.confidenceScore >= 70 && c.confidenceScore < 90).length,
      patentFree: all.filter(c => c.patentStatus === 'patent-free').length,
      patentOpp: all.filter(c => c.patentStatus === 'patent-opportunity').length,
    };
  }, [registryQuery.data]);

  const handleCopySmiles = (id: number, smiles: string) => {
    navigator.clipboard.writeText(smiles);
    setCopied(id);
    setTimeout(() => setCopied(null), 2000);
    toast.success('SMILES copied');
  };

  const handleExportCSV = () => {
    const headers = ['ID', 'Compound Name', 'SMILES', 'Confidence', 'Tier', 'Patent Status', 'Safety Score', 'Efficacy Score', 'Drug-likeness', 'Parent Compound', 'Generation', 'Discovery Date', 'PubChem CID', 'ChEMBL ID'];
    const rows = compounds.map(c => [
      c.compoundId,
      c.compoundName,
      c.smiles,
      c.confidenceScore,
      CONFIDENCE_TIER(c.confidenceScore).label,
      c.patentStatus,
      c.safetyScore,
      c.efficacyScore,
      c.drugLikenessScore,
      c.parentCompound,
      c.optimizationGeneration,
      new Date(c.createdAt).toISOString().split('T')[0],
      c.pubchemCid ?? '',
      c.chemblId ?? '',
    ]);
    const csv = [headers, ...rows].map(r => r.map(v => `"${String(v ?? '').replace(/"/g, '""')}"`).join(',')).join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `pharmasight-discovery-registry-${new Date().toISOString().split('T')[0]}.csv`;
    a.click();
    URL.revokeObjectURL(url);
    toast.success(`Exported ${compounds.length} compounds`);
  };

  return (
    <DashboardLayout>
      <div className="p-6 space-y-6">
        {/* Header */}
        <div className="flex items-start justify-between gap-4">
          <div className="flex items-center gap-3">
            <Award className="h-7 w-7 text-amber-500" />
            <div>
              <h1 className="text-2xl font-bold">Discovery Registry</h1>
              <p className="text-sm text-muted-foreground">
                Confidence-scored analog registry with IP status, audit trail, and filterable tiers
              </p>
            </div>
          </div>
          <Button onClick={handleExportCSV} variant="outline" size="sm" disabled={compounds.length === 0}>
            <Download className="h-4 w-4 mr-2" /> Export CSV
          </Button>
        </div>

        {/* Stats row */}
        <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
          {[
            { label: 'Total Compounds', value: stats.total, color: 'text-foreground' },
            { label: 'Tier 1 (≥90%)', value: stats.tier1, color: 'text-green-600' },
            { label: 'Tier 2 (70–89%)', value: stats.tier2, color: 'text-yellow-600' },
            { label: 'Patent-Free', value: stats.patentFree, color: 'text-blue-600' },
            { label: 'IP Opportunity', value: stats.patentOpp, color: 'text-purple-600' },
          ].map(({ label, value, color }) => (
            <Card key={label}>
              <CardContent className="pt-4 text-center">
                <p className={`text-2xl font-bold ${color}`}>{registryQuery.isLoading ? '–' : value}</p>
                <p className="text-xs text-muted-foreground mt-0.5">{label}</p>
              </CardContent>
            </Card>
          ))}
        </div>

        {/* Filters */}
        <Card>
          <CardContent className="pt-4">
            <div className="flex flex-wrap items-center gap-3">
              <div className="relative flex-1 min-w-48">
                <Search className="absolute left-3 top-2.5 h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="Search compounds, SMILES, parent..."
                  value={search}
                  onChange={(e) => setSearch(e.target.value)}
                  className="pl-9"
                />
              </div>
              <Select value={tierFilter} onValueChange={setTierFilter}>
                <SelectTrigger className="w-44">
                  <Filter className="h-3.5 w-3.5 mr-1.5 text-muted-foreground" />
                  <SelectValue placeholder="Confidence tier" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">All tiers</SelectItem>
                  <SelectItem value="tier1">Tier 1 — ≥ 90%</SelectItem>
                  <SelectItem value="tier2">Tier 2 — 70–89%</SelectItem>
                  <SelectItem value="tier3">Tier 3 — 50–69%</SelectItem>
                  <SelectItem value="tier4">Tier 4 — &lt; 50%</SelectItem>
                </SelectContent>
              </Select>
              <Select value={patentFilter} onValueChange={setPatentFilter}>
                <SelectTrigger className="w-40">
                  <SelectValue placeholder="IP status" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="all">All IP status</SelectItem>
                  <SelectItem value="patent-free">Patent-free</SelectItem>
                  <SelectItem value="patent-opportunity">IP Opportunity</SelectItem>
                  <SelectItem value="patented">Patented</SelectItem>
                  <SelectItem value="unknown">Unknown</SelectItem>
                </SelectContent>
              </Select>
              <Select value={sortBy} onValueChange={(v) => setSortBy(v as typeof sortBy)}>
                <SelectTrigger className="w-36">
                  <SelectValue placeholder="Sort by" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="confidence">By Confidence</SelectItem>
                  <SelectItem value="date">By Date</SelectItem>
                  <SelectItem value="name">By Name</SelectItem>
                </SelectContent>
              </Select>
              <p className="text-xs text-muted-foreground whitespace-nowrap">
                {registryQuery.isLoading ? 'Loading...' : `${compounds.length} result${compounds.length !== 1 ? 's' : ''}`}
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Registry table */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium">Compound Registry</CardTitle>
            <CardDescription className="text-xs">
              All discovered analogs, sorted by confidence score. Click a row to copy SMILES.
            </CardDescription>
          </CardHeader>
          <CardContent>
            {registryQuery.isLoading && (
              <div className="flex items-center justify-center py-12">
                <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
              </div>
            )}
            {!registryQuery.isLoading && compounds.length === 0 && (
              <div className="text-center py-12">
                <Award className="h-10 w-10 text-muted-foreground mx-auto mb-2 opacity-30" />
                <p className="text-muted-foreground text-sm">No compounds match the current filters.</p>
                <p className="text-xs text-muted-foreground mt-1">Generate analogs to populate the registry.</p>
              </div>
            )}
            {compounds.length > 0 && (
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b text-xs text-muted-foreground">
                      <th className="text-left py-2 pr-4 font-medium">Compound</th>
                      <th className="text-left py-2 pr-4 font-medium">Confidence</th>
                      <th className="text-left py-2 pr-4 font-medium">IP Status</th>
                      <th className="text-left py-2 pr-4 font-medium">Scores</th>
                      <th className="text-left py-2 pr-4 font-medium">Parent</th>
                      <th className="text-left py-2 pr-4 font-medium">Registered</th>
                      <th className="text-left py-2 font-medium">Links</th>
                    </tr>
                  </thead>
                  <tbody>
                    {compounds.map((c) => {
                      const tier = CONFIDENCE_TIER(c.confidenceScore);
                      const patent = PATENT_CONFIG[c.patentStatus] ?? PATENT_CONFIG.unknown;
                      const PatentIcon = patent.icon;
                      return (
                        <tr key={c.id} className="border-b last:border-0 hover:bg-muted/30 transition-colors">
                          <td className="py-2.5 pr-4">
                            <div className="flex items-start gap-2">
                              <span className={`w-2 h-2 rounded-full mt-1.5 shrink-0 ${tier.dot}`} />
                              <div className="min-w-0">
                                <p className="font-medium truncate max-w-[160px]">{c.compoundName}</p>
                                <div className="flex items-center gap-1 mt-0.5">
                                  <code className="text-xs text-muted-foreground font-mono truncate max-w-[140px]">
                                    {c.smiles?.slice(0, 30)}{(c.smiles?.length ?? 0) > 30 ? '…' : ''}
                                  </code>
                                  <button
                                    onClick={() => handleCopySmiles(c.id, c.smiles)}
                                    className="text-muted-foreground hover:text-foreground shrink-0"
                                  >
                                    {copied === c.id
                                      ? <CheckCircle2 className="h-3 w-3 text-green-500" />
                                      : <Copy className="h-3 w-3" />}
                                  </button>
                                </div>
                              </div>
                            </div>
                          </td>
                          <td className="py-2.5 pr-4">
                            <div>
                              <div className="flex items-center gap-1.5">
                                <div className="w-16 bg-gray-200 dark:bg-gray-700 rounded-full h-1.5">
                                  <div
                                    className={`h-1.5 rounded-full ${c.confidenceScore >= 90 ? 'bg-green-500' : c.confidenceScore >= 70 ? 'bg-yellow-500' : 'bg-orange-500'}`}
                                    style={{ width: `${c.confidenceScore}%` }}
                                  />
                                </div>
                                <span className="font-medium">{c.confidenceScore}%</span>
                              </div>
                              <p className="text-xs text-muted-foreground mt-0.5">{tier.label.split(' — ')[0]}</p>
                            </div>
                          </td>
                          <td className="py-2.5 pr-4">
                            <Badge variant="outline" className={`text-xs gap-1 ${patent.color}`}>
                              <PatentIcon className="h-2.5 w-2.5" />
                              {c.patentStatus}
                            </Badge>
                          </td>
                          <td className="py-2.5 pr-4">
                            <div className="flex gap-1.5 text-xs">
                              <span className="text-muted-foreground">S:{c.safetyScore}</span>
                              <span className="text-muted-foreground">E:{c.efficacyScore}</span>
                              <span className="text-muted-foreground">D:{c.drugLikenessScore}</span>
                            </div>
                          </td>
                          <td className="py-2.5 pr-4">
                            <p className="text-xs text-muted-foreground truncate max-w-[100px]">{c.parentCompound}</p>
                            <p className="text-xs text-muted-foreground">Gen {c.optimizationGeneration}</p>
                          </td>
                          <td className="py-2.5 pr-4">
                            <p className="text-xs text-muted-foreground">
                              {new Date(c.createdAt).toLocaleDateString()}
                            </p>
                          </td>
                          <td className="py-2.5">
                            <div className="flex gap-1">
                              {c.pubchemCid && (
                                <a href={`https://pubchem.ncbi.nlm.nih.gov/compound/${c.pubchemCid}`} target="_blank" rel="noopener noreferrer">
                                  <Badge variant="outline" className="text-xs h-5 px-1 cursor-pointer hover:bg-muted">
                                    <ExternalLink className="h-2.5 w-2.5 mr-0.5" />PC
                                  </Badge>
                                </a>
                              )}
                              {c.chemblId && (
                                <a href={`https://www.ebi.ac.uk/chembl/compound_report_card/${c.chemblId}/`} target="_blank" rel="noopener noreferrer">
                                  <Badge variant="outline" className="text-xs h-5 px-1 cursor-pointer hover:bg-muted">
                                    <ExternalLink className="h-2.5 w-2.5 mr-0.5" />CB
                                  </Badge>
                                </a>
                              )}
                            </div>
                          </td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
            )}
          </CardContent>
        </Card>
      </div>
    </DashboardLayout>
  );
}
