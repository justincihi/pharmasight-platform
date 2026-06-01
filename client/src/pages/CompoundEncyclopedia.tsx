import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import DashboardLayout from '@/components/DashboardLayout';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Textarea } from '@/components/ui/textarea';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Loader2, Search, BookOpen, FlaskConical, Dna, Activity, ChevronRight, Edit3, Save, X, ExternalLink, Copy, CheckCircle2, AlertCircle } from 'lucide-react';
import { toast } from 'sonner';

const CONFIDENCE_COLOR = (score: number) =>
  score >= 80 ? 'text-green-600 bg-green-50 border-green-200' :
  score >= 60 ? 'text-yellow-600 bg-yellow-50 border-yellow-200' :
  'text-red-600 bg-red-50 border-red-200';

const PATENT_COLOR: Record<string, string> = {
  'patent-free': 'text-green-700 bg-green-50 border-green-200',
  'patent-opportunity': 'text-blue-700 bg-blue-50 border-blue-200',
  patented: 'text-red-700 bg-red-50 border-red-200',
  unknown: 'text-gray-600 bg-gray-50 border-gray-200',
};

export default function CompoundEncyclopedia() {
  const [query, setQuery] = useState('');
  const [patentFilter, setPatentFilter] = useState<string>('all');
  const [minConfidence, setMinConfidence] = useState<number>(0);
  const [selectedId, setSelectedId] = useState<number | null>(null);
  const [editing, setEditing] = useState(false);
  const [editFields, setEditFields] = useState<Record<string, string>>({});
  const [copied, setCopied] = useState(false);

  const searchQuery = trpc.encyclopedia.search.useQuery({
    query: query || undefined,
    patentStatus: patentFilter !== 'all' ? patentFilter : undefined,
    minConfidence: minConfidence > 0 ? minConfidence : undefined,
    limit: 100,
  });

  const detailQuery = trpc.encyclopedia.getCompound.useQuery(
    { analogId: selectedId! },
    { enabled: selectedId !== null }
  );

  const updateNotesMutation = trpc.encyclopedia.updateNotes.useMutation({
    onSuccess: () => {
      toast.success('Notes saved');
      setEditing(false);
      detailQuery.refetch();
    },
    onError: (e) => toast.error(`Save failed: ${e.message}`),
  });

  const compounds = searchQuery.data ?? [];
  const detail = detailQuery.data;

  const handleCopySmiles = (smiles: string) => {
    navigator.clipboard.writeText(smiles);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
    toast.success('SMILES copied');
  };

  const handleSaveNotes = () => {
    if (!selectedId) return;
    updateNotesMutation.mutate({ analogId: selectedId, ...editFields });
  };

  const startEdit = () => {
    if (!detail?.analog) return;
    setEditFields({
      mechanismOfAction: detail.analog.mechanismOfAction ?? '',
      therapeuticPotential: detail.analog.therapeuticPotential ?? '',
      optimizationNotes: detail.analog.optimizationNotes ?? '',
      keyDifferences: detail.analog.keyDifferences ?? '',
    });
    setEditing(true);
  };

  return (
    <DashboardLayout>
      <div className="p-6 space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3">
          <BookOpen className="h-7 w-7 text-indigo-500" />
          <div>
            <h1 className="text-2xl font-bold">Compound Encyclopedia</h1>
            <p className="text-sm text-muted-foreground">
              Persistent knowledge base for all discovered analogs — research history, ADMET data, docking results, and notes
            </p>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left: Search + List */}
          <div className="lg:col-span-1 space-y-4">
            {/* Filters */}
            <Card>
              <CardContent className="pt-4 space-y-3">
                <div className="relative">
                  <Search className="absolute left-3 top-2.5 h-4 w-4 text-muted-foreground" />
                  <Input
                    placeholder="Search compounds, SMILES, MoA..."
                    value={query}
                    onChange={(e) => setQuery(e.target.value)}
                    className="pl-9"
                  />
                </div>
                <div className="grid grid-cols-2 gap-2">
                  <Select value={patentFilter} onValueChange={setPatentFilter}>
                    <SelectTrigger className="text-xs">
                      <SelectValue placeholder="Patent status" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="all">All statuses</SelectItem>
                      <SelectItem value="patent-free">Patent-free</SelectItem>
                      <SelectItem value="patent-opportunity">Opportunity</SelectItem>
                      <SelectItem value="patented">Patented</SelectItem>
                      <SelectItem value="unknown">Unknown</SelectItem>
                    </SelectContent>
                  </Select>
                  <Select value={String(minConfidence)} onValueChange={(v) => setMinConfidence(Number(v))}>
                    <SelectTrigger className="text-xs">
                      <SelectValue placeholder="Min confidence" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="0">Any confidence</SelectItem>
                      <SelectItem value="50">≥ 50%</SelectItem>
                      <SelectItem value="70">≥ 70%</SelectItem>
                      <SelectItem value="85">≥ 85%</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <p className="text-xs text-muted-foreground">
                  {searchQuery.isLoading ? 'Loading...' : `${compounds.length} compound${compounds.length !== 1 ? 's' : ''} found`}
                </p>
              </CardContent>
            </Card>

            {/* Compound list */}
            <div className="space-y-2 max-h-[calc(100vh-320px)] overflow-y-auto pr-1">
              {searchQuery.isLoading && (
                <div className="flex items-center justify-center py-8">
                  <Loader2 className="h-5 w-5 animate-spin text-muted-foreground" />
                </div>
              )}
              {!searchQuery.isLoading && compounds.length === 0 && (
                <div className="text-center py-8 text-muted-foreground text-sm">
                  No compounds found. Generate analogs to populate the encyclopedia.
                </div>
              )}
              {compounds.map((c) => (
                <div
                  key={c.id}
                  onClick={() => { setSelectedId(c.id); setEditing(false); }}
                  className={`p-3 border rounded-lg cursor-pointer transition-all hover:border-indigo-300 ${selectedId === c.id ? 'border-indigo-500 bg-indigo-50/50 dark:bg-indigo-950/20' : ''}`}
                >
                  <div className="flex items-start justify-between gap-2">
                    <div className="flex-1 min-w-0">
                      <p className="font-medium text-sm truncate">{c.compoundName}</p>
                      <p className="text-xs text-muted-foreground font-mono truncate">{c.smiles?.slice(0, 40)}{(c.smiles?.length ?? 0) > 40 ? '…' : ''}</p>
                    </div>
                    <ChevronRight className="h-4 w-4 text-muted-foreground shrink-0 mt-0.5" />
                  </div>
                  <div className="flex items-center gap-1.5 mt-1.5 flex-wrap">
                    <Badge variant="outline" className={`text-xs h-4 px-1 ${CONFIDENCE_COLOR(c.confidenceScore)}`}>
                      {c.confidenceScore}% conf
                    </Badge>
                    <Badge variant="outline" className={`text-xs h-4 px-1 ${PATENT_COLOR[c.patentStatus] ?? ''}`}>
                      {c.patentStatus}
                    </Badge>
                    {c.approvalStatus === 'approved' && (
                      <Badge variant="outline" className="text-xs h-4 px-1 text-green-600 border-green-300">
                        <CheckCircle2 className="h-2.5 w-2.5 mr-0.5" />approved
                      </Badge>
                    )}
                  </div>
                </div>
              ))}
            </div>
          </div>

          {/* Right: Detail view */}
          <div className="lg:col-span-2">
            {!selectedId && (
              <Card className="h-full flex items-center justify-center">
                <CardContent className="text-center py-16">
                  <BookOpen className="h-12 w-12 text-muted-foreground mx-auto mb-3 opacity-30" />
                  <p className="text-muted-foreground">Select a compound to view its encyclopedia entry</p>
                </CardContent>
              </Card>
            )}

            {selectedId && detailQuery.isLoading && (
              <Card className="h-full flex items-center justify-center">
                <CardContent className="text-center py-16">
                  <Loader2 className="h-8 w-8 animate-spin text-indigo-500 mx-auto" />
                </CardContent>
              </Card>
            )}

            {selectedId && detail && (
              <div className="space-y-4">
                {/* Identity card */}
                <Card>
                  <CardHeader className="pb-3">
                    <div className="flex items-start justify-between gap-3">
                      <div className="flex-1 min-w-0">
                        <CardTitle className="text-lg">{detail.analog.compoundName}</CardTitle>
                        <CardDescription className="text-xs mt-0.5">
                          ID: {detail.analog.compoundId} · Parent: {detail.analog.parentCompound} · Gen {detail.analog.optimizationGeneration}
                        </CardDescription>
                      </div>
                      <div className="flex items-center gap-2 shrink-0">
                        {!editing ? (
                          <Button variant="outline" size="sm" onClick={startEdit}>
                            <Edit3 className="h-3.5 w-3.5 mr-1" /> Edit Notes
                          </Button>
                        ) : (
                          <>
                            <Button variant="outline" size="sm" onClick={() => setEditing(false)}>
                              <X className="h-3.5 w-3.5 mr-1" /> Cancel
                            </Button>
                            <Button size="sm" onClick={handleSaveNotes} disabled={updateNotesMutation.isPending}>
                              {updateNotesMutation.isPending ? <Loader2 className="h-3.5 w-3.5 mr-1 animate-spin" /> : <Save className="h-3.5 w-3.5 mr-1" />}
                              Save
                            </Button>
                          </>
                        )}
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    {/* SMILES */}
                    <div>
                      <p className="text-xs font-medium text-muted-foreground mb-1">SMILES</p>
                      <div className="flex items-center gap-2">
                        <code className="text-xs font-mono bg-muted px-2 py-1 rounded flex-1 truncate">{detail.analog.smiles}</code>
                        <Button variant="ghost" size="sm" className="h-7 w-7 p-0 shrink-0" onClick={() => handleCopySmiles(detail.analog.smiles)}>
                          {copied ? <CheckCircle2 className="h-3.5 w-3.5 text-green-500" /> : <Copy className="h-3.5 w-3.5" />}
                        </Button>
                      </div>
                    </div>

                    {/* Scores */}
                    <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
                      {[
                        { label: 'Confidence', value: detail.analog.confidenceScore, color: CONFIDENCE_COLOR(detail.analog.confidenceScore) },
                        { label: 'Safety', value: detail.analog.safetyScore, color: CONFIDENCE_COLOR(detail.analog.safetyScore) },
                        { label: 'Efficacy', value: detail.analog.efficacyScore, color: CONFIDENCE_COLOR(detail.analog.efficacyScore) },
                        { label: 'Drug-likeness', value: detail.analog.drugLikenessScore, color: CONFIDENCE_COLOR(detail.analog.drugLikenessScore) },
                      ].map(({ label, value, color }) => (
                        <div key={label} className={`p-2 rounded-lg border text-center ${color}`}>
                          <p className="text-xs font-medium opacity-70">{label}</p>
                          <p className="text-lg font-bold">{value}</p>
                        </div>
                      ))}
                    </div>

                    {/* External IDs */}
                    {(detail.analog.pubchemCid || detail.analog.chemblId) && (
                      <div className="flex gap-2 flex-wrap">
                        {detail.analog.pubchemCid && (
                          <a href={`https://pubchem.ncbi.nlm.nih.gov/compound/${detail.analog.pubchemCid}`} target="_blank" rel="noopener noreferrer">
                            <Badge variant="outline" className="text-xs cursor-pointer hover:bg-muted">
                              <ExternalLink className="h-2.5 w-2.5 mr-1" />PubChem {detail.analog.pubchemCid}
                            </Badge>
                          </a>
                        )}
                        {detail.analog.chemblId && (
                          <a href={`https://www.ebi.ac.uk/chembl/compound_report_card/${detail.analog.chemblId}/`} target="_blank" rel="noopener noreferrer">
                            <Badge variant="outline" className="text-xs cursor-pointer hover:bg-muted">
                              <ExternalLink className="h-2.5 w-2.5 mr-1" />ChEMBL {detail.analog.chemblId}
                            </Badge>
                          </a>
                        )}
                      </div>
                    )}

                    {/* Editable notes */}
                    {editing ? (
                      <div className="space-y-3">
                        {[
                          { key: 'mechanismOfAction', label: 'Mechanism of Action' },
                          { key: 'therapeuticPotential', label: 'Therapeutic Potential' },
                          { key: 'keyDifferences', label: 'Key Differences from Parent' },
                          { key: 'optimizationNotes', label: 'Optimization Notes' },
                        ].map(({ key, label }) => (
                          <div key={key}>
                            <p className="text-xs font-medium text-muted-foreground mb-1">{label}</p>
                            <Textarea
                              value={editFields[key] ?? ''}
                              onChange={(e) => setEditFields(f => ({ ...f, [key]: e.target.value }))}
                              rows={2}
                              className="text-sm"
                              placeholder={`Enter ${label.toLowerCase()}...`}
                            />
                          </div>
                        ))}
                      </div>
                    ) : (
                      <div className="space-y-2">
                        {[
                          { label: 'Mechanism of Action', value: detail.analog.mechanismOfAction },
                          { label: 'Therapeutic Potential', value: detail.analog.therapeuticPotential },
                          { label: 'Key Differences', value: detail.analog.keyDifferences },
                          { label: 'Optimization Notes', value: detail.analog.optimizationNotes },
                        ].filter(({ value }) => value).map(({ label, value }) => (
                          <div key={label}>
                            <p className="text-xs font-medium text-muted-foreground">{label}</p>
                            <p className="text-sm mt-0.5">{value}</p>
                          </div>
                        ))}
                      </div>
                    )}
                  </CardContent>
                </Card>

                {/* ADMET results */}
                {detail.admet.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <Activity className="h-4 w-4 text-blue-500" />
                        Latest ADMET Results
                      </CardTitle>
                    </CardHeader>
                    <CardContent>
                      {(() => {
                        const a = detail.admet[detail.admet.length - 1];
                        return (
                          <div className="grid grid-cols-3 sm:grid-cols-5 gap-2">
                            {[
                              { label: 'BBB', value: a.bbbPermeability },
                              { label: 'hERG', value: a.herg },
                              { label: 'AMES', value: a.ames },
                              { label: 'DILI', value: a.dili },
                              { label: 'QED', value: a.qed },
                              { label: 'LogP', value: a.logp },
                              { label: 'TPSA', value: a.tpsa },
                              { label: 'MW', value: a.molecularWeight },
                              { label: 'Caco-2', value: a.caco2 },
                              { label: 'CYP3A4', value: a.cyp3a4 },
                            ].filter(({ value }) => value != null).map(({ label, value }) => (
                              <div key={label} className="p-2 bg-muted/50 rounded text-center">
                                <p className="text-xs text-muted-foreground">{label}</p>
                                <p className="text-sm font-medium">{typeof value === 'number' ? (value as number).toFixed(2) : String(value ?? '')}</p>
                              </div>
                            ))}
                          </div>
                        );
                      })()}
                    </CardContent>
                  </Card>
                )}

                {/* Docking results */}
                {detail.dockings.filter(d => d.status === 'completed').length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <FlaskConical className="h-4 w-4 text-purple-500" />
                        Docking Results ({detail.dockings.filter(d => d.status === 'completed').length} targets)
                      </CardTitle>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-1">
                        {detail.dockings.filter(d => d.status === 'completed').map((d, i) => (
                          <div key={i} className="flex items-center justify-between text-xs p-2 bg-muted/50 rounded">
                            <span className="font-medium">{d.target}</span>
                            <div className="flex items-center gap-3">
                              <span className="text-muted-foreground">{d.bindingAffinity} kcal/mol</span>
                              <Badge variant="outline" className="text-xs h-4 px-1">Score: {d.dockingScore}</Badge>
                            </div>
                          </div>
                        ))}
                      </div>
                    </CardContent>
                  </Card>
                )}

                {/* Metabolites */}
                {detail.metabolites.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <Dna className="h-4 w-4 text-green-500" />
                        Predicted Metabolites ({detail.metabolites.length})
                      </CardTitle>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-1">
                        {detail.metabolites.slice(0, 8).map((m, i) => (
                          <div key={i} className="flex items-center justify-between text-xs p-2 bg-muted/50 rounded">
                            <div>
                              <span className="font-medium">{m.transformation}</span>
                              <span className="text-muted-foreground ml-2">{m.enzyme}</span>
                            </div>
                            <div className="flex items-center gap-2">
                              <Badge variant="outline" className="text-xs h-4 px-1">{m.phase}</Badge>
                              <span className="text-muted-foreground">{m.probability}</span>
                            </div>
                          </div>
                        ))}
                      </div>
                    </CardContent>
                  </Card>
                )}

                {/* Test history */}
                {detail.tests.length > 0 && (
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm font-medium flex items-center gap-2">
                        <AlertCircle className="h-4 w-4 text-orange-500" />
                        Test History ({detail.tests.length} runs)
                      </CardTitle>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-1">
                        {detail.tests.slice(-10).reverse().map((t, i) => (
                          <div key={i} className="flex items-center justify-between text-xs p-2 bg-muted/50 rounded">
                            <div className="flex items-center gap-2">
                              <Badge variant="outline" className="text-xs h-4 px-1 uppercase">{t.testType}</Badge>
                              <Badge
                                variant="outline"
                                className={`text-xs h-4 px-1 ${t.testStatus === 'completed' ? 'text-green-600 border-green-300' : t.testStatus === 'failed' ? 'text-red-600 border-red-300' : 'text-yellow-600 border-yellow-300'}`}
                              >
                                {t.testStatus}
                              </Badge>
                            </div>
                            <span className="text-muted-foreground">{new Date(t.createdAt).toLocaleDateString()}</span>
                          </div>
                        ))}
                      </div>
                    </CardContent>
                  </Card>
                )}
              </div>
            )}
          </div>
        </div>
      </div>
    </DashboardLayout>
  );
}
