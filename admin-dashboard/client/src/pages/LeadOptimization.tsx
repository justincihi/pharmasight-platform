import { useState } from 'react';
import { trpc } from '@/lib/trpc';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { GitBranch, TrendingUp, Search, Filter, ArrowRight } from 'lucide-react';
import { Link } from 'wouter';

export function LeadOptimization() {
  const [searchTerm, setSearchTerm] = useState('');
  const [generationFilter, setGenerationFilter] = useState<string>('all');
  
  const { data: analogs, isLoading } = trpc.analog.list.useQuery({
    limit: 100,
    offset: 0,
  });

  // Group analogs by optimization lineage
  const groupedAnalogs = analogs?.analogs.reduce((acc: Record<number, any[]>, analog: any) => {
    const parentId = analog.parentAnalogId || analog.id;
    if (!acc[parentId]) {
      acc[parentId] = [];
    }
    acc[parentId].push(analog);
    return acc;
  }, {} as Record<number, any[]>) || {};

  // Filter lineages that have optimizations
  const optimizationLineages = Object.entries(groupedAnalogs)
    .filter(([_, lineage]) => (lineage as any[]).length > 1 || (lineage as any[])[0].optimizationGeneration > 1)
    .map(([parentId, lineage]) => ({
      parentId: parseInt(parentId),
      lineage: (lineage as any[]).sort((a: any, b: any) => (a.optimizationGeneration || 1) - (b.optimizationGeneration || 1)),
    }));

  const filteredLineages = optimizationLineages.filter(({ lineage }) => {
    const matchesSearch = lineage.some((analog: any) =>
      analog.compoundName.toLowerCase().includes(searchTerm.toLowerCase()) ||
      analog.compoundId.toLowerCase().includes(searchTerm.toLowerCase())
    );
    
    const matchesGeneration = generationFilter === 'all' || 
      lineage.some((analog: any) => (analog.optimizationGeneration || 1).toString() === generationFilter);
    
    return matchesSearch && matchesGeneration;
  });

  if (isLoading) {
    return (
      <div className="container py-8">
        <div className="flex items-center justify-center py-12">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto mb-4"></div>
            <p className="text-muted-foreground">Loading optimization history...</p>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="container py-8 space-y-6">
      <div>
        <h1 className="text-3xl font-bold flex items-center gap-2 mb-2">
          <GitBranch className="h-8 w-8" />
          Lead Optimization Tracker
        </h1>
        <p className="text-muted-foreground">
          Track molecular optimization lineages and compare generations
        </p>
      </div>

      <Card>
        <CardHeader>
          <CardTitle>Filters</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex gap-4">
            <div className="flex-1">
              <div className="relative">
                <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="Search by compound name or ID..."
                  value={searchTerm}
                  onChange={(e) => setSearchTerm(e.target.value)}
                  className="pl-10"
                />
              </div>
            </div>
            <Select value={generationFilter} onValueChange={setGenerationFilter}>
              <SelectTrigger className="w-48">
                <SelectValue placeholder="Filter by generation" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="all">All Generations</SelectItem>
                <SelectItem value="1">Generation 1 (Original)</SelectItem>
                <SelectItem value="2">Generation 2</SelectItem>
                <SelectItem value="3">Generation 3</SelectItem>
                <SelectItem value="4">Generation 4+</SelectItem>
              </SelectContent>
            </Select>
          </div>
        </CardContent>
      </Card>

      {filteredLineages.length === 0 ? (
        <Card>
          <CardContent className="py-12">
            <div className="text-center">
              <GitBranch className="h-12 w-12 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-lg font-semibold mb-2">No Optimization Lineages Found</h3>
              <p className="text-muted-foreground mb-4">
                Start optimizing analogs to track their improvement over generations
              </p>
              <Button asChild>
                <Link href="/discoveries">
                  <ArrowRight className="h-4 w-4 mr-2" />
                  View Discoveries
                </Link>
              </Button>
            </div>
          </CardContent>
        </Card>
      ) : (
        <div className="space-y-6">
          {filteredLineages.map(({ parentId, lineage }) => (
            <Card key={parentId}>
              <CardHeader>
                <div className="flex items-start justify-between">
                  <div>
                    <CardTitle className="flex items-center gap-2">
                      Optimization Lineage: {lineage[0].compoundName}
                      <Badge variant="outline">{lineage.length} Generation{lineage.length > 1 ? 's' : ''}</Badge>
                    </CardTitle>
                    <CardDescription>
                      {lineage[0].optimizationTarget && (
                        <span>Target: {lineage[0].optimizationTarget}</span>
                      )}
                    </CardDescription>
                  </div>
                </div>
              </CardHeader>
              <CardContent>
                <div className="space-y-4">
                  {lineage.map((analog: any, index: number) => (
                    <div key={analog.id} className="flex items-center gap-4">
                      <div className="flex items-center gap-2 min-w-[120px]">
                        <Badge variant={index === 0 ? "default" : "secondary"}>
                          Gen {analog.optimizationGeneration || 1}
                        </Badge>
                        {index < lineage.length - 1 && (
                          <ArrowRight className="h-4 w-4 text-muted-foreground" />
                        )}
                      </div>
                      
                      <div className="flex-1 border rounded-lg p-4">
                        <div className="flex items-start justify-between">
                          <div className="flex-1">
                            <Link href={`/analog/${analog.id}`}>
                              <h4 className="font-semibold hover:underline cursor-pointer">
                                {analog.compoundName}
                              </h4>
                            </Link>
                            <p className="text-sm text-muted-foreground">{analog.compoundId}</p>
                            {analog.optimizationNotes && (
                              <p className="text-sm mt-1">{analog.optimizationNotes}</p>
                            )}
                          </div>
                          
                          <div className="flex gap-4 ml-4">
                            <div className="text-center">
                              <div className="text-xs text-muted-foreground">Safety</div>
                              <div className="text-lg font-bold">{analog.safetyScore}</div>
                            </div>
                            <div className="text-center">
                              <div className="text-xs text-muted-foreground">Efficacy</div>
                              <div className="text-lg font-bold">{analog.efficacyScore}</div>
                            </div>
                            <div className="text-center">
                              <div className="text-xs text-muted-foreground">Confidence</div>
                              <div className="text-lg font-bold">{analog.confidenceScore}%</div>
                            </div>
                          </div>
                        </div>
                        
                        {index > 0 && (
                          <div className="mt-3 pt-3 border-t flex items-center gap-2 text-sm">
                            <TrendingUp className="h-4 w-4 text-green-600" />
                            <span className="text-muted-foreground">Improvements:</span>
                            <span className="text-green-600 font-medium">
                              Safety +{analog.safetyScore - lineage[0].safetyScore}, 
                              Efficacy +{analog.efficacyScore - lineage[0].efficacyScore}
                            </span>
                          </div>
                        )}
                      </div>
                    </div>
                  ))}
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
      )}
    </div>
  );
}

export default LeadOptimization;
