import { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';
import { ChevronDown, ChevronUp, ExternalLink, Loader2 } from 'lucide-react';
import { trpc } from '@/lib/trpc';

interface EnrichmentCardProps {
  smiles: string;
  source: 'pubchem' | 'chembl';
  title: string;
  description: string;
}

export function EnrichmentCard({ smiles, source, title, description }: EnrichmentCardProps) {
  const [isOpen, setIsOpen] = useState(false);
  const [hasEnriched, setHasEnriched] = useState(false);

  const pubchemMutation = trpc.enrichment.pubchem.useMutation();
  const chemblMutation = trpc.enrichment.chembl.useMutation();

  const mutation = source === 'pubchem' ? pubchemMutation : chemblMutation;

  const handleEnrich = async () => {
    setHasEnriched(true);
    setIsOpen(true);
    
    if (source === 'pubchem') {
      await pubchemMutation.mutateAsync({ smiles, cid: undefined });
    } else {
      await chemblMutation.mutateAsync({ smiles, chemblId: undefined, includeActivityData: true });
    }
  };

  const data = mutation.data;
  const isLoading = mutation.isPending;
  const error = mutation.error;

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <div>
            <CardTitle className="flex items-center gap-2">
              {title}
              {source === 'pubchem' && (
                <a
                  href={`https://pubchem.ncbi.nlm.nih.gov/#query=${encodeURIComponent(smiles)}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 hover:text-blue-800"
                >
                  <ExternalLink className="h-4 w-4" />
                </a>
              )}
              {source === 'chembl' && (
                <a
                  href={`https://www.ebi.ac.uk/chembl/`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 hover:text-blue-800"
                >
                  <ExternalLink className="h-4 w-4" />
                </a>
              )}
            </CardTitle>
            <CardDescription>{description}</CardDescription>
          </div>
          {!hasEnriched && (
            <Button onClick={handleEnrich} disabled={isLoading}>
              {isLoading ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Loading...
                </>
              ) : (
                'Enrich'
              )}
            </Button>
          )}
        </div>
      </CardHeader>

      {hasEnriched && (
        <Collapsible open={isOpen} onOpenChange={setIsOpen}>
          <CardContent>
            <CollapsibleTrigger asChild>
              <Button variant="ghost" className="w-full justify-between">
                {isOpen ? 'Hide Details' : 'Show Details'}
                {isOpen ? <ChevronUp className="h-4 w-4" /> : <ChevronDown className="h-4 w-4" />}
              </Button>
            </CollapsibleTrigger>

            <CollapsibleContent className="mt-4 space-y-4">
              {isLoading && (
                <div className="flex items-center justify-center py-8">
                  <Loader2 className="h-8 w-8 animate-spin text-blue-600" />
                </div>
              )}

              {error && (
                <div className="rounded-lg bg-red-50 p-4 text-red-800">
                  <p className="font-semibold">Error loading data</p>
                  <p className="text-sm">{error.message}</p>
                </div>
              )}

              {data && !data.success && (
                <div className="rounded-lg bg-yellow-50 p-4 text-yellow-800">
                  <p className="font-semibold">No data found</p>
                  <p className="text-sm">{data.error || 'Compound not found in database'}</p>
                </div>
              )}

              {data && data.success && source === 'pubchem' && data.data && 'cid' in data.data && (
                <div className="space-y-4">
                  <div>
                    <h4 className="font-semibold text-sm text-gray-700 mb-2">Compound Information</h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>
                        <span className="text-gray-600">PubChem CID:</span>
                        <span className="ml-2 font-medium">{data.data.cid}</span>
                      </div>
                      {data.data.iupacName && (
                        <div>
                          <span className="text-gray-600">IUPAC Name:</span>
                          <span className="ml-2 font-medium">{data.data.iupacName}</span>
                        </div>
                      )}
                      {data.data.molecularFormula && (
                        <div>
                          <span className="text-gray-600">Formula:</span>
                          <span className="ml-2 font-medium">{data.data.molecularFormula}</span>
                        </div>
                      )}
                      {data.data.molecularWeight && (
                        <div>
                          <span className="text-gray-600">Molecular Weight:</span>
                          <span className="ml-2 font-medium">{data.data.molecularWeight.toFixed(2)} g/mol</span>
                        </div>
                      )}
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold text-sm text-gray-700 mb-2">Physicochemical Properties</h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      {data.data.xlogp !== undefined && (
                        <div>
                          <span className="text-gray-600">XLogP:</span>
                          <span className="ml-2 font-medium">{data.data.xlogp.toFixed(2)}</span>
                        </div>
                      )}
                      {data.data.tpsa !== undefined && (
                        <div>
                          <span className="text-gray-600">TPSA:</span>
                          <span className="ml-2 font-medium">{data.data.tpsa.toFixed(2)} Ų</span>
                        </div>
                      )}
                      {data.data.hBondDonorCount !== undefined && (
                        <div>
                          <span className="text-gray-600">H-Bond Donors:</span>
                          <span className="ml-2 font-medium">{data.data.hBondDonorCount}</span>
                        </div>
                      )}
                      {data.data.hBondAcceptorCount !== undefined && (
                        <div>
                          <span className="text-gray-600">H-Bond Acceptors:</span>
                          <span className="ml-2 font-medium">{data.data.hBondAcceptorCount}</span>
                        </div>
                      )}
                      {data.data.rotatableBondCount !== undefined && (
                        <div>
                          <span className="text-gray-600">Rotatable Bonds:</span>
                          <span className="ml-2 font-medium">{data.data.rotatableBondCount}</span>
                        </div>
                      )}
                      {data.data.complexity !== undefined && (
                        <div>
                          <span className="text-gray-600">Complexity:</span>
                          <span className="ml-2 font-medium">{data.data.complexity.toFixed(0)}</span>
                        </div>
                      )}
                    </div>
                  </div>

                  {data.data.synonyms && data.data.synonyms.length > 0 && (
                    <div>
                      <h4 className="font-semibold text-sm text-gray-700 mb-2">Synonyms</h4>
                      <div className="flex flex-wrap gap-1">
                        {data.data.synonyms.slice(0, 8).map((synonym, idx) => (
                          <Badge key={idx} variant="secondary" className="text-xs">
                            {synonym}
                          </Badge>
                        ))}
                      </div>
                    </div>
                  )}

                  {data.data.bioactivity && (
                    <div>
                      <h4 className="font-semibold text-sm text-gray-700 mb-2">Bioactivity Summary</h4>
                      <div className="grid grid-cols-2 gap-2 text-sm">
                        <div>
                          <span className="text-gray-600">Total Assays:</span>
                          <span className="ml-2 font-medium">{data.data.bioactivity.assayCount || 0}</span>
                        </div>
                        <div>
                          <span className="text-gray-600">Active Assays:</span>
                          <span className="ml-2 font-medium">{data.data.bioactivity.activeAssayCount || 0}</span>
                        </div>
                      </div>
                    </div>
                  )}

                  {'description' in data.data && data.data.description && (
                    <div>
                      <h4 className="font-semibold text-sm text-gray-700 mb-2">Description</h4>
                      <p className="text-sm text-gray-600">{'description' in data.data ? data.data.description : ''}</p>
                    </div>
                  )}
                </div>
              )}

              {data && data.success && source === 'chembl' && data.data && 'chemblId' in data.data && (
                <div className="space-y-4">
                  <div>
                    <h4 className="font-semibold text-sm text-gray-700 mb-2">Compound Information</h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>
                        <span className="text-gray-600">ChEMBL ID:</span>
                        <span className="ml-2 font-medium">{data.data.chemblId}</span>
                      </div>
                      {data.data.preferredName && (
                        <div>
                          <span className="text-gray-600">Name:</span>
                          <span className="ml-2 font-medium">{data.data.preferredName}</span>
                        </div>
                      )}
                      {data.data.molecularFormula && (
                        <div>
                          <span className="text-gray-600">Formula:</span>
                          <span className="ml-2 font-medium">{data.data.molecularFormula}</span>
                        </div>
                      )}
                      {data.data.molecularWeight && (
                        <div>
                          <span className="text-gray-600">Molecular Weight:</span>
                          <span className="ml-2 font-medium">{data.data.molecularWeight.toFixed(2)} g/mol</span>
                        </div>
                      )}
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold text-sm text-gray-700 mb-2">Drug-likeness Properties</h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      {data.data.alogp !== undefined && (
                        <div>
                          <span className="text-gray-600">ALogP:</span>
                          <span className="ml-2 font-medium">{data.data.alogp.toFixed(2)}</span>
                        </div>
                      )}
                      {data.data.psa !== undefined && (
                        <div>
                          <span className="text-gray-600">PSA:</span>
                          <span className="ml-2 font-medium">{data.data.psa.toFixed(2)} Ų</span>
                        </div>
                      )}
                      {data.data.hba !== undefined && (
                        <div>
                          <span className="text-gray-600">H-Bond Acceptors:</span>
                          <span className="ml-2 font-medium">{data.data.hba}</span>
                        </div>
                      )}
                      {data.data.hbd !== undefined && (
                        <div>
                          <span className="text-gray-600">H-Bond Donors:</span>
                          <span className="ml-2 font-medium">{data.data.hbd}</span>
                        </div>
                      )}
                      {data.data.rotatableBonds !== undefined && (
                        <div>
                          <span className="text-gray-600">Rotatable Bonds:</span>
                          <span className="ml-2 font-medium">{data.data.rotatableBonds}</span>
                        </div>
                      )}
                      {data.data.numRo5Violations !== undefined && (
                        <div>
                          <span className="text-gray-600">Ro5 Violations:</span>
                          <span className="ml-2 font-medium">{data.data.numRo5Violations}</span>
                        </div>
                      )}
                    </div>
                  </div>

                  {data.data.activities && data.data.activities.length > 0 && (
                    <div>
                      <h4 className="font-semibold text-sm text-gray-700 mb-2">Bioactivity Data</h4>
                      <div className="space-y-2 max-h-64 overflow-y-auto">
                        {data.data.activities.slice(0, 10).map((activity, idx) => (
                          <div key={idx} className="border rounded-lg p-3 text-sm">
                            <div className="font-medium text-gray-900">{activity.targetName}</div>
                            <div className="text-gray-600 text-xs mt-1">
                              {activity.activityType}: {activity.activityRelation} {activity.activityValue} {activity.activityUnits}
                              {activity.pChEMBL && ` (pChEMBL: ${activity.pChEMBL.toFixed(2)})`}
                            </div>
                            <div className="text-gray-500 text-xs mt-1">
                              {activity.assayType} • {activity.organism}
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </CollapsibleContent>
          </CardContent>
        </Collapsible>
      )}
    </Card>
  );
}
