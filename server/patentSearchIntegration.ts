// Use native fetch (Node 18+)

interface PatentSearchResult {
  patentNumber: string;
  title: string;
  abstract?: string;
  assignee?: string;
  filingDate?: string;
  publicationDate?: string;
  expirationDate?: string;
  status: 'active' | 'expired' | 'pending' | 'abandoned';
  url: string;
}

interface PatentSearchResponse {
  found: boolean;
  patents: PatentSearchResult[];
  freedomToOperate: boolean;
  riskLevel: 'low' | 'medium' | 'high';
  notes: string;
}

/**
 * Search USPTO for patents related to a compound
 * Note: USPTO API has limited public access. This uses the PatentsView API.
 */
export async function searchUSPTO(smiles: string, compoundName: string): Promise<PatentSearchResult[]> {
  try {
    // PatentsView API - free public access
    // Search by title/abstract keywords
    const keywords = compoundName.split(/[-\s]+/).filter(k => k.length > 3).join(' OR ');
    
    const response = await fetch('https://api.patentsview.org/patents/query', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        q: {
          _or: [
            { _text_any: { patent_title: keywords } },
            { _text_any: { patent_abstract: keywords } },
          ],
        },
        f: [
          'patent_number',
          'patent_title',
          'patent_abstract',
          'patent_date',
          'assignee_organization',
        ],
        o: { per_page: 10 },
      }),
    });

    if (!response.ok) {
      console.warn('[USPTO] API request failed:', response.statusText);
      return [];
    }

    const data: any = await response.json();
    
    if (!data.patents || data.patents.length === 0) {
      return [];
    }

    return data.patents.map((p: any) => ({
      patentNumber: p.patent_number,
      title: p.patent_title,
      abstract: p.patent_abstract,
      assignee: p.assignees?.[0]?.assignee_organization,
      publicationDate: p.patent_date,
      expirationDate: calculateExpirationDate(p.patent_date),
      status: determinePatentStatus(p.patent_date),
      url: `https://patents.google.com/patent/US${p.patent_number}`,
    }));
  } catch (error) {
    console.error('[USPTO] Search error:', error);
    return [];
  }
}

/**
 * Search EPO (Espacenet) for patents
 * Note: EPO Open Patent Services requires registration. This is a simplified version.
 */
export async function searchEPO(smiles: string, compoundName: string): Promise<PatentSearchResult[]> {
  try {
    // EPO's Open Patent Services API
    // For production, you would need to register and get API credentials
    // This is a placeholder that returns empty results
    
    console.log('[EPO] Patent search not fully implemented - requires API registration');
    console.log('[EPO] Would search for:', compoundName);
    
    // In production, you would call:
    // https://ops.epo.org/3.2/rest-services/published-data/search
    
    return [];
  } catch (error) {
    console.error('[EPO] Search error:', error);
    return [];
  }
}

/**
 * Comprehensive patent search across multiple databases
 */
export async function searchPatents(smiles: string, compoundName: string): Promise<PatentSearchResponse> {
  console.log(`[Patent Search] Searching for: ${compoundName}`);
  
  // Search both databases in parallel
  const [usptoResults, epoResults] = await Promise.all([
    searchUSPTO(smiles, compoundName),
    searchEPO(smiles, compoundName),
  ]);

  const allPatents = [...usptoResults, ...epoResults];
  
  // Analyze freedom to operate
  const activePatents = allPatents.filter(p => p.status === 'active');
  const expiringPatents = allPatents.filter(p => {
    if (!p.expirationDate) return false;
    const expDate = new Date(p.expirationDate);
    const now = new Date();
    const yearsUntilExpiration = (expDate.getTime() - now.getTime()) / (1000 * 60 * 60 * 24 * 365);
    return yearsUntilExpiration > 0 && yearsUntilExpiration < 2;
  });

  let freedomToOperate = true;
  let riskLevel: 'low' | 'medium' | 'high' = 'low';
  let notes = '';

  if (activePatents.length === 0) {
    notes = 'No active patents found. Compound appears to be patent-free.';
    riskLevel = 'low';
    freedomToOperate = true;
  } else if (activePatents.length <= 2) {
    notes = `${activePatents.length} active patent(s) found. Review required for freedom to operate.`;
    riskLevel = 'medium';
    freedomToOperate = false;
  } else {
    notes = `${activePatents.length} active patents found. High risk - extensive patent landscape.`;
    riskLevel = 'high';
    freedomToOperate = false;
  }

  if (expiringPatents.length > 0) {
    notes += ` ${expiringPatents.length} patent(s) expiring within 2 years.`;
  }

  return {
    found: allPatents.length > 0,
    patents: allPatents,
    freedomToOperate,
    riskLevel,
    notes,
  };
}

/**
 * Calculate patent expiration date (20 years from filing for US patents)
 */
function calculateExpirationDate(publicationDate: string): string {
  try {
    const pubDate = new Date(publicationDate);
    // US patents typically expire 20 years from filing date
    // Using publication date as approximation
    pubDate.setFullYear(pubDate.getFullYear() + 20);
    return pubDate.toISOString().split('T')[0];
  } catch {
    return '';
  }
}

/**
 * Determine patent status based on dates
 */
function determinePatentStatus(publicationDate: string): 'active' | 'expired' | 'pending' | 'abandoned' {
  try {
    const pubDate = new Date(publicationDate);
    const expDate = new Date(pubDate);
    expDate.setFullYear(expDate.getFullYear() + 20);
    
    const now = new Date();
    
    if (now > expDate) {
      return 'expired';
    } else {
      return 'active';
    }
  } catch {
    return 'pending';
  }
}

/**
 * Monitor patents for expiration and send alerts
 */
export async function checkExpiringPatents(analogId: number, patents: PatentSearchResult[]): Promise<string[]> {
  const alerts: string[] = [];
  const now = new Date();
  
  for (const patent of patents) {
    if (!patent.expirationDate) continue;
    
    const expDate = new Date(patent.expirationDate);
    const yearsUntilExpiration = (expDate.getTime() - now.getTime()) / (1000 * 60 * 60 * 24 * 365);
    
    if (yearsUntilExpiration > 0 && yearsUntilExpiration < 2) {
      alerts.push(
        `Patent ${patent.patentNumber} expires in ${Math.round(yearsUntilExpiration * 12)} months. ` +
        `Consider monitoring for freedom to operate opportunities.`
      );
    }
  }
  
  return alerts;
}
