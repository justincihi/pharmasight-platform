/**
 * PubChem Plugin
 * 
 * Fetches compound metadata, properties, and bioactivity data from PubChem
 * https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
 */

import { DataSourcePlugin, PluginMetadata, PluginConfig, PluginResult } from './base';

export interface PubChemQuery {
  smiles?: string;
  cid?: number; // PubChem Compound ID
  name?: string;
  properties?: string[]; // e.g., ['MolecularWeight', 'XLogP', 'TPSA']
}

export interface PubChemCompoundData {
  cid: number;
  iupacName?: string;
  molecularFormula?: string;
  molecularWeight?: number;
  canonicalSMILES?: string;
  isomericSMILES?: string;
  inchi?: string;
  inchiKey?: string;
  xlogp?: number;
  tpsa?: number; // Topological Polar Surface Area
  complexity?: number;
  hBondDonorCount?: number;
  hBondAcceptorCount?: number;
  rotatableBondCount?: number;
  heavyAtomCount?: number;
  charge?: number;
  synonyms?: string[];
  description?: string;
  bioactivity?: {
    assayCount?: number;
    activeAssayCount?: number;
    targets?: string[];
  };
}

export class PubChemPlugin implements DataSourcePlugin<PubChemQuery, PubChemCompoundData> {
  readonly metadata: PluginMetadata = {
    name: 'pubchem',
    version: '1.0.0',
    description: 'PubChem compound data source',
    author: 'PharmaSight',
    license: 'MIT',
  };

  config: PluginConfig = {
    enabled: true,
    baseUrl: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug',
    timeout: 30000,
    retryAttempts: 3,
    cacheEnabled: true,
    cacheTTL: 86400000, // 24 hours
  };

  private cache: Map<string, { data: PubChemCompoundData; timestamp: number }> = new Map();

  async initialize(): Promise<void> {
    console.log('[PubChem Plugin] Initialized');
  }

  async isAvailable(): Promise<boolean> {
    try {
      const response = await fetch(`${this.config.baseUrl}/compound/cid/2244/property/MolecularWeight/JSON`, {
        signal: AbortSignal.timeout(5000),
      });
      return response.ok;
    } catch {
      return false;
    }
  }

  async getRateLimit(): Promise<{ limit: number; remaining: number; reset: Date }> {
    // PubChem doesn't have strict rate limits for reasonable use
    // But recommends no more than 5 requests per second
    return {
      limit: 300, // 5 req/sec * 60 sec
      remaining: 300,
      reset: new Date(Date.now() + 60000),
    };
  }

  async validate(input: PubChemQuery): Promise<{ valid: boolean; error?: string }> {
    if (!input.smiles && !input.cid && !input.name) {
      return { valid: false, error: 'Must provide smiles, cid, or name' };
    }
    return { valid: true };
  }

  async run(input: PubChemQuery): Promise<PluginResult<PubChemCompoundData>> {
    const startTime = Date.now();
    
    try {
      // Check cache first
      const cacheKey = JSON.stringify(input);
      if (this.config.cacheEnabled) {
        const cached = this.cache.get(cacheKey);
        if (cached && Date.now() - cached.timestamp < this.config.cacheTTL) {
          return {
            success: true,
            data: cached.data,
            metadata: {
              executionTime: Date.now() - startTime,
              source: 'pubchem',
              cached: true,
            },
          };
        }
      }

      // Fetch from PubChem
      let cid: number | undefined = input.cid;
      
      // If no CID provided, search by SMILES or name
      if (!cid) {
        if (input.smiles) {
          cid = await this.searchBySMILES(input.smiles);
        } else if (input.name) {
          cid = await this.searchByName(input.name);
        }
      }

      if (!cid) {
        return {
          success: false,
          error: 'Compound not found in PubChem',
          metadata: {
            executionTime: Date.now() - startTime,
            source: 'pubchem',
          },
        };
      }

      // Fetch compound data
      const data = await this.fetchCompoundData(cid);

      // Cache the result
      if (this.config.cacheEnabled) {
        this.cache.set(cacheKey, { data, timestamp: Date.now() });
      }

      return {
        success: true,
        data,
        metadata: {
          executionTime: Date.now() - startTime,
          source: 'pubchem',
          cached: false,
        },
      };
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
        metadata: {
          executionTime: Date.now() - startTime,
          source: 'pubchem',
        },
      };
    }
  }

  private async searchBySMILES(smiles: string): Promise<number | undefined> {
    const url = `${this.config.baseUrl}/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`;
    const response = await fetch(url, { signal: AbortSignal.timeout(this.config.timeout) });
    
    if (!response.ok) return undefined;
    
    const json = await response.json();
    return json.IdentifierList?.CID?.[0];
  }

  private async searchByName(name: string): Promise<number | undefined> {
    const url = `${this.config.baseUrl}/compound/name/${encodeURIComponent(name)}/cids/JSON`;
    const response = await fetch(url, { signal: AbortSignal.timeout(this.config.timeout) });
    
    if (!response.ok) return undefined;
    
    const json = await response.json();
    return json.IdentifierList?.CID?.[0];
  }

  private async fetchCompoundData(cid: number): Promise<PubChemCompoundData> {
    // Fetch properties
    const propsUrl = `${this.config.baseUrl}/compound/cid/${cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,XLogP,TPSA,Complexity,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,Charge/JSON`;
    const propsResponse = await fetch(propsUrl, { signal: AbortSignal.timeout(this.config.timeout) });
    const propsJson = await propsResponse.json();
    const props = propsJson.PropertyTable?.Properties?.[0] || {};

    // Fetch synonyms
    const synonymsUrl = `${this.config.baseUrl}/compound/cid/${cid}/synonyms/JSON`;
    const synonymsResponse = await fetch(synonymsUrl, { signal: AbortSignal.timeout(this.config.timeout) });
    const synonymsJson = await synonymsResponse.json();
    const synonyms = synonymsJson.InformationList?.Information?.[0]?.Synonym || [];

    // Fetch description (if available)
    let description: string | undefined;
    try {
      const descUrl = `${this.config.baseUrl}/compound/cid/${cid}/description/JSON`;
      const descResponse = await fetch(descUrl, { signal: AbortSignal.timeout(this.config.timeout) });
      const descJson = await descResponse.json();
      description = descJson.InformationList?.Information?.[0]?.Description;
    } catch {
      // Description not always available
    }

    // Fetch bioactivity summary
    let bioactivity: PubChemCompoundData['bioactivity'];
    try {
      const bioUrl = `${this.config.baseUrl}/compound/cid/${cid}/assaysummary/JSON`;
      const bioResponse = await fetch(bioUrl, { signal: AbortSignal.timeout(this.config.timeout) });
      const bioJson = await bioResponse.json();
      const summary = bioJson.Table?.Row?.[0] || {};
      bioactivity = {
        assayCount: parseInt(summary.assaycount) || 0,
        activeAssayCount: parseInt(summary.activeassaycount) || 0,
      };
    } catch {
      // Bioactivity not always available
    }

    return {
      cid,
      iupacName: synonyms[0],
      molecularFormula: props.MolecularFormula,
      molecularWeight: props.MolecularWeight,
      canonicalSMILES: props.CanonicalSMILES,
      isomericSMILES: props.IsomericSMILES,
      inchi: props.InChI,
      inchiKey: props.InChIKey,
      xlogp: props.XLogP,
      tpsa: props.TPSA,
      complexity: props.Complexity,
      hBondDonorCount: props.HBondDonorCount,
      hBondAcceptorCount: props.HBondAcceptorCount,
      rotatableBondCount: props.RotatableBondCount,
      heavyAtomCount: props.HeavyAtomCount,
      charge: props.Charge,
      synonyms: synonyms.slice(0, 10), // Limit to 10 synonyms
      description,
      bioactivity,
    };
  }

  summarizeResults(result: PluginResult<PubChemCompoundData>): string | object {
    if (!result.success || !result.data) {
      return `PubChem lookup failed: ${result.error}`;
    }

    const data = result.data;
    return {
      summary: `Found compound CID ${data.cid}: ${data.iupacName || 'Unknown name'}`,
      properties: {
        'Molecular Formula': data.molecularFormula,
        'Molecular Weight': data.molecularWeight ? `${data.molecularWeight.toFixed(2)} g/mol` : 'N/A',
        'XLogP': data.xlogp?.toFixed(2),
        'TPSA': data.tpsa ? `${data.tpsa.toFixed(2)} Ų` : 'N/A',
        'H-Bond Donors': data.hBondDonorCount,
        'H-Bond Acceptors': data.hBondAcceptorCount,
      },
      bioactivity: data.bioactivity ? {
        'Total Assays': data.bioactivity.assayCount,
        'Active Assays': data.bioactivity.activeAssayCount,
      } : 'No bioactivity data',
      cached: result.metadata?.cached ? 'Yes' : 'No',
    };
  }

  async cleanup(): Promise<void> {
    this.cache.clear();
    console.log('[PubChem Plugin] Cleaned up');
  }
}
