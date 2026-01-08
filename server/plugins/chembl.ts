/**
 * ChEMBL Plugin
 * 
 * Fetches bioactivity data, target information, and compound data from ChEMBL
 * https://www.ebi.ac.uk/chembl/api/data/docs
 */

import { DataSourcePlugin, PluginMetadata, PluginConfig, PluginResult } from './base';

export interface ChEMBLQuery {
  smiles?: string;
  chemblId?: string;
  targetName?: string;
  includeActivityData?: boolean;
  includeSimilarCompounds?: boolean;
  similarityThreshold?: number; // 0-100, default 70
}

export interface ChEMBLActivityData {
  assayId: string;
  assayType: string;
  targetName: string;
  targetType: string;
  organism: string;
  activityType: string; // IC50, Ki, EC50, etc.
  activityValue: number;
  activityUnits: string;
  activityRelation: string; // =, <, >, etc.
  pChEMBL?: number; // Standardized activity value
}

export interface ChEMBLCompoundData {
  chemblId: string;
  preferredName?: string;
  smiles?: string;
  inchi?: string;
  inchiKey?: string;
  molecularWeight?: number;
  alogp?: number;
  psa?: number;
  hba?: number; // H-bond acceptors
  hbd?: number; // H-bond donors
  rotatableBonds?: number;
  molecularFormula?: string;
  numRo5Violations?: number; // Lipinski's Rule of 5 violations
  activities?: ChEMBLActivityData[];
  similarCompounds?: Array<{
    chemblId: string;
    similarity: number;
    smiles: string;
  }>;
}

export class ChEMBLPlugin implements DataSourcePlugin<ChEMBLQuery, ChEMBLCompoundData> {
  readonly metadata: PluginMetadata = {
    name: 'chembl',
    version: '1.0.0',
    description: 'ChEMBL bioactivity and compound data source',
    author: 'PharmaSight',
    license: 'MIT',
  };

  config: PluginConfig = {
    enabled: true,
    baseUrl: 'https://www.ebi.ac.uk/chembl/api/data',
    timeout: 30000,
    retryAttempts: 3,
    cacheEnabled: true,
    cacheTTL: 86400000, // 24 hours
  };

  private cache: Map<string, { data: ChEMBLCompoundData; timestamp: number }> = new Map();

  async initialize(): Promise<void> {
    console.log('[ChEMBL Plugin] Initialized');
  }

  async isAvailable(): Promise<boolean> {
    try {
      const response = await fetch(`${this.config.baseUrl}/status.json`, {
        signal: AbortSignal.timeout(5000),
      });
      return response.ok;
    } catch {
      return false;
    }
  }

  async validate(input: ChEMBLQuery): Promise<{ valid: boolean; error?: string }> {
    if (!input.smiles && !input.chemblId && !input.targetName) {
      return { valid: false, error: 'Must provide smiles, chemblId, or targetName' };
    }
    return { valid: true };
  }

  async run(input: ChEMBLQuery): Promise<PluginResult<ChEMBLCompoundData>> {
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
              source: 'chembl',
              cached: true,
            },
          };
        }
      }

      // Fetch from ChEMBL
      let chemblId: string | undefined = input.chemblId;
      
      // If no ChEMBL ID provided, search by SMILES
      if (!chemblId && input.smiles) {
        chemblId = await this.searchBySMILES(input.smiles);
      }

      if (!chemblId) {
        return {
          success: false,
          error: 'Compound not found in ChEMBL',
          metadata: {
            executionTime: Date.now() - startTime,
            source: 'chembl',
          },
        };
      }

      // Fetch compound data
      const data = await this.fetchCompoundData(chemblId, input);

      // Cache the result
      if (this.config.cacheEnabled) {
        this.cache.set(cacheKey, { data, timestamp: Date.now() });
      }

      return {
        success: true,
        data,
        metadata: {
          executionTime: Date.now() - startTime,
          source: 'chembl',
          cached: false,
        },
      };
    } catch (error) {
      return {
        success: false,
        error: error instanceof Error ? error.message : 'Unknown error',
        metadata: {
          executionTime: Date.now() - startTime,
          source: 'chembl',
        },
      };
    }
  }

  private async searchBySMILES(smiles: string): Promise<string | undefined> {
    // ChEMBL similarity search
    const url = `${this.config.baseUrl}/similarity/${encodeURIComponent(smiles)}/70.json`;
    const response = await fetch(url, { signal: AbortSignal.timeout(this.config.timeout) });
    
    if (!response.ok) return undefined;
    
    const json = await response.json();
    return json.molecules?.[0]?.molecule_chembl_id;
  }

  private async fetchCompoundData(chemblId: string, query: ChEMBLQuery): Promise<ChEMBLCompoundData> {
    // Fetch molecule data
    const molUrl = `${this.config.baseUrl}/molecule/${chemblId}.json`;
    const molResponse = await fetch(molUrl, { signal: AbortSignal.timeout(this.config.timeout) });
    const molJson = await molResponse.json();

    const data: ChEMBLCompoundData = {
      chemblId,
      preferredName: molJson.pref_name,
      smiles: molJson.molecule_structures?.canonical_smiles,
      inchi: molJson.molecule_structures?.standard_inchi,
      inchiKey: molJson.molecule_structures?.standard_inchi_key,
      molecularWeight: molJson.molecule_properties?.full_mwt,
      alogp: molJson.molecule_properties?.alogp,
      psa: molJson.molecule_properties?.psa,
      hba: molJson.molecule_properties?.hba,
      hbd: molJson.molecule_properties?.hbd,
      rotatableBonds: molJson.molecule_properties?.rtb,
      molecularFormula: molJson.molecule_properties?.full_molformula,
      numRo5Violations: molJson.molecule_properties?.num_ro5_violations,
    };

    // Fetch activity data if requested
    if (query.includeActivityData) {
      data.activities = await this.fetchActivityData(chemblId);
    }

    // Fetch similar compounds if requested
    if (query.includeSimilarCompounds && data.smiles) {
      data.similarCompounds = await this.fetchSimilarCompounds(
        data.smiles,
        query.similarityThreshold || 70
      );
    }

    return data;
  }

  private async fetchActivityData(chemblId: string): Promise<ChEMBLActivityData[]> {
    try {
      const url = `${this.config.baseUrl}/activity.json?molecule_chembl_id=${chemblId}&limit=100`;
      const response = await fetch(url, { signal: AbortSignal.timeout(this.config.timeout) });
      const json = await response.json();

      return (json.activities || []).map((activity: any) => ({
        assayId: activity.assay_chembl_id,
        assayType: activity.assay_type,
        targetName: activity.target_pref_name,
        targetType: activity.target_organism,
        organism: activity.target_organism,
        activityType: activity.standard_type,
        activityValue: parseFloat(activity.standard_value),
        activityUnits: activity.standard_units,
        activityRelation: activity.standard_relation,
        pChEMBL: activity.pchembl_value ? parseFloat(activity.pchembl_value) : undefined,
      }));
    } catch {
      return [];
    }
  }

  private async fetchSimilarCompounds(
    smiles: string,
    threshold: number
  ): Promise<Array<{ chemblId: string; similarity: number; smiles: string }>> {
    try {
      const url = `${this.config.baseUrl}/similarity/${encodeURIComponent(smiles)}/${threshold}.json`;
      const response = await fetch(url, { signal: AbortSignal.timeout(this.config.timeout) });
      const json = await response.json();

      return (json.molecules || []).slice(0, 10).map((mol: any) => ({
        chemblId: mol.molecule_chembl_id,
        similarity: mol.similarity,
        smiles: mol.molecule_structures?.canonical_smiles,
      }));
    } catch {
      return [];
    }
  }

  summarizeResults(result: PluginResult<ChEMBLCompoundData>): string | object {
    if (!result.success || !result.data) {
      return `ChEMBL lookup failed: ${result.error}`;
    }

    const data = result.data;
    return {
      summary: `Found compound ${data.chemblId}: ${data.preferredName || 'Unknown name'}`,
      properties: {
        'Molecular Formula': data.molecularFormula,
        'Molecular Weight': data.molecularWeight ? `${data.molecularWeight.toFixed(2)} g/mol` : 'N/A',
        'ALogP': data.alogp?.toFixed(2),
        'PSA': data.psa ? `${data.psa.toFixed(2)} Ų` : 'N/A',
        'H-Bond Donors': data.hbd,
        'H-Bond Acceptors': data.hba,
        'Ro5 Violations': data.numRo5Violations,
      },
      bioactivity: data.activities ? {
        'Total Activities': data.activities.length,
        'Unique Targets': new Set(data.activities.map(a => a.targetName)).size,
        'Activity Types': Array.from(new Set(data.activities.map(a => a.activityType))).join(', '),
      } : 'No activity data',
      similarCompounds: data.similarCompounds ? `Found ${data.similarCompounds.length} similar compounds` : undefined,
      cached: result.metadata?.cached ? 'Yes' : 'No',
    };
  }

  async cleanup(): Promise<void> {
    this.cache.clear();
    console.log('[ChEMBL Plugin] Cleaned up');
  }
}
