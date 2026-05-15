/**
 * PubChem API Wrapper
 * Provides functions for SMILES validation, similarity search, and patent screening
 * Using PubChem's free REST API with rate limiting and caching
 */

import { cache } from "./cache";

const PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
const RATE_LIMIT_DELAY = 300; // 300ms between requests to respect rate limits

// Simple in-memory rate limiter
let lastRequestTime = 0;

async function rateLimit() {
  const now = Date.now();
  const timeSinceLastRequest = now - lastRequestTime;
  if (timeSinceLastRequest < RATE_LIMIT_DELAY) {
    await new Promise((resolve) =>
      setTimeout(resolve, RATE_LIMIT_DELAY - timeSinceLastRequest)
    );
  }
  lastRequestTime = Date.now();
}

/**
 * Validate SMILES string using PubChem
 */
export async function validateSmiles(smiles: string): Promise<{
  valid: boolean;
  canonicalSmiles?: string;
  molecularFormula?: string;
  molecularWeight?: number;
  error?: string;
}> {
  try {
    const cacheKey = `smiles_valid_${smiles}`;
    const cached = cache.get<{ valid: boolean; canonicalSmiles?: string; molecularFormula?: string; molecularWeight?: number; error?: string }>(cacheKey);
    if (cached) return cached;

    await rateLimit();

    const response = await fetch(
      `${PUBCHEM_BASE}/compound/smiles/${encodeURIComponent(smiles)}/property/CanonicalSMILES,MolecularFormula,MolecularWeight/JSON`
    );

    if (!response.ok) {
      return {
        valid: false,
        error: `SMILES validation failed: ${response.statusText}`,
      } as const;
    }

    const data = (await response.json()) as {
      properties?: Array<{
        CanonicalSMILES?: string;
        MolecularFormula?: string;
        MolecularWeight?: number;
      }>;
    };

    if (!data.properties || data.properties.length === 0) {
      return { valid: false, error: "Invalid SMILES string" } as const;
    }

    const result = {
      valid: true,
      canonicalSmiles: data.properties[0].CanonicalSMILES,
      molecularFormula: data.properties[0].MolecularFormula,
      molecularWeight: data.properties[0].MolecularWeight,
    };

    cache.set(cacheKey, result, 86400); // Cache for 24 hours
    return result;
  } catch (error) {
    return {
      valid: false,
      error: `SMILES validation error: ${error instanceof Error ? error.message : String(error)}`,
    } as const;
  }
}

/**
 * Search for similar compounds using PubChem
 */
export async function searchSimilarCompounds(
  smiles: string,
  threshold: number = 0.9,
  maxResults: number = 10
): Promise<{
  success: boolean;
  compounds?: Array<{
    cid: number;
    smiles: string;
    name: string;
    similarity: number;
    molecularWeight: number;
  }>;
  error?: string;
}> {
  try {
    const cacheKey = `similar_${smiles}_${threshold}_${maxResults}`;
    const cached = cache.get<{ success: boolean; compounds?: Array<{ cid: number; smiles: string; name: string; similarity: number; molecularWeight: number }>; error?: string }>(cacheKey);
    if (cached) return cached;

    await rateLimit();

    // First, get CID for the query SMILES
    const cidResponse = await fetch(
      `${PUBCHEM_BASE}/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`
    );

    if (!cidResponse.ok) {
      return { success: false, error: "Could not find compound" } as const;
    }

    const cidData = (await cidResponse.json()) as { IdentifierList?: { CID?: number[] } };
    const queryCid = cidData.IdentifierList?.CID?.[0];

    if (!queryCid) {
      return { success: false, error: "Could not find compound CID" } as const;
    }

    // Search for similar compounds
    await rateLimit();
    const similarResponse = await fetch(
      `${PUBCHEM_BASE}/compound/cid/${queryCid}/cids/JSON?cids_type=similar&Threshold=${Math.round(threshold * 100)}&MaxResults=${maxResults}`
    );

    if (!similarResponse.ok) {
      return { success: false, error: "Similarity search failed" } as const;
    }

    const similarData = (await similarResponse.json()) as {
      IdentifierList?: { CID?: number[] };
    };
    const similarCids = similarData.IdentifierList?.CID || [];

    // Get properties for similar compounds
    const compounds = [];
    for (const cid of similarCids.slice(0, maxResults)) {
      await rateLimit();
      const propResponse = await fetch(
        `${PUBCHEM_BASE}/compound/cid/${cid}/property/CanonicalSMILES,IUPACName,MolecularWeight/JSON`
      );

      if (propResponse.ok) {
        const propData = (await propResponse.json()) as {
          properties?: Array<{
            CanonicalSMILES?: string;
            IUPACName?: string;
            MolecularWeight?: number;
          }>;
        };
        if (propData.properties && propData.properties.length > 0) {
          compounds.push({
            cid,
            smiles: propData.properties[0].CanonicalSMILES || "",
            name: propData.properties[0].IUPACName || `CID ${cid}`,
            similarity: threshold,
            molecularWeight: propData.properties[0].MolecularWeight || 0,
          });
        }
      }
    }

    const result = { success: true, compounds };
    cache.set(cacheKey, result, 86400);
    return result;
  } catch (error) {
    return {
      success: false,
      error: `Similarity search error: ${error instanceof Error ? error.message : String(error)}`,
    } as const;
  }
}

/**
 * Check patent status for a compound
 */
export async function checkPatentStatus(
  cid: number
): Promise<{
  success: boolean;
  hasPatents: boolean;
  patentCount?: number;
  error?: string;
}> {
  try {
    const cacheKey = `patent_${cid}`;
    const cached = cache.get<{ success: boolean; hasPatents: boolean; patentCount?: number; error?: string }>(cacheKey);
    if (cached) return cached;

    await rateLimit();

    const response = await fetch(
      `${PUBCHEM_BASE}/compound/cid/${cid}/data/JSON?data_type=patent`
    );

    if (!response.ok) {
      return { success: false, hasPatents: false, patentCount: 0 } as const;
    }

    const data = (await response.json()) as {
      InformationList?: { Information?: Array<{ PatentID?: string[] }> };
    };
    const patents = data.InformationList?.Information?.[0]?.PatentID || [];

    const result = {
      success: true,
      hasPatents: patents.length > 0,
      patentCount: patents.length,
    } as const;

    cache.set(cacheKey, result, 86400);
    return result;
  } catch (error) {
    return {
      success: false,
      hasPatents: false,
      error: `Patent check error: ${error instanceof Error ? error.message : String(error)}`,
    } as const;
  }
}

/**
 * Get compound properties from PubChem
 */
export async function getCompoundProperties(
  smiles: string
): Promise<{
  success: boolean;
  properties?: {
    cid: number;
    smiles: string;
    name: string;
    molecularWeight: number;
    molecularFormula: string;
    logP?: number;
    hbondDonor?: number;
    hbondAcceptor?: number;
  };
  error?: string;
}> {
  try {
    const cacheKey = `props_${smiles}`;
    const cached = cache.get<{ success: boolean; properties?: { cid: number; smiles: string; name: string; molecularWeight: number; molecularFormula: string; logP?: number; hbondDonor?: number; hbondAcceptor?: number }; error?: string }>(cacheKey);
    if (cached) return cached;

    await rateLimit();

    const response = await fetch(
      `${PUBCHEM_BASE}/compound/smiles/${encodeURIComponent(smiles)}/property/CanonicalSMILES,IUPACName,MolecularWeight,MolecularFormula,LogP,HBondDonorCount,HBondAcceptorCount/JSON`
    );

    if (!response.ok) {
      return { success: false, error: "Could not fetch properties" } as const;
    }

    const data = (await response.json()) as {
      properties?: Array<{
        CID?: number;
        CanonicalSMILES?: string;
        IUPACName?: string;
        MolecularWeight?: number;
        MolecularFormula?: string;
        LogP?: number;
        HBondDonorCount?: number;
        HBondAcceptorCount?: number;
      }>;
    };

    if (!data.properties || data.properties.length === 0) {
      return { success: false, error: "No properties found" } as const;
    }

    const prop = data.properties[0];
    const result = {
      success: true,
      properties: {
        cid: prop.CID || 0,
        smiles: prop.CanonicalSMILES || "",
        name: prop.IUPACName || "Unknown",
        molecularWeight: prop.MolecularWeight || 0,
        molecularFormula: prop.MolecularFormula || "",
        logP: prop.LogP,
        hbondDonor: prop.HBondDonorCount,
        hbondAcceptor: prop.HBondAcceptorCount,
      },
    };

    cache.set(cacheKey, result, 86400);
    return result;
  } catch (error) {
    return {
      success: false,
      error: `Property fetch error: ${error instanceof Error ? error.message : String(error)}`,
    } as const;
  }
}
