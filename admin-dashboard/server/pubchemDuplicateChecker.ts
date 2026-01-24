/**
 * PubChem Duplicate Checking
 * 
 * Prevents rediscovery of known compounds by checking PubChem database
 * Returns CID/SID and compound information if found
 * 
 * API Documentation: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
 * No authentication required - public API
 */

interface PubChemCompound {
  cid: number;
  iupac_name?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  canonical_smiles?: string;
  inchi_key?: string;
  synonyms?: string[];
  pubmed_count?: number;
  patent_count?: number;
  bioassay_count?: number;
}

interface DuplicateCheckResult {
  is_known_compound: boolean;
  pubchem_cid?: number;
  compound_info?: PubChemCompound;
  novelty_score: number; // 0-100, higher is more novel
  analysis: string;
}

/**
 * Check if a compound exists in PubChem database
 */
export async function checkPubChemDuplicate(smiles: string): Promise<DuplicateCheckResult> {
  try {
    // Step 1: Search PubChem by SMILES
    const compound = await searchPubChemBySmiles(smiles);
    
    if (!compound) {
      return {
        is_known_compound: false,
        novelty_score: 100,
        analysis: 'This compound is not found in PubChem database. It appears to be a novel chemical entity with no prior documentation.',
      };
    }

    // Step 2: Get detailed compound information
    const compoundInfo = await getPubChemCompoundInfo(compound.cid);
    
    // Step 3: Calculate novelty score based on literature/patent presence
    const noveltyScore = calculateNoveltyScore(compoundInfo);
    
    // Step 4: Generate analysis
    const analysis = generateDuplicateAnalysis(compoundInfo, noveltyScore);

    return {
      is_known_compound: true,
      pubchem_cid: compound.cid,
      compound_info: compoundInfo,
      novelty_score: noveltyScore,
      analysis,
    };
  } catch (error) {
    console.error('PubChem duplicate check failed:', error);
    return {
      is_known_compound: false,
      novelty_score: 50,
      analysis: 'Unable to verify compound novelty due to PubChem API error. Manual verification recommended.',
    };
  }
}

/**
 * Search PubChem by SMILES structure
 */
async function searchPubChemBySmiles(smiles: string): Promise<{ cid: number } | null> {
  try {
    // PubChem PUG REST API - Structure search
    const searchUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`;
    
    const response = await fetch(searchUrl, {
      headers: {
        'Accept': 'application/json',
      },
    });

    if (!response.ok) {
      if (response.status === 404) {
        // Compound not found - this is expected for novel compounds
        return null;
      }
      throw new Error(`PubChem API error: ${response.statusText}`);
    }

    const data = await response.json();
    
    if (data.IdentifierList?.CID && data.IdentifierList.CID.length > 0) {
      return { cid: data.IdentifierList.CID[0] };
    }
    
    return null;
  } catch (error) {
    console.error('PubChem SMILES search failed:', error);
    return null;
  }
}

/**
 * Get detailed compound information from PubChem
 */
async function getPubChemCompoundInfo(cid: number): Promise<PubChemCompound> {
  try {
    // Get compound properties
    const propsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON`;
    const propsResponse = await fetch(propsUrl);
    const propsData = await propsResponse.json();
    const props = propsData.PropertyTable?.Properties?.[0] || {};

    // Get synonyms
    const synonymsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/synonyms/JSON`;
    const synonymsResponse = await fetch(synonymsUrl);
    const synonymsData = await synonymsResponse.json();
    const synonyms = synonymsData.InformationList?.Information?.[0]?.Synonym || [];

    // Get literature/patent counts
    const countsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/xrefs/PubMedID,PatentID,BioAssayID/JSON`;
    const countsResponse = await fetch(countsUrl);
    const countsData = await countsResponse.json();
    const xrefs = countsData.InformationList?.Information?.[0] || {};

    return {
      cid,
      iupac_name: props.IUPACName,
      molecular_formula: props.MolecularFormula,
      molecular_weight: props.MolecularWeight,
      canonical_smiles: props.CanonicalSMILES,
      inchi_key: props.InChIKey,
      synonyms: synonyms.slice(0, 10), // Top 10 synonyms
      pubmed_count: xrefs.PubMedID?.length || 0,
      patent_count: xrefs.PatentID?.length || 0,
      bioassay_count: xrefs.BioAssayID?.length || 0,
    };
  } catch (error) {
    console.error('PubChem compound info retrieval failed:', error);
    return {
      cid,
    };
  }
}

/**
 * Calculate novelty score based on literature/patent presence
 * 100 = completely novel, 0 = extensively documented
 */
function calculateNoveltyScore(compound: PubChemCompound): number {
  let score = 100;

  // Penalize based on literature presence
  if (compound.pubmed_count) {
    score -= Math.min(40, compound.pubmed_count * 2);
  }

  // Penalize based on patent presence
  if (compound.patent_count) {
    score -= Math.min(30, compound.patent_count * 3);
  }

  // Penalize based on bioassay data
  if (compound.bioassay_count) {
    score -= Math.min(20, compound.bioassay_count);
  }

  // Penalize if it has common synonyms (indicates well-known compound)
  if (compound.synonyms && compound.synonyms.length > 5) {
    score -= 10;
  }

  return Math.max(0, Math.min(100, score));
}

/**
 * Generate duplicate analysis text
 */
function generateDuplicateAnalysis(compound: PubChemCompound, noveltyScore: number): string {
  const name = compound.iupac_name || compound.synonyms?.[0] || `CID ${compound.cid}`;
  
  let analysis = `Found in PubChem as ${name} (CID: ${compound.cid}). `;

  if (noveltyScore >= 80) {
    analysis += 'This compound is documented but has minimal research literature. It may still represent a novel therapeutic opportunity.';
  } else if (noveltyScore >= 50) {
    analysis += `Moderate documentation exists (${compound.pubmed_count || 0} PubMed articles, ${compound.patent_count || 0} patents). Review existing literature before proceeding.`;
  } else if (noveltyScore >= 20) {
    analysis += `Well-documented compound with ${compound.pubmed_count || 0} PubMed articles and ${compound.patent_count || 0} patents. Significant prior art exists.`;
  } else {
    analysis += `Extensively studied compound with ${compound.pubmed_count || 0} PubMed articles, ${compound.patent_count || 0} patents, and ${compound.bioassay_count || 0} bioassays. This is a known pharmaceutical entity.`;
  }

  return analysis;
}

/**
 * Batch check multiple compounds for duplicates
 */
export async function batchCheckPubChemDuplicates(
  smilesList: string[]
): Promise<Map<string, DuplicateCheckResult>> {
  const results = new Map<string, DuplicateCheckResult>();

  // Process in batches of 5 to avoid rate limiting
  for (let i = 0; i < smilesList.length; i += 5) {
    const batch = smilesList.slice(i, i + 5);
    const batchPromises = batch.map((smiles) => checkPubChemDuplicate(smiles));
    const batchResults = await Promise.all(batchPromises);

    batch.forEach((smiles, index) => {
      results.set(smiles, batchResults[index]);
    });

    // Add delay between batches to respect rate limits
    if (i + 5 < smilesList.length) {
      await new Promise((resolve) => setTimeout(resolve, 1000));
    }
  }

  return results;
}
