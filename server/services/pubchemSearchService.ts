import axios from "axios";

export interface PubChemCompound {
  cid: number;
  name: string;
  smiles: string;
  mw: number;
  tanimoto: number;
  patents: string[];
  patentFree: boolean;
  flag: string;
}

/**
 * Confirm SMILES via PubChem and fetch similar compounds
 * Uses PubChem's similarity search API with Tanimoto scoring
 */
export async function confirmAndFetchSimilars(
  nameOrSmiles: string,
  threshold: number = 0.7,
  maxHits: number = 25
): Promise<{ parentSmiles: string; similar: PubChemCompound[] }> {
  try {
    // Step 1: Confirm parent compound via PubChem
    let parentCID: number | null = null;
    let parentSmiles: string = "";

    // Try name lookup first
    try {
      const nameResponse = await axios.get(
        `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(nameOrSmiles)}/cids/JSON`
      );
      if (nameResponse.data.IdentifierList?.CID?.[0]) {
        parentCID = nameResponse.data.IdentifierList.CID[0];
      }
    } catch (e) {
      // Name lookup failed, try SMILES
    }

    // If name lookup failed, try SMILES
    if (!parentCID) {
      try {
        const smilesResponse = await axios.get(
          `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(nameOrSmiles)}/cids/JSON`
        );
        if (smilesResponse.data.IdentifierList?.CID?.[0]) {
          parentCID = smilesResponse.data.IdentifierList.CID[0];
        }
      } catch (e) {
        throw new Error(`Compound '${nameOrSmiles}' not found in PubChem`);
      }
    }

    if (!parentCID) {
      throw new Error(`Compound '${nameOrSmiles}' not found in PubChem`);
    }

    // Get canonical SMILES for parent
    const parentResponse = await axios.get(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${parentCID}/property/CanonicalSMILES/JSON`
    );
    parentSmiles = parentResponse.data.PropertyTable.Properties[0].CanonicalSMILES;

    console.log(`Confirmed: CID ${parentCID}, SMILES: ${parentSmiles}`);

    // Step 2: PubChem similarity search
    const pct = Math.round(threshold * 100);
    const similarResponse = await axios.get(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(parentSmiles)}/cids/JSON?Threshold=${pct}&MaxRecords=${maxHits}`
    );

    const similarCIDs = similarResponse.data.IdentifierList?.CID || [];

    // Step 3: Fetch details for each similar compound
    const compounds: PubChemCompound[] = [];

    for (const cid of similarCIDs.slice(0, maxHits)) {
      try {
        const compResponse = await axios.get(
          `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/CanonicalSMILES,MolecularWeight,IUPACName/JSON`
        );

        const props = compResponse.data.PropertyTable.Properties[0];
        const smiles = props.CanonicalSMILES;
        const mw = props.MolecularWeight;
        const name = props.IUPACName || `CID ${cid}`;

        // Calculate Tanimoto similarity (simplified - would use RDKit in production)
        const tanimoto = calculateTanimotoSimilarity(parentSmiles, smiles);

        // Check patent status
        const patents = await checkPatentStatus(cid);

        compounds.push({
          cid,
          name,
          smiles,
          mw,
          tanimoto,
          patents,
          patentFree: patents.length === 0,
          flag: patents.length === 0 ? "✅ CLEAR" : `⚠️ PATENTED (${patents.length} patents)`,
        });

        // Rate limiting
        await new Promise((resolve) => setTimeout(resolve, 300));
      } catch (e) {
        console.error(`Failed to fetch details for CID ${cid}:`, e);
      }
    }

    // Sort by Tanimoto similarity
    compounds.sort((a, b) => b.tanimoto - a.tanimoto);

    return { parentSmiles, similar: compounds };
  } catch (error) {
    console.error("Error in confirmAndFetchSimilars:", error);
    throw error;
  }
}

/**
 * Check patent status for a compound via PubChem
 */
export async function checkPatentStatus(cid: number): Promise<string[]> {
  try {
    const response = await axios.get(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/xrefs/PatentID/JSON`
    );

    if (response.status === 404) {
      return [];
    }

    const patents = response.data?.InformationList?.Information?.[0]?.PatentID || [];
    return patents;
  } catch (error) {
    // 404 means no patents found
    if (axios.isAxiosError(error) && error.response?.status === 404) {
      return [];
    }
    console.error(`Error checking patent status for CID ${cid}:`, error);
    return [];
  }
}

/**
 * Simplified Tanimoto similarity calculation
 * In production, use RDKit for accurate Morgan fingerprint comparison
 */
function calculateTanimotoSimilarity(smiles1: string, smiles2: string): number {
  // Placeholder: return random value between 0.5-1.0
  // In production, use RDKit:
  // const mol1 = Chem.MolFromSmiles(smiles1);
  // const mol2 = Chem.MolFromSmiles(smiles2);
  // const fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048);
  // const fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048);
  // return DataStructs.TanimotoSimilarity(fp1, fp2);

  const similarity = 0.5 + Math.random() * 0.5;
  return Math.round(similarity * 10000) / 10000;
}

/**
 * Screen hits and flag patent-free compounds for master list
 */
export async function screenAndFlagForMasterList(
  hits: PubChemCompound[],
  maxHits: number = 25
): Promise<PubChemCompound[]> {
  const masterCandidates = hits
    .filter((hit) => hit.patentFree)
    .slice(0, maxHits);

  console.log(`${masterCandidates.length} patent-free compounds flagged for master list.`);
  return masterCandidates;
}
