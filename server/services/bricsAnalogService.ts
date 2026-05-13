/**
 * BRICS Analog Generation Service
 * Implements BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)
 * for generating novel analogs from parent compounds
 */

export interface BRICSAnalog {
  smiles: string;
  name: string;
  tanimotoSimilarity: number;
  molecularWeight: number;
  logP: number;
  hbdCount: number;
  hbaCount: number;
  rotBondCount: number;
  qed: number;
  pains: boolean;
  confidence: number;
}

/**
 * Generate BRICS analogs from a parent compound SMILES
 * Uses BRICS fragmentation and reassembly strategy
 */
export async function generateBRICSAnalogs(
  parentSmiles: string,
  numAnalogs: number = 25,
  similarityThreshold: number = 0.7
): Promise<BRICSAnalog[]> {
  try {
    // Import pythonBridge dynamically to avoid circular dependencies
    const { runPythonJson } = await import('../pythonBridge.js');
    
    console.log(`Generating ${numAnalogs} BRICS analogs from ${parentSmiles}`);
    console.log(`Similarity threshold: ${similarityThreshold}`);

    const result = await runPythonJson({
      modulePath: 'brics_generator.py',
      args: ['generate_brics_analogs', parentSmiles, numAnalogs, similarityThreshold],
      timeoutMs: 60000
    });

    if (!result.ok) {
      throw new Error(result.error?.message || 'Failed to generate BRICS analogs');
    }

    return result.data || [];
  } catch (error) {
    console.error("Error generating BRICS analogs:", error);
    throw error;
  }
}

/**
 * Generate analogs using substituent replacement strategy
 * Replaces functional groups with common drug-like substituents
 */
export async function generateSubstituentAnalogs(
  parentSmiles: string,
  numAnalogs: number = 25
): Promise<BRICSAnalog[]> {
  try {
    const { runPythonJson } = await import('../pythonBridge.js');
    
    console.log(`Generating ${numAnalogs} substituent analogs from ${parentSmiles}`);

    const result = await runPythonJson({
      modulePath: 'substituent_generator.py',
      args: ['generate_substituent_analogs', parentSmiles, numAnalogs],
      timeoutMs: 60000
    });

    if (!result.ok) {
      throw new Error(result.error?.message || 'Failed to generate substituent analogs');
    }

    return result.data || [];
  } catch (error) {
    console.error("Error generating substituent analogs:", error);
    throw error;
  }
}

/**
 * Generate analogs using scaffold hopping
 * Replaces core scaffold with bioisosteric alternatives
 */
export async function generateScaffoldHoppingAnalogs(
  parentSmiles: string,
  numAnalogs: number = 25
): Promise<BRICSAnalog[]> {
  try {
    const { runPythonJson } = await import('../pythonBridge.js');
    
    console.log(`Generating ${numAnalogs} scaffold-hopped analogs from ${parentSmiles}`);

    const result = await runPythonJson({
      modulePath: 'scaffold_hopper.py',
      args: ['generate_scaffold_hops', parentSmiles, numAnalogs],
      timeoutMs: 60000
    });

    if (!result.ok) {
      throw new Error(result.error?.message || 'Failed to generate scaffold-hopped analogs');
    }

    return result.data || [];
  } catch (error) {
    console.error("Error generating scaffold-hopped analogs:", error);
    throw error;
  }
}

/**
 * Score analog based on drug-likeness criteria
 * Uses Lipinski's Rule of Five, QED, and PAINS filters
 */
export function scoreAnalog(analog: Partial<BRICSAnalog>): number {
  let score = 100;

  // Lipinski's Rule of Five penalties
  if ((analog.molecularWeight || 0) > 500) score -= 10;
  if ((analog.logP || 0) > 5) score -= 10;
  if ((analog.hbdCount || 0) > 5) score -= 5;
  if ((analog.hbaCount || 0) > 10) score -= 5;

  // Rotatable bonds penalty (flexibility)
  if ((analog.rotBondCount || 0) > 10) score -= 5;

  // QED score (0-1, higher is better)
  if ((analog.qed || 0) < 0.5) score -= 15;

  // PAINS filter (if true, it's a PAINS compound - bad)
  if (analog.pains) score -= 20;

  return Math.max(0, Math.min(100, score));
}

/**
 * Filter analogs by drug-likeness criteria
 */
export function filterDrugLikeAnalogs(
  analogs: BRICSAnalog[],
  minQED: number = 0.5,
  maxMW: number = 500,
  excludePAINS: boolean = true
): BRICSAnalog[] {
  return analogs.filter((analog) => {
    if (analog.molecularWeight > maxMW) return false;
    if (analog.qed < minQED) return false;
    if (excludePAINS && analog.pains) return false;
    return true;
  });
}

/**
 * Rank analogs by multiple criteria
 */
export function rankAnalogs(
  analogs: BRICSAnalog[],
  weights: {
    similarity?: number;
    qed?: number;
    mw?: number;
    confidence?: number;
  } = {}
): BRICSAnalog[] {
  const w = {
    similarity: weights.similarity || 0.4,
    qed: weights.qed || 0.3,
    mw: weights.mw || 0.2,
    confidence: weights.confidence || 0.1,
  };

  const scored = analogs.map((analog) => {
    // Normalize scores to 0-1
    const simScore = analog.tanimotoSimilarity; // already 0-1
    const qedScore = analog.qed; // already 0-1
    const mwScore = Math.max(0, 1 - analog.molecularWeight / 500); // lower MW is better
    const confScore = analog.confidence / 100; // convert 0-100 to 0-1

    const totalScore =
      simScore * w.similarity +
      qedScore * w.qed +
      mwScore * w.mw +
      confScore * w.confidence;

    return { ...analog, score: totalScore };
  });

  return scored.sort((a, b) => (b as any).score - (a as any).score);
}
