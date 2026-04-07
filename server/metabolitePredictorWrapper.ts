import { spawn } from 'child_process';
import path from 'path';

export interface Metabolite {
  smiles: string;
  parent_smiles: string;
  transformation: string;
  phase: 'Phase I' | 'Phase II';
  enzyme: string;
  probability: number;
  molecular_weight: number;
  logp: number;
}

export interface MetabolicStability {
  stability_score: number;
  classification: 'High' | 'Moderate' | 'Low';
  num_metabolites: number;
  avg_probability: number;
  analysis: string;
}

export interface MetabolitePredictionResult {
  parent_smiles: string;
  metabolic_stability: MetabolicStability;
  metabolites: Metabolite[];
}

/**
 * Predict metabolites for a given SMILES string using Python RDKit module
 */
export async function predictMetabolites(
  smiles: string,
  maxMetabolites: number = 10
): Promise<MetabolitePredictionResult> {
  return new Promise((resolve, reject) => {
    const pythonScript = path.join(__dirname, 'python_modules', 'metabolite_predictor.py');
    const python = spawn('python3', [pythonScript, smiles]);

    let stdout = '';
    let stderr = '';

    python.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    python.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    python.on('close', (code) => {
      if (code !== 0) {
        reject(new Error(`Metabolite prediction failed: ${stderr}`));
        return;
      }

      try {
        const result = JSON.parse(stdout);
        resolve(result);
      } catch (error) {
        reject(new Error(`Failed to parse metabolite prediction result: ${error}`));
      }
    });

    python.on('error', (error) => {
      reject(new Error(`Failed to spawn Python process: ${error.message}`));
    });
  });
}

/**
 * Store metabolites in database
 */
export async function storeMetabolites(
  analogId: number,
  metabolites: Metabolite[]
): Promise<void> {
  const { getDb } = await import('./db');
  const db = await getDb();
  const { metabolites: metabolitesTable } = await import('../drizzle/schema');

  const insertData = metabolites.map((met) => ({
    parentAnalogId: analogId,
    smiles: met.smiles,
    transformation: met.transformation,
    phase: met.phase,
    enzyme: met.enzyme,
    probability: met.probability.toString(),
    molecularWeight: met.molecular_weight.toString(),
    logP: met.logp.toString(),
  }));

  if (!db) throw new Error('Database connection failed');
  await db.insert(metabolitesTable).values(insertData);
}

/**
 * Get metabolites for an analog from database
 */
export async function getMetabolitesForAnalog(analogId: number) {
  const { getDb } = await import('./db');
  const db = await getDb();
  const { metabolites } = await import('../drizzle/schema');
  const { eq } = await import('drizzle-orm');

  if (!db) throw new Error('Database connection failed');
  const results = await db
    .select()
    .from(metabolites)
    .where(eq(metabolites.parentAnalogId, analogId))
    .orderBy(metabolites.probability);

  return results;
}
