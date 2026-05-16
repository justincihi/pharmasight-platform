import { executePythonScriptSafe } from './pythonBridgeSafe';

/**
 * Wrapper for cheminformatics Python workflows
 * Executes Python scripts with graceful fallback to mock responses
 */

export interface CheminformaticsResult {
  success: boolean;
  data?: any;
  error?: string;
}

/**
 * Execute a Python cheminformatics script with safe wrapper
 */
async function executePythonScript(
  scriptName: string,
  functionName: string,
  args: unknown[]
): Promise<CheminformaticsResult> {
  try {
    const result = await executePythonScriptSafe(scriptName, functionName, args);
    return {
      success: true,
      data: result,
    };
  } catch (error) {
    return {
      success: false,
      error: error instanceof Error ? error.message : 'Unknown error',
    };
  }
}

/**
 * Confirm SMILES and fetch similar compounds from PubChem
 */
export async function confirmAndFetchSimilars(
  nameOrSmiles: string,
  threshold: number = 0.70,
  maxHits: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'confirm_and_fetch_similars', [
    nameOrSmiles,
    threshold,
    maxHits,
  ]);
}

/**
 * Check patent status for a compound
 */
export async function checkPatentStatus(cid: number): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'check_patent_status', [cid]);
}

/**
 * Screen compounds and flag patent-free candidates
 */
export async function screenAndFlagForMasterlist(
  hits: any[],
  maxHits: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'screen_and_flag_for_masterlist', [
    hits,
    maxHits,
  ]);
}

/**
 * Generate BRICS analogs from a parent SMILES
 */
export async function generateBricsAnalogs(
  smiles: string,
  n: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'generate_brics_analogs', [
    smiles,
    n,
  ]);
}

/**
 * Enumerate substituent analogs with R-group substitution
 */
export async function enumerateSubstituentAnalogs(
  baseSmiles: string,
  attachmentIdx: number,
  n: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'enumerate_substituent_analogs', [
    baseSmiles,
    attachmentIdx,
    n,
  ]);
}

/**
 * Full integrated analog pipeline
 */
export async function fullAnalogPipeline(
  inputSmiles: string,
  threshold: number = 0.70,
  maxHits: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'full_analog_pipeline', [
    inputSmiles,
    threshold,
    maxHits,
  ]);
}

/**
 * Validate SMILES string
 */
export async function validateSmiles(smiles: string): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', 'validate_smiles', [smiles]);
}
