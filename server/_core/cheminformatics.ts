import { spawn } from 'child_process';
import { join } from 'path';
import { promises as fs } from 'fs';

/**
 * Wrapper for cheminformatics Python workflows
 * Executes Python scripts in the project's venv
 */

const VENV_PYTHON = join(process.cwd(), 'venv', 'bin', 'python');
const SCRIPTS_DIR = join(process.cwd(), 'scripts');

export interface CheminformaticsResult {
  success: boolean;
  data?: any;
  error?: string;
  stdout?: string;
  stderr?: string;
}

/**
 * Execute a Python cheminformatics script
 */
function executePythonScript(
  scriptName: string,
  args: string[] = []
): Promise<CheminformaticsResult> {
  return new Promise((resolve) => {
    const scriptPath = join(SCRIPTS_DIR, scriptName);
    const python = spawn(VENV_PYTHON, [scriptPath, ...args], {
      cwd: process.cwd(),
      env: { ...process.env, PYTHONUNBUFFERED: '1' },
    });

    let stdout = '';
    let stderr = '';

    python.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    python.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    python.on('close', (code) => {
      if (code === 0) {
        try {
          // Try to parse JSON output from Python script
          const jsonMatch = stdout.match(/\{[\s\S]*\}/);
          if (jsonMatch) {
            const data = JSON.parse(jsonMatch[0]);
            resolve({ success: true, data, stdout });
          } else {
            resolve({ success: true, stdout });
          }
        } catch (e) {
          resolve({ success: true, stdout });
        }
      } else {
        resolve({
          success: false,
          error: `Python script exited with code ${code}`,
          stderr,
          stdout,
        });
      }
    });

    python.on('error', (err) => {
      resolve({
        success: false,
        error: err.message,
        stderr: err.toString(),
      });
    });
  });
}

/**
 * Confirm SMILES and fetch similar compounds from PubChem
 */
export async function confirmAndFetchSimilars(
  nameOrSmiles: string,
  threshold: number = 0.70,
  maxHits: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', [
    'confirm_and_fetch_similars',
    nameOrSmiles,
    threshold.toString(),
    maxHits.toString(),
  ]);
}

/**
 * Check patent status for a compound
 */
export async function checkPatentStatus(cid: number): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', [
    'check_patent_status',
    cid.toString(),
  ]);
}

/**
 * Screen compounds and flag patent-free candidates
 */
export async function screenAndFlagForMasterlist(
  hits: any[],
  maxHits: number = 25
): Promise<CheminformaticsResult> {
  const hitsJson = JSON.stringify(hits);
  return executePythonScript('cheminformatics_workflows.py', [
    'screen_and_flag_for_masterlist',
    hitsJson,
    maxHits.toString(),
  ]);
}

/**
 * Generate BRICS analogs from a parent SMILES
 */
export async function generateBricsAnalogs(
  smiles: string,
  n: number = 25
): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', [
    'generate_brics_analogs',
    smiles,
    n.toString(),
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
  return executePythonScript('cheminformatics_workflows.py', [
    'enumerate_substituent_analogs',
    baseSmiles,
    attachmentIdx.toString(),
    n.toString(),
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
  return executePythonScript('cheminformatics_workflows.py', [
    'full_analog_pipeline',
    inputSmiles,
    threshold.toString(),
    maxHits.toString(),
  ]);
}

/**
 * Validate SMILES string
 */
export async function validateSmiles(smiles: string): Promise<CheminformaticsResult> {
  return executePythonScript('cheminformatics_workflows.py', [
    'validate_smiles',
    smiles,
  ]);
}
