import { spawn } from "child_process";
import { join, dirname } from "path";
import { fileURLToPath } from "url";

/**
 * Python Bridge for Cheminformatics
 * Executes Python scripts from the python_modules directory
 */

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const PYTHON_MODULES_PATH = "/home/ubuntu/pharmasight-admin-dashboard/server/python_modules";
const VENV_PYTHON = "/home/ubuntu/pharmasight-admin-dashboard/server/python_modules/venv/bin/python3";

interface PythonResult {
  success: boolean;
  data?: any;
  error?: string;
}

/**
 * Execute a Python script and return the result
 */
export async function executePythonScript(
  scriptName: string,
  functionName: string,
  args: any[] = []
): Promise<PythonResult> {
  return new Promise((resolve) => {
    const scriptPath = join(PYTHON_MODULES_PATH, scriptName);
    
    // Create a Python wrapper script that imports and calls the function
    const pythonCode = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    from ${scriptName.replace('.py', '')} import ${functionName}
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else []
    result = ${functionName}(*args)
    print(json.dumps({"success": True, "data": result}))
except Exception as e:
    print(json.dumps({"success": False, "error": str(e)}))
`;

    const python = spawn(VENV_PYTHON, ["-c", pythonCode, JSON.stringify(args)]);
    
    let stdout = "";
    let stderr = "";

    python.stdout.on("data", (data) => {
      stdout += data.toString();
    });

    python.stderr.on("data", (data) => {
      stderr += data.toString();
    });

    python.on("close", (code) => {
      if (code !== 0) {
        resolve({
          success: false,
          error: stderr || `Python process exited with code ${code}`,
        });
        return;
      }

      try {
        const result = JSON.parse(stdout);
        resolve(result);
      } catch (error) {
        resolve({
          success: false,
          error: `Failed to parse Python output: ${stdout}`,
        });
      }
    });
  });
}

/**
 * ChEMBL Validation
 */
export async function validateWithChEMBL(smiles: string): Promise<PythonResult> {
  return executePythonScript("chembl_validation.py", "validate_compound", [smiles]);
}

/**
 * ADMET Prediction
 */
export async function predictADMET(smiles: string): Promise<PythonResult> {
  return executePythonScript("admet_predictor_advanced.py", "predict_admet", [smiles]);
}

/**
 * Molecular Docking
 */
export async function runMolecularDocking(
  ligandSmiles: string,
  receptorPDB: string
): Promise<PythonResult> {
  return executePythonScript("molecular_docking.py", "dock_ligand", [
    ligandSmiles,
    receptorPDB,
  ]);
}

/**
 * Toxicity Prediction
 */
export async function predictToxicity(smiles: string): Promise<PythonResult> {
  return executePythonScript("toxicity_prediction.py", "predict_toxicity", [smiles]);
}

/**
 * PK/PD Simulation
 */
export async function simulatePKPD(
  smiles: string,
  dose: number,
  route: string
): Promise<PythonResult> {
  return executePythonScript("pkpd_pbpk_simulator.py", "simulate_pkpd", [
    smiles,
    dose,
    route,
  ]);
}

/**
 * Generate Analogs using RDKit
 */
export async function generateAnalogs(
  parentSmiles: string,
  numAnalogs: number = 10
): Promise<PythonResult> {
  return executePythonScript("rdkit_analog_generator.py", "generate_analogs", [
    parentSmiles,
    numAnalogs,
  ]);
}

/**
 * Query External Databases (PubChem, ChEMBL, etc.)
 */
export async function queryExternalDatabase(
  database: string,
  query: string
): Promise<PythonResult> {
  return executePythonScript("external_database_apis.py", "query_database", [
    database,
    query,
  ]);
}
