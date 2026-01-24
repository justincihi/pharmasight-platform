import { spawn } from "child_process";
import { join, dirname } from "path";
import { fileURLToPath } from "url";

/**
 * Python Bridge for Cheminformatics
 * Executes Python scripts from the python_modules directory
 */

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Use relative path or environment variable for Python modules
const PYTHON_MODULES_PATH = process.env.PYTHON_MODULES_PATH || join(__dirname, "python_modules");

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

    const python = spawn("python3", ["-c", pythonCode, JSON.stringify(args)]);
    
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

/**
 * BioTransformer - Metabolite Prediction
 */
export async function predictMetabolites(
  smiles: string,
  metabolismType: string = 'human',
  steps: number = 1
): Promise<PythonResult> {
  return new Promise((resolve) => {
    const pythonCode = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    from biotransformer_wrapper import BioTransformerPredictor
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else []

    predictor = BioTransformerPredictor()
    result = predictor.predict_metabolites(
        smiles=args[0],
        metabolism_type=args[1] if len(args) > 1 else 'human',
        steps=args[2] if len(args) > 2 else 1
    )
    print(json.dumps({"success": True, "data": result}))
except Exception as e:
    import traceback
    print(json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}))
`;

    const python = spawn("python3", ["-c", pythonCode, JSON.stringify([smiles, metabolismType, steps])]);

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
 * BioTransformer - Batch Prediction
 */
export async function batchPredictMetabolites(
  compounds: Array<{ smiles: string; id?: string }>,
  metabolismType: string = 'human',
  steps: number = 1
): Promise<PythonResult> {
  return new Promise((resolve) => {
    const pythonCode = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    from biotransformer_wrapper import BioTransformerPredictor
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else []

    predictor = BioTransformerPredictor()
    result = predictor.batch_predict(
        compounds=args[0],
        metabolism_type=args[1] if len(args) > 1 else 'human',
        steps=args[2] if len(args) > 2 else 1
    )
    print(json.dumps({"success": True, "data": result}))
except Exception as e:
    import traceback
    print(json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}))
`;

    const python = spawn("python3", ["-c", pythonCode, JSON.stringify([compounds, metabolismType, steps])]);

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
 * Get available metabolism types from BioTransformer
 */
export async function getMetabolismTypes(): Promise<PythonResult> {
  return new Promise((resolve) => {
    const pythonCode = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    from biotransformer_wrapper import BioTransformerPredictor
    metabolism_types = BioTransformerPredictor.METABOLISM_TYPES
    print(json.dumps({"success": True, "data": metabolism_types}))
except Exception as e:
    print(json.dumps({"success": False, "error": str(e)}))
`;

    const python = spawn("python3", ["-c", pythonCode]);

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
