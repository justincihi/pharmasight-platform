import { spawn, execSync } from 'child_process';
import { existsSync } from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Use the project venv Python which has rdkit, admet-ai, chemprop installed.
// System python3 does NOT have rdkit — always prefer the venv.
const VENV_PYTHON_PATH = path.join(__dirname, 'python_modules', 'venv', 'bin', 'python');
let PYTHON_EXECUTABLE = existsSync(VENV_PYTHON_PATH)
  ? VENV_PYTHON_PATH
  : '/usr/bin/python3';
console.log('[advancedAnalysis] Python executable:', PYTHON_EXECUTABLE, '| venv found:', existsSync(VENV_PYTHON_PATH));

const ANALYSIS_SCRIPT = path.join(__dirname, 'python_modules', 'comprehensive_analysis.py');

interface ToxicityProfile {
  smiles: string;
  hERG: {
    risk_score: number;
    risk_level: string;
    recommendation: string;
    logP: number;
    molecular_weight: number;
    tpsa: number;
    aromatic_rings: number;
    basic_groups: number;
    structural_alerts: number;
  };
  hepatotoxicity: {
    risk_score: number;
    risk_level: string;
    recommendation: string;
    logP: number;
    molecular_weight: number;
    structural_alerts: number;
  };
  mutagenicity: {
    risk_score: number;
    risk_level: string;
    prediction: string;
    recommendation: string;
    structural_alerts: number;
    alert_types: string[];
  };
  carcinogenicity: {
    risk_score: number;
    risk_level: string;
    prediction: string;
    recommendation: string;
    structural_alerts: number;
    alert_types: string[];
    aromatic_rings: number;
  };
}

interface SAScore {
  sa_score: number;
  difficulty: string;
  estimated_steps: string;
  recommendation: string;
  num_atoms: number;
  num_rings: number;
  num_stereocenters: number;
  num_rotatable_bonds: number;
  complexity_penalty: number;
  fragment_score: number;
}

interface OptimizationSuggestion {
  modification: string;
  category: string;
  rationale: string;
  original_smiles: string;
  optimized_smiles: string;
  property_changes: {
    mw_change: number;
    logP_change: number;
    tpsa_change: number;
  };
  original_properties: any;
  new_properties: any;
  priority: number;
}

interface ComprehensiveAnalysisResult {
  smiles: string;
  toxicity_profile: ToxicityProfile;
  synthetic_accessibility: SAScore;
  optimization_suggestions: OptimizationSuggestion[];
  status: string;
  error?: string;
}

function runPythonScript(command: string, smiles: string, ...args: string[]): Promise<any> {
  return new Promise((resolve, reject) => {
    const python = spawn(PYTHON_EXECUTABLE, [ANALYSIS_SCRIPT, command, smiles, ...args], { 
      env: { ...process.env, PYTHONUNBUFFERED: '1' },
      shell: false,
    });
    
    let stdout = '';
    let stderr = '';
    
    python.stdout.on('data', (data: Buffer) => {
      stdout += data.toString();
    });
    
    python.stderr.on('data', (data: Buffer) => {
      stderr += data.toString();
    });
    
    python.on('close', (code: number | null) => {
      if (code !== 0) {
        reject(new Error(`Python script failed: ${stderr}`));
        return;
      }
      
      try {
        const result = JSON.parse(stdout);
        resolve(result);
      } catch (e) {
        reject(new Error(`Failed to parse Python output: ${stdout}`));
      }
    });
    
    python.on('error', (err: Error) => {
      reject(new Error(`Failed to spawn Python process: ${err.message}`));
    });
  });
}

export async function analyzeToxicity(smiles: string): Promise<ToxicityProfile> {
  return runPythonScript('toxicity', smiles);
}

export async function analyzeSyntheticAccessibility(smiles: string): Promise<SAScore> {
  return runPythonScript('sascore', smiles);
}

export async function generateOptimizationSuggestions(smiles: string): Promise<OptimizationSuggestion[]> {
  return runPythonScript('optimize', smiles);
}

export async function runComprehensiveAnalysis(smiles: string): Promise<ComprehensiveAnalysisResult> {
  return runPythonScript('comprehensive', smiles);
}


export async function runSAAnalysis(smiles: string): Promise<SAScore> {
  return analyzeSyntheticAccessibility(smiles);
}

export async function runOptimizationAnalysis(smiles: string, targetProperty?: string): Promise<OptimizationSuggestion[]> {
  return generateOptimizationSuggestions(smiles);
}

export async function runIterativeOptimization(smiles: string, targetProperty: string, iterations: number): Promise<any> {
  // Run comprehensive analysis multiple times for iterative optimization
  const results = [];
  for (let i = 0; i < iterations; i++) {
    const result = await runComprehensiveAnalysis(smiles);
    results.push(result);
  }
  return { iterations: results, final_smiles: smiles };
}

export async function checkPythonEnvironment(): Promise<boolean> {
  try {
    // Test if Python can be executed and has required modules
    const result = execSync(`${PYTHON_EXECUTABLE} -c "import rdkit, meeko, vina; print('OK')"`, { 
      encoding: 'utf-8',
      timeout: 5000,
    });
    return result.includes('OK');
  } catch (e) {
    console.error('Python environment check failed:', e);
    return false;
  }
}


export async function runToxicityAnalysis(smiles: string): Promise<ToxicityProfile> {
  return analyzeToxicity(smiles);
}
