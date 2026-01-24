import { spawn } from 'child_process';
import path from 'path';

const PYTHON_VENV = path.join(__dirname, 'python_modules', 'venv', 'bin', 'python');
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
    const python = spawn(PYTHON_VENV, [ANALYSIS_SCRIPT, command, smiles, ...args]);
    
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
    
    python.on('error', (err) => {
      reject(new Error(`Failed to spawn Python process: ${err.message}`));
    });
  });
}

export async function runComprehensiveAnalysis(smiles: string): Promise<ComprehensiveAnalysisResult> {
  return runPythonScript('comprehensive', smiles);
}

export async function runToxicityAnalysis(smiles: string): Promise<ToxicityProfile> {
  return runPythonScript('toxicity', smiles);
}

export async function runSAAnalysis(smiles: string): Promise<SAScore> {
  return runPythonScript('sa_score', smiles);
}

export async function runOptimizationAnalysis(
  smiles: string,
  targetProperty?: string
): Promise<OptimizationSuggestion[]> {
  if (targetProperty) {
    return runPythonScript('optimize', smiles, targetProperty);
  }
  return runPythonScript('optimize', smiles);
}

export async function runIterativeOptimization(
  smiles: string,
  targetProperty: string,
  iterations: number = 3
): Promise<any[]> {
  return runPythonScript('iterative_optimize', smiles, targetProperty, iterations.toString());
}

// Helper to check if Python environment is available
export async function checkPythonEnvironment(): Promise<boolean> {
  try {
    const result = await runPythonScript('sa_score', 'CCO'); // Test with ethanol
    return !result.error;
  } catch (e) {
    return false;
  }
}
