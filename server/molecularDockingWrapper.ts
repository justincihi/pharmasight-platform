import { spawn } from 'child_process';
import { join } from 'path';
import { fileURLToPath } from 'url';
import { dirname } from 'path';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

interface DockingResult {
  success: boolean;
  binding_affinity?: number;
  rmsd_lb?: number;
  rmsd_ub?: number;
  num_poses?: number;
  top_poses?: Array<{
    rank: number;
    affinity: number;
    rmsd_lb: number;
    rmsd_ub: number;
  }>;
  error?: string;
  note?: string;
}

/**
 * Run molecular docking for a compound against specified receptor target
 */
export async function runMolecularDocking(smiles: string, compoundId: string, target: string = 'NMDA'): Promise<DockingResult> {
  return new Promise((resolve, reject) => {
    // Use the existing Python module which has mock docking
    const pythonScript = join(__dirname, 'python_modules', 'molecular_docking.py');
    
    // For now, use mock docking results since the Python module already has this functionality
    // In production, you would call the actual docking script
    const python = spawn('python3', ['-c', `
import sys
import json
sys.path.insert(0, '${join(__dirname, 'python_modules')}')

from molecular_docking import get_molecular_docking

docking = get_molecular_docking()
result = docking.dock_ligand('${smiles.replace(/'/g, "\\'")}')
print(json.dumps(result))
`]);

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
        console.error('[Molecular Docking] Python error:', stderr);
        reject(new Error(stderr || 'Docking failed'));
        return;
      }

      try {
        const lines = stdout.trim().split('\n');
        const lastLine = lines[lines.length - 1];
        const result = JSON.parse(lastLine);
        resolve(result);
      } catch (error: any) {
        reject(new Error(`Failed to parse docking output: ${error.message}`));
      }
    });
  });
}
