import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { spawn } from 'child_process';
import { join } from 'path';

/**
 * Test suite for cheminformatics workflows
 * Tests Python integration and tRPC procedures
 */

const VENV_PYTHON = join(process.cwd(), 'venv', 'bin', 'python');
const SCRIPTS_DIR = join(process.cwd(), 'scripts');

function executePythonScript(
  scriptName: string,
  args: string[] = []
): Promise<any> {
  return new Promise((resolve, reject) => {
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
          const jsonMatch = stdout.match(/\{[\s\S]*\}/);
          if (jsonMatch) {
            resolve(JSON.parse(jsonMatch[0]));
          } else {
            resolve({ stdout });
          }
        } catch (e) {
          reject(new Error(`Failed to parse JSON: ${stdout}`));
        }
      } else {
        reject(new Error(`Script failed: ${stderr}`));
      }
    });

    python.on('error', (err) => {
      reject(err);
    });
  });
}

describe('Cheminformatics Workflows', () => {
  describe('SMILES Validation', () => {
    it('should validate a valid SMILES string', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'COc1ccc2[nH]cc(CCN(C)C)c2c1',
      ]);

      expect(result.success).toBe(true);
      expect(result.valid).toBe(true);
      expect(result.canonical_smiles).toBeDefined();
      expect(result.mw).toBeGreaterThan(0);
      expect(result.num_atoms).toBeGreaterThan(0);
    });

    it('should reject an invalid SMILES string', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'INVALID_SMILES_XYZ',
      ]);

      expect(result.success).toBe(false);
      expect(result.valid).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should canonicalize SMILES correctly', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'c1ccccc1',
      ]);

      expect(result.success).toBe(true);
      expect(result.canonical_smiles).toBe('c1ccccc1');
      expect(result.num_atoms).toBe(6);
    });
  });

  describe('BRICS Analog Generation', () => {
    it('should generate BRICS analogs from a valid SMILES', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'generate_brics_analogs',
        'COc1ccc2[nH]cc(CCN(C)C)c2c1',
        '10',
      ]);

      expect(result.success).toBe(true);
      expect(result.parent_smiles).toBeDefined();
      expect(result.fragments).toBeDefined();
      expect(Array.isArray(result.analogs)).toBe(true);
      expect(result.total_generated).toBeGreaterThanOrEqual(0);
    });

    it('should handle invalid SMILES in BRICS generation', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'generate_brics_analogs',
        'INVALID',
        '10',
      ]);

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should respect max hits parameter', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'generate_brics_analogs',
        'c1ccccc1',
        '5',
      ]);

      expect(result.success).toBe(true);
      if (result.analogs) {
        expect(result.analogs.length).toBeLessThanOrEqual(5);
      }
    });

    it('should calculate Tanimoto similarity for analogs', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'generate_brics_analogs',
        'COc1ccc2[nH]cc(CCN(C)C)c2c1',
        '5',
      ]);

      expect(result.success).toBe(true);
      if (result.analogs && result.analogs.length > 0) {
        result.analogs.forEach((analog: any) => {
          expect(analog.tanimoto).toBeGreaterThanOrEqual(0);
          expect(analog.tanimoto).toBeLessThanOrEqual(1);
          expect(analog.smiles).toBeDefined();
          expect(analog.mw).toBeGreaterThan(0);
        });
      }
    });
  });

  describe('Patent Status Checking', () => {
    it('should check patent status for a valid CID', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'check_patent_status',
        '5289',
      ]);

      expect(result.success).toBe(true);
      expect(result.cid).toBe(5289);
      expect(Array.isArray(result.patents)).toBe(true);
      expect(result.patent_free).toBeDefined();
    });

    it('should return empty patents for non-existent CID', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'check_patent_status',
        '999999999',
      ]);

      expect(result.success).toBe(true);
      expect(result.patents).toBeDefined();
    });
  });

  describe('Similarity Search', () => {
    it('should search for similar compounds by SMILES', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'confirm_and_fetch_similars',
        'c1ccccc1',
        '0.70',
        '5',
      ]);

      expect(result.success).toBe(true);
      expect(result.parent_smiles).toBeDefined();
      expect(result.parent_cid).toBeDefined();
      expect(Array.isArray(result.hits)).toBe(true);
      expect(result.total_hits).toBeGreaterThanOrEqual(0);
    });

    it('should search for similar compounds by name', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'confirm_and_fetch_similars',
        'benzene',
        '0.70',
        '5',
      ]);

      expect(result.success).toBe(true);
      expect(result.parent_smiles).toBeDefined();
    });

    it('should respect similarity threshold', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'confirm_and_fetch_similars',
        'c1ccccc1',
        '0.90',
        '10',
      ]);

      expect(result.success).toBe(true);
      if (result.hits && result.hits.length > 0) {
        result.hits.forEach((hit: any) => {
          expect(hit.tanimoto).toBeGreaterThanOrEqual(0.90);
        });
      }
    });

    it('should handle invalid compound names', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'confirm_and_fetch_similars',
        'NONEXISTENT_COMPOUND_XYZ_123',
        '0.70',
        '5',
      ]);

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should return hits with patent information', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'confirm_and_fetch_similars',
        'c1ccccc1',
        '0.70',
        '3',
      ]);

      expect(result.success).toBe(true);
      if (result.hits && result.hits.length > 0) {
        result.hits.forEach((hit: any) => {
          expect(hit.cid).toBeDefined();
          expect(hit.smiles).toBeDefined();
          expect(hit.mw).toBeGreaterThan(0);
          expect(hit.tanimoto).toBeGreaterThanOrEqual(0);
          expect(hit.tanimoto).toBeLessThanOrEqual(1);
        });
      }
    });
  });

  describe('Full Pipeline', () => {
    it('should execute full analog pipeline', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'full_analog_pipeline',
        'COc1ccc2[nH]cc(CCN(C)C)c2c1',
        '0.70',
        '10',
      ]);

      expect(result.success).toBe(true);
      expect(result.input_smiles).toBeDefined();
      expect(result.canonical_smiles).toBeDefined();
      expect(result.in_pubchem).toBeDefined();
      expect(Array.isArray(result.master_list)).toBe(true);
      expect(result.total_candidates).toBeGreaterThanOrEqual(0);
      expect(result.patent_free_count).toBeGreaterThanOrEqual(0);
    });

    it('should flag patent-free candidates correctly', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'full_analog_pipeline',
        'c1ccccc1',
        '0.70',
        '5',
      ]);

      expect(result.success).toBe(true);
      if (result.master_list && result.master_list.length > 0) {
        result.master_list.forEach((candidate: any) => {
          expect(candidate.patent_free).toBe(true);
          expect(candidate.flag).toContain('CLEAR');
        });
      }
    });

    it('should handle invalid SMILES in full pipeline', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'full_analog_pipeline',
        'INVALID_SMILES',
        '0.70',
        '10',
      ]);

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should respect threshold parameter in full pipeline', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'full_analog_pipeline',
        'c1ccccc1',
        '0.95',
        '10',
      ]);

      expect(result.success).toBe(true);
      if (result.master_list && result.master_list.length > 0) {
        result.master_list.forEach((candidate: any) => {
          expect(candidate.tanimoto).toBeGreaterThanOrEqual(0.95);
        });
      }
    });
  });

  describe('Error Handling', () => {
    it('should handle missing workflow argument', async () => {
      try {
        await executePythonScript('cheminformatics_workflows.py', []);
        expect.fail('Should have thrown an error');
      } catch (error) {
        expect(error).toBeDefined();
      }
    });

    it('should handle unknown workflow', async () => {
      try {
        await executePythonScript('cheminformatics_workflows.py', [
          'unknown_workflow',
          'arg1',
        ]);
        expect.fail('Should have thrown an error');
      } catch (error) {
        expect(error).toBeDefined();
      }
    });

    it('should handle missing required arguments', async () => {
      try {
        await executePythonScript('cheminformatics_workflows.py', [
          'validate_smiles',
        ]);
        expect.fail('Should have thrown an error');
      } catch (error) {
        expect(error).toBeDefined();
      }
    });
  });

  describe('Data Integrity', () => {
    it('should return consistent SMILES canonicalization', async () => {
      const result1 = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'c1ccccc1',
      ]);

      const result2 = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'C1=CC=CC=C1',
      ]);

      expect(result1.canonical_smiles).toBe(result2.canonical_smiles);
    });

    it('should calculate consistent molecular weights', async () => {
      const result1 = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'c1ccccc1',
      ]);

      const result2 = await executePythonScript('cheminformatics_workflows.py', [
        'validate_smiles',
        'c1ccccc1',
      ]);

      expect(result1.mw).toBe(result2.mw);
    });

    it('should return sorted analogs by similarity', async () => {
      const result = await executePythonScript('cheminformatics_workflows.py', [
        'generate_brics_analogs',
        'COc1ccc2[nH]cc(CCN(C)C)c2c1',
        '10',
      ]);

      expect(result.success).toBe(true);
      if (result.analogs && result.analogs.length > 1) {
        for (let i = 0; i < result.analogs.length - 1; i++) {
          expect(result.analogs[i].tanimoto).toBeGreaterThanOrEqual(
            result.analogs[i + 1].tanimoto
          );
        }
      }
    });
  });
});
