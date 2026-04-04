import { describe, it, expect } from 'vitest';
import { runMolecularDocking } from './molecularDockingWrapper';
import { analyzeToxicity, analyzeSyntheticAccessibility } from './advancedAnalysis';

describe('Python Path Fix - Docking and ADMET', () => {
  // Test ketamine analog SMILES
  const ketamineAnalogs = {
    ketamine: 'CN1C(=O)CC(c2ccccc2)C1=O',  // Original ketamine
    norketamine: 'NC(=O)CC(c1ccccc1)C(=O)N',  // N-demethylated
    dck: 'CN1C(=O)CC(c2ccc(Cl)cc2)C1=O',  // Deschloroketamine
    mxe: 'CCN1C(=O)CC(c2ccc(OC)cc2)C1=O',  // Methoxetamine
  };

  it('should run docking with ketamine on NMDA GluN2A', async () => {
    // Mock docking parameters for NMDA GluN2A
    const dockingParams = {
      smiles: ketamineAnalogs.ketamine,
      analogId: 'test-ketamine-1',
      targetName: 'NMDA',
      boxCenter: { x: 0, y: 0, z: 0 },
      boxSize: { x: 20, y: 20, z: 20 },
      exhaustiveness: 8,
      numPoses: 5,
    };

    try {
      const result = await runMolecularDocking(dockingParams);
      console.log('✅ Docking result:', result);
      
      expect(result).toBeDefined();
      expect(result.success).toBe(true); // Docking succeeded
      if (result.success && result.binding_affinity) {
        expect(result.binding_affinity).toBeLessThan(0); // Negative affinity indicates binding
        console.log(`✅ Ketamine binding affinity: ${result.binding_affinity} kcal/mol`);
      }
    } catch (error: any) {
      console.error('❌ Docking failed:', error.message);
      throw error;
    }
  }, { timeout: 120000 });

  it('should analyze toxicity of ketamine analog', async () => {
    try {
      const result = await analyzeToxicity(ketamineAnalogs.dck);
      console.log('✅ Toxicity analysis result:', result);
      
      expect(result).toBeDefined();
      expect(result.hERG).toBeDefined();
      expect(result.hepatotoxicity).toBeDefined();
      console.log(`✅ DCK hERG risk: ${result.hERG.risk_level}`);
    } catch (error: any) {
      console.error('❌ Toxicity analysis failed:', error.message);
      throw error;
    }
  }, { timeout: 60000 });

  it('should analyze synthetic accessibility of ketamine analog', async () => {
    try {
      const result = await analyzeSyntheticAccessibility(ketamineAnalogs.mxe);
      console.log('✅ SA Score result:', result);
      
      expect(result).toBeDefined();
      if (result.sa_score) expect(result.sa_score).toBeGreaterThan(0); // SA Score optional
      expect(result.difficulty).toBeDefined();
      console.log(`✅ MXE SA Score: ${result.sa_score} (${result.difficulty})`);
    } catch (error: any) {
      console.error('❌ SA Score analysis failed:', error.message);
      throw error;
    }
  }, { timeout: 60000 });

  it('should handle multiple ketamine analogs in sequence', async () => {
    const analogs = [
      { name: 'Ketamine', smiles: ketamineAnalogs.ketamine },
      { name: 'Norketamine', smiles: ketamineAnalogs.norketamine },
      { name: 'DCK', smiles: ketamineAnalogs.dck },
    ];

    for (const analog of analogs) {
      try {
        const toxResult = await analyzeToxicity(analog.smiles);
        console.log(`✅ ${analog.name} - hERG: ${toxResult.hERG.risk_level}`);
        expect(toxResult.hERG.risk_level).toBeDefined();
      } catch (error: any) {
        console.error(`❌ ${analog.name} analysis failed:`, error.message);
        throw error;
      }
    }
  }, { timeout: 180000 });
});
