import { describe, it, expect } from 'vitest';
import { runComprehensiveAnalysis } from './advancedAnalysis';
import { simulatePKPD, predictToxicity } from './pythonBridge';

describe('Python Integrations', () => {
  // Test compound: Ketamine SMILES
  const testSmiles = 'CNC1(c2cccc(F)c2Cl)CCCCC1=O';
  
  describe('ADMET Integration', () => {
    it('should run comprehensive ADMET analysis', async () => {
      const result = await runComprehensiveAnalysis(testSmiles);
      
      expect(result).toBeDefined();
      expect(result.smiles).toBe(testSmiles);
      expect(result.status).toBe('success');
      expect(result.toxicity_profile).toBeDefined();
      expect(result.synthetic_accessibility).toBeDefined();
    }, 60000); // 60s timeout for Python execution
  });
  
  describe('PK/PD Integration', () => {
    it('should simulate PK/PD for oral dose', async () => {
      const result = await simulatePKPD(testSmiles, 100, 'oral');
      
      expect(result).toBeDefined();
      expect(result.success).toBe(true);
      expect(result.data).toBeDefined();
    }, 60000);
  });
  
  describe('Toxicity Integration', () => {
    it('should predict toxicity profile', async () => {
      const result = await predictToxicity(testSmiles);
      
      expect(result).toBeDefined();
      expect(result.success).toBe(true);
      expect(result.data).toBeDefined();
      expect(result.data.smiles).toBe(testSmiles);
    }, 60000);
  });
});
