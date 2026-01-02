import { describe, it, expect, beforeAll } from 'vitest';
import { importFromSDF } from './sdfImporter';
import { exportAnalogs } from './batchExporter';
import { getDb } from './db';
import { analogDiscoveries } from '../drizzle/schema';
import { writeFile, unlink } from 'fs/promises';
import { join } from 'path';
import { tmpdir } from 'os';

describe('SDF Import and Export', () => {
  let testSdfPath: string;
  
  beforeAll(async () => {
    // Create a test SDF file
    testSdfPath = join(tmpdir(), 'test-ketamine-analog.sdf');
    const sdfContent = `KETAMINE-TEST-001
  PharmaSight Test
  
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
> <COMPOUND_NAME>
KETAMINE-TEST-001
> <PARENT_COMPOUND>
Ketamine
> <SMILES>
CNC1(c2cccc(F)c2Cl)CCCCC1=O
> <CONFIDENCE_SCORE>
90
> <SIMILARITY_SCORE>
88
> <SAFETY_SCORE>
90
> <EFFICACY_SCORE>
88
> <MARKET_VALUE>
$65M
> <PATENT_STATUS>
patent-free
> <THERAPEUTIC_POTENTIAL>
Enhanced NMDA receptor antagonism for treatment-resistant depression
> <KEY_DIFFERENCES>
Fluorinated derivative with stronger binding affinity
$$$$
`;
    await writeFile(testSdfPath, sdfContent, 'utf-8');
  });
  
  it('should import SDF file with ADMET analysis', async () => {
    const result = await importFromSDF(testSdfPath);
    
    expect(result.success).toBe(true);
    expect(result.imported).toBeGreaterThan(0);
    expect(result.analogs).toHaveLength(result.imported);
    
    // Check first analog has required fields
    const analog = result.analogs[0];
    expect(analog.compoundName).toBe('KETAMINE-TEST-001');
    expect(analog.smiles).toBe('CNC1(c2cccc(F)c2Cl)CCCCC1=O');
    expect(analog.confidenceScore).toBeGreaterThan(0);
    expect(analog.safetyScore).toBeGreaterThan(0);
  }, 60000);
  
  it('should export analogs to CSV format', async () => {
    const result = await exportAnalogs({
      format: 'csv',
      filters: { minConfidence: 50 },
    });
    
    expect(result.success).toBe(true);
    expect(result.count).toBeGreaterThan(0);
    expect(result.filePath).toBeDefined();
    expect(result.filePath).toContain('.csv');
    
    // Clean up
    if (result.filePath) {
      await unlink(result.filePath).catch(() => {});
    }
  }, 30000);
  
  it('should export analogs to SDF format', async () => {
    const result = await exportAnalogs({
      format: 'sdf',
      filters: { minConfidence: 50 },
    });
    
    expect(result.success).toBe(true);
    expect(result.count).toBeGreaterThan(0);
    expect(result.filePath).toBeDefined();
    expect(result.filePath).toContain('.sdf');
    
    // Clean up
    if (result.filePath) {
      await unlink(result.filePath).catch(() => {});
    }
  }, 60000);
  
  it('should filter exports by patent status', async () => {
    const result = await exportAnalogs({
      format: 'csv',
      filters: {
        patentStatus: ['patent-free'],
      },
    });
    
    expect(result.success).toBe(true);
    
    // Clean up
    if (result.filePath) {
      await unlink(result.filePath).catch(() => {});
    }
  }, 30000);
});
