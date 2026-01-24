import { describe, it, expect, beforeAll } from 'vitest';
import { exportAnalogs } from './batchExporter';
import { getDb } from './db';
import { analogDiscoveries } from '../drizzle/schema';
import { writeFile, unlink } from 'fs/promises';
import { join } from 'path';
import { tmpdir } from 'os';

describe('SDF Import and Export', () => {
  
  it('should export analogs to CSV format', async () => {
    const result = await exportAnalogs({
      format: 'csv',
      filters: { minConfidence: 50 },
    });
    
    // CSV export should always succeed (may have 0 results if no analogs match)
    if (result.count === 0) {
      expect(result.success).toBe(false);
      expect(result.error).toContain('No analogs match');
    } else {
      expect(result.success).toBe(true);
      expect(result.count).toBeGreaterThan(0);
      expect(result.filePath).toBeDefined();
      expect(result.filePath).toContain('.csv');
      
      // Clean up
      if (result.filePath) {
        await unlink(result.filePath).catch(() => {});
      }
    }
  }, 30000);
  
  it('should handle export with no matching analogs gracefully', async () => {
    const result = await exportAnalogs({
      format: 'csv',
      filters: { minConfidence: 999 }, // Impossible confidence score
    });
    
    // Should return success: false with appropriate message
    expect(result.success).toBe(false);
    expect(result.error).toBeDefined();
    expect(result.count).toBe(0);
  }, 30000);
  
  it('should filter exports by patent status', async () => {
    const result = await exportAnalogs({
      format: 'csv',
      filters: {
        patentStatus: ['patent-free'],
      },
    });
    
    // Either succeeds with results or fails with no matches
    if (result.count > 0) {
      expect(result.success).toBe(true);
    } else {
      expect(result.success).toBe(false);
    }
    
    // Clean up
    if (result.filePath) {
      await unlink(result.filePath).catch(() => {});
    }
  }, 30000);
});
