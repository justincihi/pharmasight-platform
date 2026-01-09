/**
 * Master Analog File Synchronization Service
 * Keeps master_analogs.json in sync with database
 */

import fs from 'fs/promises';
import path from 'path';
import { getDb } from './db';
import { analogDiscoveries } from '../drizzle/schema';

const MASTER_FILE_PATH = path.join(process.cwd(), 'master_analogs.json');

export interface MasterAnalog {
  compoundId: string;
  compoundName: string;
  smiles: string;
  parentCompound: string;
  mechanismOfAction: string;
  keyDifferences: string;
  confidence: number;
  similarity: number;
  safetyScore: number;
  efficacyScore: number;
  drugLikenessScore: number;
  patentStatus: string;
  marketValue: string;
  discoveryMethod: string;
  discoveredAt?: string;
}

/**
 * Export all analogs from database to master JSON file
 */
export async function exportToMasterFile(): Promise<{ success: boolean; count: number }> {
  try {
    const db = await getDb();
    if (!db) throw new Error('Database not available');
    const analogs = await db.select().from(analogDiscoveries);
    
    const masterAnalogs: MasterAnalog[] = analogs.map((analog: any) => ({
      compoundId: analog.compoundId,
      compoundName: analog.compoundName,
      smiles: analog.smiles,
      parentCompound: analog.parentCompound,
      mechanismOfAction: analog.mechanismOfAction || '',
      keyDifferences: analog.keyDifferences || '',
      confidence: analog.confidenceScore,
      similarity: analog.similarityScore,
      safetyScore: analog.safetyScore,
      efficacyScore: analog.efficacyScore,
      drugLikenessScore: analog.drugLikenessScore,
      patentStatus: analog.patentStatus,
      marketValue: analog.marketValue || '$0',
      discoveryMethod: analog.discoveryMethod || 'unknown',
      discoveredAt: analog.discoveredAt?.toISOString()
    }));
    
    await fs.writeFile(
      MASTER_FILE_PATH,
      JSON.stringify(masterAnalogs, null, 2),
      'utf-8'
    );
    
    console.log(`[MasterFileSync] Exported ${masterAnalogs.length} analogs to master file`);
    
    return { success: true, count: masterAnalogs.length };
  } catch (error) {
    console.error('[MasterFileSync] Export error:', error);
    return { success: false, count: 0 };
  }
}

/**
 * Import analogs from master JSON file to database
 */
export async function importFromMasterFile(): Promise<{ success: boolean; imported: number; skipped: number }> {
  try {
    const fileContent = await fs.readFile(MASTER_FILE_PATH, 'utf-8');
    const masterAnalogs: MasterAnalog[] = JSON.parse(fileContent);
    
    let imported = 0;
    let skipped = 0;
    
    for (const analog of masterAnalogs) {
      try {
        // Check if already exists
        const db = await getDb();
        if (!db) continue;
        
        const existing = await db
          .select()
          .from(analogDiscoveries)
          .where((fields: any) => fields.compoundId.eq(analog.compoundId))
          .limit(1);
        
        if (existing.length > 0) {
          skipped++;
          continue;
        }
        
        // Insert new analog
        await db.insert(analogDiscoveries).values({
          compoundId: analog.compoundId,
          compoundName: analog.compoundName,
          smiles: analog.smiles,
          parentCompound: analog.parentCompound,
          mechanismOfAction: analog.mechanismOfAction,
          keyDifferences: analog.keyDifferences,
          confidenceScore: analog.confidence,
          similarityScore: typeof analog.similarity === 'number' && analog.similarity <= 1 ? Math.round(analog.similarity * 100) : analog.similarity,
          safetyScore: analog.safetyScore,
          efficacyScore: analog.efficacyScore,
          drugLikenessScore: analog.drugLikenessScore,
          patentStatus: analog.patentStatus as any,
          marketValue: analog.marketValue,
          discoveryMethod: analog.discoveryMethod,
          discoveredBy: 'master-file-import',
          discoveredAt: analog.discoveredAt ? new Date(analog.discoveredAt) : new Date()
        });
        
        imported++;
      } catch (err) {
        console.error(`[MasterFileSync] Error importing ${analog.compoundId}:`, err);
        skipped++;
      }
    }
    
    console.log(`[MasterFileSync] Import complete: ${imported} imported, ${skipped} skipped`);
    
    return { success: true, imported, skipped };
  } catch (error) {
    console.error('[MasterFileSync] Import error:', error);
    return { success: false, imported: 0, skipped: 0 };
  }
}

/**
 * Sync database changes to master file (called after DB updates)
 */
export async function syncToMasterFile(): Promise<void> {
  try {
    await exportToMasterFile();
  } catch (error) {
    console.error('[MasterFileSync] Sync error:', error);
  }
}

/**
 * Check if master file exists
 */
export async function masterFileExists(): Promise<boolean> {
  try {
    await fs.access(MASTER_FILE_PATH);
    return true;
  } catch {
    return false;
  }
}
