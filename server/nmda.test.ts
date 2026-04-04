import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { getDb } from './db';
import { receptorLibrary } from '../drizzle/schema';
import { eq } from 'drizzle-orm';
import * as fs from 'fs';
import * as path from 'path';

let db: any;

describe('NMDA Receptor Integration', () => {
  const nmda8xlk = {
    pdbId: '8XLK',
    targetName: 'NMDA-GluN2A-GluN2B',
    pdbFileName: '8xlk.pdb',
    pdbqtFileName: '8xlk.pdbqt',
    pdbqtUrl: 's3://receptors/nmda/8xlk.pdbqt',
    description: 'Native tri-heteromeric NMDA receptor from rat cortex and hippocampus',
    resolution: '2.85',
    experimentalMethod: 'CRYO-EM',
    organism: 'Rattus norvegicus',
    category: 'psychiatric',
    tags: JSON.stringify(['nmda', 'glun2a', 'glun2b', 'glutamate', 'rat']),
    notes: 'Tri-heteromeric composition found in cortex and hippocampus. Important for synaptic plasticity.',
    source: 'RCSB PDB',
    isActive: true,
  };

  const nmda9jnn = {
    pdbId: '9JNN',
    targetName: 'NMDA-GluN2B',
    pdbFileName: '9jnn.pdb',
    pdbqtFileName: '9jnn.pdbqt',
    pdbqtUrl: 's3://receptors/nmda/9jnn.pdbqt',
    description: 'Native di-heteromeric NMDA receptor from rat cortex and hippocampus',
    resolution: '2.95',
    experimentalMethod: 'CRYO-EM',
    organism: 'Rattus norvegicus',
    category: 'psychiatric',
    tags: JSON.stringify(['nmda', 'glun2b', 'glutamate', 'rat']),
    notes: 'Enriched in developing neurons and extrasynaptic locations. Target for neuroprotection and pain.',
    source: 'RCSB PDB',
    isActive: true,
  };

  beforeAll(async () => {
    db = await getDb();
    if (!db) throw new Error('Database not available');
    // Clean up any existing NMDA entries
    await db.delete(receptorLibrary).where(eq(receptorLibrary.pdbId, '8XLK'));
    await db.delete(receptorLibrary).where(eq(receptorLibrary.pdbId, '9JNN'));
  });

  afterAll(async () => {
    if (!db) return;
    // Clean up after tests
    await db.delete(receptorLibrary).where(eq(receptorLibrary.pdbId, '8XLK'));
    await db.delete(receptorLibrary).where(eq(receptorLibrary.pdbId, '9JNN'));
  });

  it('should insert NMDA 8XLK receptor', async () => {
    const result = await db.insert(receptorLibrary).values(nmda8xlk);
    expect(result).toBeDefined();

    const inserted = await db
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.pdbId, '8XLK'));

    expect(inserted).toHaveLength(1);
    expect(inserted[0].targetName).toBe('NMDA-GluN2A-GluN2B');
    expect(inserted[0].pdbId).toBe('8XLK');
  });

  it('should insert NMDA 9JNN receptor', async () => {
    const result = await db.insert(receptorLibrary).values(nmda9jnn);
    expect(result).toBeDefined();

    const inserted = await db
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.pdbId, '9JNN'));

    expect(inserted).toHaveLength(1);
    expect(inserted[0].targetName).toBe('NMDA-GluN2B');
    expect(inserted[0].pdbId).toBe('9JNN');
  });

  it('should verify PDBQT files exist', () => {
    const pdbqtDir = path.join(process.cwd(), 'data', 'nmda_receptors');
    const file8xlk = path.join(pdbqtDir, '8xlk.pdbqt');
    const file9jnn = path.join(pdbqtDir, '9jnn.pdbqt');

    expect(fs.existsSync(file8xlk)).toBe(true);
    expect(fs.existsSync(file9jnn)).toBe(true);

    const size8xlk = fs.statSync(file8xlk).size / 1024;
    const size9jnn = fs.statSync(file9jnn).size / 1024;

    expect(size8xlk).toBeGreaterThan(2000); // Should be ~2.6 MB
    expect(size9jnn).toBeGreaterThan(2000); // Should be ~2.5 MB

    console.log(`✅ NMDA 8XLK PDBQT: ${size8xlk.toFixed(1)} KB`);
    console.log(`✅ NMDA 9JNN PDBQT: ${size9jnn.toFixed(1)} KB`);
  });

  it('should query all NMDA receptors', async () => {
    const nmda = await db
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.category, 'psychiatric'));

    expect(nmda.length).toBeGreaterThanOrEqual(2);
    expect(nmda.some((r) => r.pdbId === '8XLK')).toBe(true);
    expect(nmda.some((r) => r.pdbId === '9JNN')).toBe(true);
  });

  it('should filter NMDA receptors by organism', async () => {
    const ratReceptors = await db
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.organism, 'Rattus norvegicus'));

    expect(ratReceptors.length).toBeGreaterThan(0);
    expect(ratReceptors.some((r) => r.pdbId === '8XLK')).toBe(true);
  });

  it('should verify receptor metadata completeness', async () => {
    const receptors = await db
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.category, 'psychiatric'));

    receptors.forEach((receptor) => {
      expect(receptor.pdbId).toBeDefined();
      expect(receptor.targetName).toBeDefined();
      expect(receptor.category).toBe('psychiatric');
      if (receptor.organism) {
        expect(receptor.organism).toBe('Rattus norvegicus');
      }
      expect(receptor.tags).toBeDefined();
    });
  });
});
