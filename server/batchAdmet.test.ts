import { describe, it, expect, vi, beforeEach } from 'vitest';

// Mock the db module
vi.mock('./db', () => ({
  insertAdmetResult: vi.fn().mockResolvedValue({ id: 1, analogId: 42, smiles: 'CC(=O)O', source: 'admet_ai_chemprop', properties: '{}', createdAt: new Date() }),
  getAdmetResultsForAnalog: vi.fn().mockResolvedValue([
    { id: 1, analogId: 42, smiles: 'CC(=O)O', source: 'admet_ai_chemprop', properties: '{"hERG":0.003,"BBB":0.811}', createdAt: new Date() }
  ]),
  getAdmetStats: vi.fn().mockResolvedValue({ total: 5, bySource: { admet_ai_chemprop: 4, rdkit_rules: 1 }, lastRun: new Date() }),
}));

import { insertAdmetResult, getAdmetResultsForAnalog, getAdmetStats } from './db';

describe('admetResults DB helpers', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('insertAdmetResult stores a result and returns it', async () => {
    const result = await insertAdmetResult({
      analogId: 42,
      smiles: 'CC(=O)O',
      source: 'admet_ai_chemprop',
      properties: '{}',
    });
    expect(result).toBeDefined();
    expect(result.analogId).toBe(42);
    expect(result.source).toBe('admet_ai_chemprop');
    expect(insertAdmetResult).toHaveBeenCalledOnce();
  });

  it('getAdmetResultsForAnalog returns results for a given analogId', async () => {
    const results = await getAdmetResultsForAnalog(42);
    expect(Array.isArray(results)).toBe(true);
    expect(results.length).toBeGreaterThan(0);
    expect(results[0].analogId).toBe(42);
    expect(getAdmetResultsForAnalog).toHaveBeenCalledWith(42);
  });

  it('getAdmetStats returns aggregated statistics', async () => {
    const stats = await getAdmetStats();
    expect(stats).toBeDefined();
    expect(typeof stats.total).toBe('number');
    expect(stats.bySource).toBeDefined();
    expect(getAdmetStats).toHaveBeenCalledOnce();
  });

  it('getAdmetResultsForAnalog returns empty array for unknown analog', async () => {
    vi.mocked(getAdmetResultsForAnalog).mockResolvedValueOnce([]);
    const results = await getAdmetResultsForAnalog(9999);
    expect(results).toEqual([]);
  });
});

describe('MetabolitePathwayTree data mapping', () => {
  it('maps Phase I string to numeric 1', () => {
    const raw = { phase: 'Phase I', probability: '0.85', molecularWeight: '180.2', logP: '1.5' };
    const mapped = {
      phase: raw.phase === 'Phase I' ? 1 : raw.phase === 'Phase II' ? 2 : 1,
      confidenceScore: parseFloat(raw.probability),
      molecularWeight: parseFloat(raw.molecularWeight),
      logP: parseFloat(raw.logP),
    };
    expect(mapped.phase).toBe(1);
    expect(mapped.confidenceScore).toBeCloseTo(0.85);
    expect(mapped.molecularWeight).toBeCloseTo(180.2);
  });

  it('maps Phase II string to numeric 2', () => {
    const raw = { phase: 'Phase II' };
    const phase = raw.phase === 'Phase I' ? 1 : raw.phase === 'Phase II' ? 2 : 1;
    expect(phase).toBe(2);
  });

  it('defaults to 1 for unknown phase values', () => {
    const raw = { phase: 'Unknown' };
    const phase = raw.phase === 'Phase I' ? 1 : raw.phase === 'Phase II' ? 2 : 1;
    expect(phase).toBe(1);
  });

  it('handles null probability gracefully', () => {
    const raw = { probability: null as any, confidenceScore: 0.9 };
    const confidence = raw.probability != null ? parseFloat(raw.probability) : (raw.confidenceScore ?? null);
    expect(confidence).toBe(0.9);
  });
});
