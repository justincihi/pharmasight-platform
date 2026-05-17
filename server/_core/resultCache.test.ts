import { describe, it, expect, beforeEach } from 'vitest';
import {
  ResultCacheManager,
  dockingCache,
  toxicityCache,
  getCacheForAnalysisType,
  getAllCacheStats,
  clearAllCaches,
} from './resultCache';

describe('Result Cache System', () => {
  let cache: ResultCacheManager<Record<string, unknown>>;

  beforeEach(() => {
    cache = new ResultCacheManager(3600, 100);
  });

  describe('Basic Cache Operations', () => {
    it('should set and get cached result', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O'; // Aspirin
      const result = { binding_affinity: -6.5, score: 75 };

      cache.set(smiles, result);
      const cached = cache.get(smiles);

      expect(cached).toEqual(result);
    });

    it('should return null for uncached SMILES', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const cached = cache.get(smiles);

      expect(cached).toBeNull();
    });

    it('should check if result is cached', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      expect(cache.has(smiles)).toBe(false);

      cache.set(smiles, result);

      expect(cache.has(smiles)).toBe(true);
    });

    it('should clear cache', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles, result);
      expect(cache.has(smiles)).toBe(true);

      cache.clear();

      expect(cache.has(smiles)).toBe(false);
    });
  });

  describe('SMILES Normalization', () => {
    it('should treat uppercase and lowercase SMILES as same', () => {
      const smiles1 = 'CC(=O)Oc1ccccc1C(=O)O';
      const smiles2 = 'cc(=o)oc1ccccc1c(=o)o';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles1, result);
      const cached = cache.get(smiles2);

      expect(cached).toEqual(result);
    });

    it('should handle SMILES with whitespace', () => {
      const smiles1 = 'CC(=O)Oc1ccccc1C(=O)O';
      const smiles2 = '  CC(=O)Oc1ccccc1C(=O)O  ';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles1, result);
      const cached = cache.get(smiles2);

      expect(cached).toEqual(result);
    });
  });

  describe('Cache Statistics', () => {
    it('should track cache hits', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles, result);

      // Access 3 times
      cache.get(smiles);
      cache.get(smiles);
      cache.get(smiles);

      const stats = cache.getStats();

      expect(stats.totalHits).toBe(3);
      expect(stats.totalEntries).toBe(1);
    });

    it('should calculate hit rate', () => {
      const smiles1 = 'CC(=O)Oc1ccccc1C(=O)O';
      const smiles2 = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles1, result);

      // 3 hits on cached item
      cache.get(smiles1);
      cache.get(smiles1);
      cache.get(smiles1);

      // 2 misses on uncached items
      cache.get(smiles2);
      cache.get(smiles2);

      const stats = cache.getStats();

      // Hit rate = 3 / (3 + 2) = 0.6
      expect(stats.hitRate).toBeCloseTo(0.6, 1);
    });

    it('should track average age of entries', () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles, result);

      const stats = cache.getStats();

      expect(stats.averageAge).toBeGreaterThanOrEqual(0);
      expect(stats.averageAge).toBeLessThan(1000); // Should be very recent
    });
  });

  describe('Cache Eviction', () => {
    it('should evict oldest entry when max size reached', () => {
      const smallCache = new ResultCacheManager(3600, 2);

      const smiles1 = 'CC(=O)Oc1ccccc1C(=O)O';
      const smiles2 = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O';
      const smiles3 = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C';

      const result = { binding_affinity: -6.5 };

      // Add 3 items to cache with max size 2
      smallCache.set(smiles1, result);
      smallCache.set(smiles2, result);
      smallCache.set(smiles3, result);

      const stats = smallCache.getStats();

      // Should only have 2 entries
      expect(stats.totalEntries).toBe(2);

      // First entry should be evicted
      expect(smallCache.has(smiles1)).toBe(false);
    });
  });

  describe('Get or Compute', () => {
    it('should return cached result without computing', async () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles, result);

      let computeCalled = false;
      const computed = await cache.getOrCompute(smiles, async () => {
        computeCalled = true;
        return { binding_affinity: -8.0 };
      });

      expect(computeCalled).toBe(false);
      expect(computed).toEqual(result);
    });

    it('should compute and cache result if not cached', async () => {
      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const computedResult = { binding_affinity: -8.0 };

      let computeCalled = false;
      const result = await cache.getOrCompute(smiles, async () => {
        computeCalled = true;
        return computedResult;
      });

      expect(computeCalled).toBe(true);
      expect(result).toEqual(computedResult);

      // Should be cached now
      expect(cache.has(smiles)).toBe(true);
    });
  });

  describe('Global Cache Instances', () => {
    it('should get correct cache for analysis type', () => {
      const docking = getCacheForAnalysisType('docking');
      const toxicity = getCacheForAnalysisType('toxicity');

      expect(docking).toBe(dockingCache);
      expect(toxicity).toBe(toxicityCache);
    });

    it('should throw error for unknown analysis type', () => {
      expect(() => {
        getCacheForAnalysisType('unknown' as any);
      }).toThrow();
    });

    it('should get all cache statistics', () => {
      const stats = getAllCacheStats();

      expect(stats).toHaveProperty('docking');
      expect(stats).toHaveProperty('toxicity');
      expect(stats).toHaveProperty('admet');
      expect(stats).toHaveProperty('pkpd');
    });

    it('should clear all caches', () => {
      dockingCache.set('CC(=O)Oc1ccccc1C(=O)O', { test: true });
      toxicityCache.set('CC(=O)Oc1ccccc1C(=O)O', { test: true });

      expect(dockingCache.has('CC(=O)Oc1ccccc1C(=O)O')).toBe(true);
      expect(toxicityCache.has('CC(=O)Oc1ccccc1C(=O)O')).toBe(true);

      clearAllCaches();

      expect(dockingCache.has('CC(=O)Oc1ccccc1C(=O)O')).toBe(false);
      expect(toxicityCache.has('CC(=O)Oc1ccccc1C(=O)O')).toBe(false);
    });
  });

  describe('Cache Configuration', () => {
    it('should allow setting TTL', () => {
      cache.setTTL(60); // 60 seconds

      const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles, result);

      // Should be cached
      expect(cache.has(smiles)).toBe(true);
    });

    it('should allow setting max entries', () => {
      cache.setMaxEntries(5);

      const result = { binding_affinity: -6.5 };

      for (let i = 0; i < 10; i++) {
        cache.set(`SMILES${i}`, result);
      }

      const stats = cache.getStats();

      // Should not exceed max entries
      expect(stats.totalEntries).toBeLessThanOrEqual(5);
    });
  });

  describe('Get All Entries', () => {
    it('should return all cached entries sorted by hits', () => {
      const smiles1 = 'CC(=O)Oc1ccccc1C(=O)O';
      const smiles2 = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O';
      const result = { binding_affinity: -6.5 };

      cache.set(smiles1, result);
      cache.set(smiles2, result);

      // Access smiles1 more times
      cache.get(smiles1);
      cache.get(smiles1);
      cache.get(smiles1);

      const entries = cache.getAll();

      expect(entries.length).toBe(2);
      expect(entries[0].hits).toBeGreaterThanOrEqual(entries[1].hits);
    });
  });
});
