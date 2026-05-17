/**
 * Result Cache System
 * Caches analysis results by SMILES hash to optimize performance
 * Reduces redundant computations for identical compounds
 */

import crypto from 'crypto';

interface CacheEntry<T> {
  hash: string;
  smiles: string;
  result: T;
  timestamp: number;
  expiresAt: number;
  hits: number;
}

interface CacheStats {
  totalEntries: number;
  totalHits: number;
  hitRate: number;
  averageAge: number;
  cacheSize: number;
}

/**
 * Generate SHA256 hash of SMILES string
 */
function generateSmilesHash(smiles: string): string {
  return crypto
    .createHash('sha256')
    .update(smiles.toLowerCase().trim())
    .digest('hex');
}

/**
 * Result Cache Manager
 * Manages caching of analysis results with TTL and statistics
 */
export class ResultCacheManager<T> {
  private cache: Map<string, CacheEntry<T>> = new Map();
  private ttl: number; // Time to live in milliseconds
  private maxEntries: number;
  private totalHits: number = 0;
  private totalRequests: number = 0;

  constructor(ttlSeconds: number = 3600, maxEntries: number = 1000) {
    this.ttl = ttlSeconds * 1000;
    this.maxEntries = maxEntries;

    // Cleanup expired entries periodically
    setInterval(() => this.cleanup(), 60000); // Every minute
  }

  /**
   * Get cached result by SMILES
   */
  get(smiles: string): T | null {
    const hash = generateSmilesHash(smiles);
    const entry = this.cache.get(hash);

    this.totalRequests++;

    if (!entry) {
      return null;
    }

    // Check if expired
    if (entry.expiresAt < Date.now()) {
      this.cache.delete(hash);
      return null;
    }

    // Update hit count
    entry.hits++;
    this.totalHits++;

    return entry.result;
  }

  /**
   * Set cached result by SMILES
   */
  set(smiles: string, result: T): void {
    const hash = generateSmilesHash(smiles);

    // Check cache size limit
    if (this.cache.size >= this.maxEntries && !this.cache.has(hash)) {
      // Remove oldest entry
      let oldestHash = '';
      let oldestTime = Date.now();

      this.cache.forEach((entry, key) => {
        if (entry.timestamp < oldestTime) {
          oldestTime = entry.timestamp;
          oldestHash = key;
        }
      });

      if (oldestHash) {
        this.cache.delete(oldestHash);
      }
    }

    const entry: CacheEntry<T> = {
      hash,
      smiles: smiles.toLowerCase().trim(),
      result,
      timestamp: Date.now(),
      expiresAt: Date.now() + this.ttl,
      hits: 0,
    };

    this.cache.set(hash, entry);
  }

  /**
   * Check if result is cached
   */
  has(smiles: string): boolean {
    const hash = generateSmilesHash(smiles);
    const entry = this.cache.get(hash);

    if (!entry) {
      return false;
    }

    // Check if expired
    if (entry.expiresAt < Date.now()) {
      this.cache.delete(hash);
      return false;
    }

    return true;
  }

  /**
   * Get or compute result
   */
  async getOrCompute(
    smiles: string,
    compute: (smiles: string) => Promise<T>
  ): Promise<T> {
    // Try to get from cache
    const cached = this.get(smiles);
    if (cached !== null) {
      return cached;
    }

    // Compute and cache
    const result = await compute(smiles);
    this.set(smiles, result);

    return result;
  }

  /**
   * Clear cache
   */
  clear(): void {
    this.cache.clear();
    this.totalHits = 0;
    this.totalRequests = 0;
  }

  /**
   * Remove expired entries
   */
  private cleanup(): void {
    const now = Date.now();
    const keysToDelete: string[] = [];

    this.cache.forEach((entry, hash) => {
      if (entry.expiresAt < now) {
        keysToDelete.push(hash);
      }
    });

    keysToDelete.forEach((hash) => this.cache.delete(hash));
  }

  /**
   * Get cache statistics
   */
  getStats(): CacheStats {
    let totalAge = 0;
    let cacheSize = 0;

    this.cache.forEach((entry) => {
      totalAge += Date.now() - entry.timestamp;
      cacheSize += JSON.stringify(entry.result).length;
    });

    const averageAge = this.cache.size > 0 ? totalAge / this.cache.size : 0;
    const hitRate = this.totalRequests > 0 ? this.totalHits / this.totalRequests : 0;

    return {
      totalEntries: this.cache.size,
      totalHits: this.totalHits,
      hitRate,
      averageAge,
      cacheSize,
    };
  }

  /**
   * Get all cached entries (for debugging)
   */
  getAll(): Array<{ smiles: string; hash: string; hits: number; age: number }> {
    const entries: Array<{ smiles: string; hash: string; hits: number; age: number }> = [];

    this.cache.forEach((entry) => {
      entries.push({
        smiles: entry.smiles,
        hash: entry.hash,
        hits: entry.hits,
        age: Date.now() - entry.timestamp,
      });
    });

    return entries.sort((a, b) => b.hits - a.hits);
  }

  /**
   * Set TTL for cache entries
   */
  setTTL(seconds: number): void {
    this.ttl = seconds * 1000;
  }

  /**
   * Set maximum number of entries
   */
  setMaxEntries(max: number): void {
    this.maxEntries = max;
  }
}

/**
 * Global cache instances for different analysis types
 */
export const dockingCache = new ResultCacheManager<Record<string, unknown>>(3600); // 1 hour
export const toxicityCache = new ResultCacheManager<Record<string, unknown>>(3600);
export const admetCache = new ResultCacheManager<Record<string, unknown>>(3600);
export const pkpdCache = new ResultCacheManager<Record<string, unknown>>(3600);

/**
 * Utility function to get cache for analysis type
 */
export function getCacheForAnalysisType(
  type: 'docking' | 'toxicity' | 'admet' | 'pkpd'
): ResultCacheManager<Record<string, unknown>> {
  switch (type) {
    case 'docking':
      return dockingCache;
    case 'toxicity':
      return toxicityCache;
    case 'admet':
      return admetCache;
    case 'pkpd':
      return pkpdCache;
    default:
      throw new Error(`Unknown analysis type: ${type}`);
  }
}

/**
 * Utility function to get all cache statistics
 */
export function getAllCacheStats(): Record<string, CacheStats> {
  return {
    docking: dockingCache.getStats(),
    toxicity: toxicityCache.getStats(),
    admet: admetCache.getStats(),
    pkpd: pkpdCache.getStats(),
  };
}

/**
 * Utility function to clear all caches
 */
export function clearAllCaches(): void {
  dockingCache.clear();
  toxicityCache.clear();
  admetCache.clear();
  pkpdCache.clear();
}
