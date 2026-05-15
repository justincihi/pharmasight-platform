/**
 * Simple in-memory cache with TTL support
 * Used for caching API responses and expensive computations
 */

interface CacheEntry<T> {
  value: T;
  expiresAt: number;
}

class Cache {
  private store: Map<string, CacheEntry<unknown>> = new Map();

  /**
   * Get value from cache
   */
  get<T>(key: string): T | null {
    const entry = this.store.get(key) as CacheEntry<T> | undefined;

    if (!entry) {
      return null;
    }

    // Check if expired
    if (Date.now() > entry.expiresAt) {
      this.store.delete(key);
      return null;
    }

    return entry.value;
  }

  /**
   * Set value in cache with TTL
   */
  set<T>(key: string, value: T, ttlSeconds: number = 3600): void {
    const expiresAt = Date.now() + ttlSeconds * 1000;
    this.store.set(key, { value, expiresAt });
  }

  /**
   * Delete value from cache
   */
  delete(key: string): boolean {
    return this.store.delete(key);
  }

  /**
   * Clear all cache
   */
  clear(): void {
    this.store.clear();
  }

  /**
   * Get cache size
   */
  size(): number {
    return this.store.size;
  }

  /**
   * Clean up expired entries
   */
  cleanup(): number {
    const now = Date.now();
    let cleaned = 0;

    const entries = Array.from(this.store.entries());
    for (const [key, entry] of entries) {
      if (now > entry.expiresAt) {
        this.store.delete(key);
        cleaned++;
      }
    }

    return cleaned;
  }
}

// Export singleton instance
export const cache = new Cache();

// Cleanup expired entries every 5 minutes
setInterval(() => {
  const cleaned = cache.cleanup();
  if (cleaned > 0) {
    console.log(`[Cache] Cleaned up ${cleaned} expired entries`);
  }
}, 5 * 60 * 1000);
