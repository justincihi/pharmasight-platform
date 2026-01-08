/**
 * Base Plugin Interface for PharmaSight Extensions
 * 
 * This provides a clean, extensible architecture for adding new
 * data sources, simulation engines, and analysis tools without
 * modifying core code.
 */

export interface PluginMetadata {
  name: string;
  version: string;
  description: string;
  author?: string;
  license?: string;
  dependencies?: string[];
}

export interface PluginConfig {
  enabled: boolean;
  [key: string]: any;
}

export interface PluginResult<T = any> {
  success: boolean;
  data?: T;
  error?: string;
  metadata?: {
    executionTime?: number;
    source?: string;
    cached?: boolean;
    [key: string]: any;
  };
}

/**
 * Base Plugin Interface
 * All plugins must implement these methods
 */
export interface Plugin<TInput = any, TOutput = any> {
  /** Plugin metadata */
  readonly metadata: PluginMetadata;
  
  /** Plugin configuration */
  config: PluginConfig;
  
  /**
   * Initialize the plugin
   * Called once when the plugin is loaded
   */
  initialize(): Promise<void>;
  
  /**
   * Validate input before processing
   * @param input - Input data to validate
   * @returns Validation result with error message if invalid
   */
  validate(input: TInput): Promise<{ valid: boolean; error?: string }>;
  
  /**
   * Execute the plugin's main functionality
   * @param input - Input data
   * @returns Plugin result with data or error
   */
  run(input: TInput): Promise<PluginResult<TOutput>>;
  
  /**
   * Summarize results in a human-readable format
   * @param result - Plugin execution result
   * @returns Summary string or object
   */
  summarizeResults(result: PluginResult<TOutput>): string | object;
  
  /**
   * Cleanup resources
   * Called when the plugin is unloaded or server shuts down
   */
  cleanup?(): Promise<void>;
}

/**
 * Data Source Plugin Interface
 * For plugins that fetch data from external databases/APIs
 */
export interface DataSourcePlugin<TQuery = any, TData = any> extends Plugin<TQuery, TData> {
  /**
   * Check if the data source is available
   */
  isAvailable(): Promise<boolean>;
  
  /**
   * Get rate limit information
   */
  getRateLimit?(): Promise<{
    limit: number;
    remaining: number;
    reset: Date;
  }>;
}

/**
 * Simulation Plugin Interface
 * For plugins that run computational simulations (docking, PBPK, etc.)
 */
export interface SimulationPlugin<TInput = any, TOutput = any> extends Plugin<TInput, TOutput> {
  /**
   * Estimate execution time for the simulation
   * @param input - Input parameters
   * @returns Estimated time in milliseconds
   */
  estimateExecutionTime(input: TInput): Promise<number>;
  
  /**
   * Cancel a running simulation
   * @param jobId - Job identifier
   */
  cancel?(jobId: string): Promise<void>;
  
  /**
   * Get simulation progress
   * @param jobId - Job identifier
   * @returns Progress percentage (0-100)
   */
  getProgress?(jobId: string): Promise<number>;
}

/**
 * Analysis Plugin Interface
 * For plugins that analyze molecular properties, toxicity, etc.
 */
export interface AnalysisPlugin<TInput = any, TOutput = any> extends Plugin<TInput, TOutput> {
  /**
   * Get supported analysis types
   */
  getSupportedAnalyses(): string[];
  
  /**
   * Batch analysis support
   * @param inputs - Array of inputs to analyze
   * @returns Array of results
   */
  runBatch?(inputs: TInput[]): Promise<PluginResult<TOutput>[]>;
}

/**
 * Plugin Registry
 * Manages all loaded plugins
 */
export class PluginRegistry {
  private plugins: Map<string, Plugin> = new Map();
  
  /**
   * Register a plugin
   */
  register(plugin: Plugin): void {
    const name = plugin.metadata.name;
    if (this.plugins.has(name)) {
      throw new Error(`Plugin "${name}" is already registered`);
    }
    this.plugins.set(name, plugin);
  }
  
  /**
   * Get a plugin by name
   */
  get<T extends Plugin = Plugin>(name: string): T | undefined {
    return this.plugins.get(name) as T | undefined;
  }
  
  /**
   * Get all registered plugins
   */
  getAll(): Plugin[] {
    return Array.from(this.plugins.values());
  }
  
  /**
   * Get plugins by type
   */
  getByType<T extends Plugin>(type: new (...args: any[]) => T): T[] {
    return this.getAll().filter(p => p instanceof type) as T[];
  }
  
  /**
   * Unregister a plugin
   */
  async unregister(name: string): Promise<void> {
    const plugin = this.plugins.get(name);
    if (plugin && plugin.cleanup) {
      await plugin.cleanup();
    }
    this.plugins.delete(name);
  }
  
  /**
   * Initialize all plugins
   */
  async initializeAll(): Promise<void> {
    const promises = Array.from(this.plugins.values()).map(p => p.initialize());
    await Promise.all(promises);
  }
  
  /**
   * Cleanup all plugins
   */
  async cleanupAll(): Promise<void> {
    const promises = Array.from(this.plugins.values())
      .filter(p => p.cleanup)
      .map(p => p.cleanup!());
    await Promise.all(promises);
  }
}

// Global plugin registry instance
export const pluginRegistry = new PluginRegistry();
