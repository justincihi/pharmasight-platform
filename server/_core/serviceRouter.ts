/**
 * Service Router
 * Routes analysis requests to appropriate backend:
 * - Python Microservice (heavy compute: docking, toxicity, PK/PD)
 * - External APIs (quick lookups: SMILES validation, similarity)
 * - Local Node.js (fallback, visualization)
 */

import { validateSmiles, searchSimilarCompounds, checkPatentStatus, getCompoundProperties } from "./pubchemApi";
import { callDockingService, callADMETService, callToxicityService, callMetaboliteService, callLeadOptimizationService, callBioNemoService } from "./pythonServiceGateway";

interface ServiceConfig {
  pythonServiceUrl?: string;
  pythonServiceTimeout?: number;
  enablePythonService: boolean;
  enableExternalApis: boolean;
}

const defaultConfig: ServiceConfig = {
  pythonServiceUrl: process.env.PYTHON_SERVICE_URL || "http://localhost:5000",
  pythonServiceTimeout: 120000, // 120 seconds
  enablePythonService: process.env.ENABLE_PYTHON_SERVICE !== "false",
  enableExternalApis: process.env.ENABLE_EXTERNAL_APIS !== "false",
};

/**
 * Check if Python microservice is available
 */
async function checkPythonServiceHealth(
  config: ServiceConfig
): Promise<boolean> {
  if (!config.enablePythonService || !config.pythonServiceUrl) {
    return false;
  }

  try {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 5000);
    try {
      const response = await fetch(`${config.pythonServiceUrl}/health`, {
        signal: controller.signal,
      });
      clearTimeout(timeoutId);
      return response.ok;
    } catch {
      clearTimeout(timeoutId);
      return false;
    }
  } catch {
    return false;
  }
}

/**
 * Route docking request to Python service
 */
export async function routeDocking(
  smiles: string,
  receptorId: string,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "python" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  // Use Python service gateway (handles HTTP calls + fallback to mock)
  try {
    const result = await callDockingService(smiles, receptorId);
    return {
      source: result.isDemo ? "fallback" : "python",
      success: true,
      data: result,
    };
  } catch (error) {
    console.error("[ServiceRouter] Docking error:", error);
    return {
      source: "fallback",
      success: false,
      error: "Docking service unavailable",
    };
  }
}

/**
 * Route toxicity prediction to Python service
 */
export async function routeToxicityPrediction(
  smiles: string,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "python" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  // Use Python service gateway (handles HTTP calls + fallback to mock)
  try {
    const result = await callToxicityService(smiles);
    return {
      source: result.isDemo ? "fallback" : "python",
      success: true,
      data: result,
    };
  } catch (error) {
    console.error("[ServiceRouter] Toxicity error:", error);
    return {
      source: "fallback",
      success: false,
      error: "Toxicity service unavailable",
    };
  }
}

/**
 * Route PK/PD simulation to Python service
 */
export async function routePkpdSimulation(
  smiles: string,
  dose: number = 100,
  route: string = "oral",
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "python" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  // Try Python service first
  if (config.enablePythonService) {
    try {
      const isHealthy = await checkPythonServiceHealth(config);
      if (isHealthy) {
        const response = await fetch(
          `${config.pythonServiceUrl}/api/pkpd/simulate`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ smiles, dose, route }),
            signal: AbortSignal.timeout(config.pythonServiceTimeout || 120000),
          }
        );

        if (response.ok) {
          const data = await response.json();
          return { source: "python", success: true, data };
        }
      }
    } catch (error) {
      console.error("[ServiceRouter] Python PK/PD error:", error);
    }
  }

  // Fallback: return error
  return {
    source: "fallback",
    success: false,
    error: "PK/PD service unavailable. Please ensure Python microservice is running.",
  };
}

/**
 * Route SMILES validation to external API
 */
export async function routeSmilesValidation(
  smiles: string,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "api" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  if (!config.enableExternalApis) {
    return {
      source: "fallback",
      success: false,
      error: "External APIs disabled",
    };
  }

  try {
    const result = await validateSmiles(smiles);
    return {
      source: "api",
      success: result.valid,
      data: result,
      error: result.error,
    };
  } catch (error) {
    return {
      source: "fallback",
      success: false,
      error: `SMILES validation error: ${error instanceof Error ? error.message : String(error)}`,
    };
  }
}

/**
 * Route similarity search to external API
 */
export async function routeSimilaritySearch(
  smiles: string,
  threshold: number = 0.9,
  maxResults: number = 10,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "api" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  if (!config.enableExternalApis) {
    return {
      source: "fallback",
      success: false,
      error: "External APIs disabled",
    };
  }

  try {
    const result = await searchSimilarCompounds(smiles, threshold, maxResults);
    return {
      source: "api",
      success: result.success,
      data: result,
      error: result.error,
    };
  } catch (error) {
    return {
      source: "fallback",
      success: false,
      error: `Similarity search error: ${error instanceof Error ? error.message : String(error)}`,
    };
  }
}

/**
 * Route patent status check to external API
 */
export async function routePatentCheck(
  cid: number,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "api" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  if (!config.enableExternalApis) {
    return {
      source: "fallback",
      success: false,
      error: "External APIs disabled",
    };
  }

  try {
    const result = await checkPatentStatus(cid);
    return {
      source: "api",
      success: result.success,
      data: result,
      error: result.error,
    };
  } catch (error) {
    return {
      source: "fallback",
      success: false,
      error: `Patent check error: ${error instanceof Error ? error.message : String(error)}`,
    };
  }
}

/**
 * Route compound properties to external API
 */
export async function routeCompoundProperties(
  smiles: string,
  config: ServiceConfig = defaultConfig
): Promise<{
  source: "api" | "fallback";
  success: boolean;
  data?: unknown;
  error?: string;
}> {
  if (!config.enableExternalApis) {
    return {
      source: "fallback",
      success: false,
      error: "External APIs disabled",
    };
  }

  try {
    const result = await getCompoundProperties(smiles);
    return {
      source: "api",
      success: result.success,
      data: result,
      error: result.error,
    };
  } catch (error) {
    return {
      source: "fallback",
      success: false,
      error: `Property fetch error: ${error instanceof Error ? error.message : String(error)}`,
    };
  }
}

/**
 * Get service health status
 */
export async function getServiceStatus(
  config: ServiceConfig = defaultConfig
): Promise<{
  pythonService: boolean;
  externalApis: boolean;
  timestamp: string;
}> {
  const pythonHealthy = await checkPythonServiceHealth(config);

  return {
    pythonService: pythonHealthy,
    externalApis: config.enableExternalApis,
    timestamp: new Date().toISOString(),
  };
}
