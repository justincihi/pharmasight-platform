/**
 * Python Microservices Gateway
 * Routes tRPC calls to persistent Python services via HTTP
 * Implements hybrid architecture: Node.js frontend + Python backend services
 * 
 * This solves the "spawn ENOENT" issue by:
 * 1. Running Python services on persistent compute (separate from Cloud Run)
 * 2. Calling them via HTTP instead of child_process.spawn
 * 3. Falling back to mock responses if services are unavailable
 */

import axios, { AxiosError } from 'axios';

// Service configuration - points to persistent Python services
const PYTHON_SERVICES = {
  docking: process.env.PYTHON_DOCKING_SERVICE_URL || 'http://localhost:5001',
  admet: process.env.PYTHON_ADMET_SERVICE_URL || 'http://localhost:5002',
  toxicity: process.env.PYTHON_TOXICITY_SERVICE_URL || 'http://localhost:5003',
  metabolites: process.env.PYTHON_METABOLITES_SERVICE_URL || 'http://localhost:5004',
  leadOptimization: process.env.PYTHON_LEAD_OPT_SERVICE_URL || 'http://localhost:5005',
  bionemo: process.env.PYTHON_BIONEMO_SERVICE_URL || 'http://localhost:5006',
};

// Service health cache (5 minute TTL)
const healthCache = new Map<string, { status: boolean; timestamp: number }>();
const HEALTH_CHECK_TTL = 5 * 60 * 1000;

/**
 * Check if a Python service is available
 */
async function isServiceAvailable(serviceName: keyof typeof PYTHON_SERVICES): Promise<boolean> {
  const cached = healthCache.get(serviceName);
  if (cached && Date.now() - cached.timestamp < HEALTH_CHECK_TTL) {
    return cached.status;
  }

  try {
    const serviceUrl = PYTHON_SERVICES[serviceName];
    const response = await axios.get(`${serviceUrl}/health`, { timeout: 2000 });
    const isAvailable = response.status === 200;
    healthCache.set(serviceName, { status: isAvailable, timestamp: Date.now() });
    return isAvailable;
  } catch (error) {
    healthCache.set(serviceName, { status: false, timestamp: Date.now() });
    return false;
  }
}

/**
 * Call Python docking service
 */
export async function callDockingService(
  smiles: string,
  targetName: string,
  receptorPdb?: string
): Promise<any> {
  const isAvailable = await isServiceAvailable('docking');

  if (!isAvailable) {
    console.warn('[Docking Gateway] Service unavailable. Using mock response.');
    return getMockDockingResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.docking}/dock`, {
      smiles,
      targetName,
      receptorPdb,
    });
    return response.data;
  } catch (error) {
    console.error('[Docking Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockDockingResult();
  }
}

/**
 * Call Python ADMET service
 */
export async function callADMETService(smiles: string): Promise<any> {
  const isAvailable = await isServiceAvailable('admet');

  if (!isAvailable) {
    console.warn('[ADMET Gateway] Service unavailable. Using mock response.');
    return getMockADMETResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.admet}/predict`, {
      smiles,
    });
    return response.data;
  } catch (error) {
    console.error('[ADMET Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockADMETResult();
  }
}

/**
 * Call Python toxicity service
 */
export async function callToxicityService(smiles: string): Promise<any> {
  const isAvailable = await isServiceAvailable('toxicity');

  if (!isAvailable) {
    console.warn('[Toxicity Gateway] Service unavailable. Using mock response.');
    return getMockToxicityResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.toxicity}/predict`, {
      smiles,
    });
    return response.data;
  } catch (error) {
    console.error('[Toxicity Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockToxicityResult();
  }
}

/**
 * Call Python metabolite prediction service (Biotransformer)
 */
export async function callMetaboliteService(smiles: string, phase: 'all' | 'phase1' | 'phase2' | 'phase3' = 'all'): Promise<any> {
  const isAvailable = await isServiceAvailable('metabolites');

  if (!isAvailable) {
    console.warn('[Metabolite Gateway] Service unavailable. Using mock response.');
    return getMockMetaboliteResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.metabolites}/predict`, {
      smiles,
      phase,
    });
    return response.data;
  } catch (error) {
    console.error('[Metabolite Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockMetaboliteResult();
  }
}

/**
 * Call Python lead optimization service (dragonfly_gen)
 */
export async function callLeadOptimizationService(
  parentSmiles: string,
  numAnalogs: number = 20,
  constraints?: Record<string, any>
): Promise<any> {
  const isAvailable = await isServiceAvailable('leadOptimization');

  if (!isAvailable) {
    console.warn('[Lead Optimization Gateway] Service unavailable. Using mock response.');
    return getMockLeadOptimizationResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.leadOptimization}/optimize`, {
      parentSmiles,
      numAnalogs,
      constraints,
    });
    return response.data;
  } catch (error) {
    console.error('[Lead Optimization Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockLeadOptimizationResult();
  }
}

/**
 * Call Python BioNemo service for protein language models
 */
export async function callBioNemoService(
  proteinSequence: string,
  task: 'embedding' | 'classification' | 'prediction' = 'embedding'
): Promise<any> {
  const isAvailable = await isServiceAvailable('bionemo');

  if (!isAvailable) {
    console.warn('[BioNemo Gateway] Service unavailable. Using mock response.');
    return getMockBioNemoResult();
  }

  try {
    const response = await axios.post(`${PYTHON_SERVICES.bionemo}/process`, {
      proteinSequence,
      task,
    });
    return response.data;
  } catch (error) {
    console.error('[BioNemo Gateway] Error:', error instanceof AxiosError ? error.message : error);
    return getMockBioNemoResult();
  }
}

/**
 * Get service status for admin dashboard
 */
export async function getServiceStatus(): Promise<Record<string, boolean>> {
  const status: Record<string, boolean> = {};

  for (const [serviceName] of Object.entries(PYTHON_SERVICES)) {
    status[serviceName] = await isServiceAvailable(serviceName as keyof typeof PYTHON_SERVICES);
  }

  return status;
}

// ============================================================================
// Mock Responses (used when services are unavailable)
// ============================================================================

function getMockDockingResult() {
  return {
    success: true,
    binding_affinity: -6.5 + Math.random() * 3,
    rmsd: 1.2 + Math.random() * 0.5,
    pose: 'mock_pose_data',
    isDemo: true,
    message: 'Demo mode: Python docking service unavailable',
  };
}

function getMockADMETResult() {
  return {
    smiles: '',
    admet_results: {
      absorption: 0.75,
      distribution: 0.68,
      metabolism: 0.82,
      excretion: 0.71,
      toxicity: 0.65,
      molecular_weight: 350,
      logp: 2.5,
      hbd: 2,
      hba: 3,
      tpsa: 45,
    },
    isDemo: true,
    message: 'Demo mode: Python ADMET service unavailable',
  };
}

function getMockToxicityResult() {
  return {
    smiles: '',
    toxicity_score: 0.3 + Math.random() * 0.3,
    toxicity_class: 'Low',
    predictions: {
      hepatotoxicity: 0.15,
      nephrotoxicity: 0.25,
      cardiotoxicity: 0.1,
      neurotoxicity: 0.35,
    },
    isDemo: true,
    message: 'Demo mode: Python toxicity service unavailable',
  };
}

function getMockMetaboliteResult() {
  return {
    total_metabolites: 5,
    metabolites: [
      {
        smiles: 'CC(C)Cc1ccc(cc1)C(C)C(O)=O',
        name: 'Metabolite 1',
        phase: 'phase1',
        probability: 0.8,
      },
      {
        smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        name: 'Metabolite 2',
        phase: 'phase2',
        probability: 0.6,
      },
    ],
    isDemo: true,
    message: 'Demo mode: Python metabolite service unavailable',
  };
}

function getMockLeadOptimizationResult() {
  return {
    total_analogs_generated: 20,
    valid_analogs: 18,
    top_leads: [
      {
        smiles: 'CC(C)Cc1ccc(cc1)C(C)C(O)=O',
        name: 'Lead 1',
        rank: 1,
        predicted_potency: 85,
        predicted_selectivity: 90,
        predicted_admet_score: 80,
        overall_score: 8.5,
        rationale: 'Improved potency with maintained selectivity',
      },
      {
        smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        name: 'Lead 2',
        rank: 2,
        predicted_potency: 80,
        predicted_selectivity: 85,
        predicted_admet_score: 75,
        overall_score: 8.0,
        rationale: 'Good balance of properties',
      },
    ],
    isDemo: true,
    message: 'Demo mode: Python lead optimization service unavailable',
  };
}

function getMockBioNemoResult() {
  return {
    task: 'embedding',
    embedding_dim: 1280,
    embedding: Array(1280).fill(0.1),
    confidence: 0.85,
    isDemo: true,
    message: 'Demo mode: Python BioNemo service unavailable',
  };
}
