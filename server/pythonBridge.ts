/**
 * Python Bridge for Cheminformatics - HTTP API Version
 * Calls Python backend microservices via HTTP instead of spawning processes
 */

interface PythonResult {
  success: boolean;
  data?: any;
  error?: string;
}

// Service URLs from environment variables
const BACKEND_PK_API_URL = process.env.BACKEND_PK_API_URL || 'http://localhost:8001';
const RESEARCH_ENGINE_URL = process.env.RESEARCH_ENGINE_URL || 'http://localhost:8006';
const ANALOG_SERVICE_URL = process.env.ANALOG_SERVICE_URL || 'http://localhost:8002';
const BIOTRANSFORMER_URL = process.env.BIOTRANSFORMER_URL || 'http://localhost:8007';

/**
 * Generic fetch wrapper with error handling
 */
async function callPythonService(
  url: string,
  method: 'GET' | 'POST' = 'POST',
  body?: any
): Promise<PythonResult> {
  try {
    const options: RequestInit = {
      method,
      headers: {
        'Content-Type': 'application/json',
      },
    };

    if (body && method === 'POST') {
      options.body = JSON.stringify(body);
    }

    const response = await fetch(url, options);

    if (!response.ok) {
      const errorText = await response.text();
      return {
        success: false,
        error: `HTTP ${response.status}: ${errorText}`,
      };
    }

    const data = await response.json();

    return {
      success: true,
      data,
    };
  } catch (error) {
    return {
      success: false,
      error: error instanceof Error ? error.message : String(error),
    };
  }
}

/**
 * ChEMBL Validation via Research Engine
 */
export async function validateWithChEMBL(smiles: string): Promise<PythonResult> {
  const url = `${RESEARCH_ENGINE_URL}/chembl/validate`;
  return callPythonService(url, 'POST', { smiles });
}

/**
 * ADMET Prediction via Compound Analysis Service
 */
export async function predictADMET(smiles: string): Promise<PythonResult> {
  const url = `${BACKEND_PK_API_URL}/admet/predict`;
  return callPythonService(url, 'POST', { smiles });
}

/**
 * Molecular Docking via Compound Analysis Service
 */
export async function runMolecularDocking(
  ligandSmiles: string,
  receptorPDB: string
): Promise<PythonResult> {
  const url = `${BACKEND_PK_API_URL}/docking/simulate`;
  return callPythonService(url, 'POST', {
    ligand_smiles: ligandSmiles,
    receptor_pdb: receptorPDB,
  });
}

/**
 * Toxicity Prediction via Compound Analysis Service
 */
export async function predictToxicity(smiles: string): Promise<PythonResult> {
  const url = `${BACKEND_PK_API_URL}/toxicity/predict`;
  return callPythonService(url, 'POST', { smiles });
}

/**
 * PK/PD Simulation via Backend PK API
 */
export async function simulatePKPD(
  smiles: string,
  dose: number,
  route: string
): Promise<PythonResult> {
  const url = `${BACKEND_PK_API_URL}/pk/simulate`;
  return callPythonService(url, 'POST', {
    drug_name: 'custom',
    dose_mg: dose,
    route,
    drug_params: {
      // Default parameters - can be enhanced with SMILES-based prediction
      cl: 10.0,
      vc: 50.0,
      ka: 0.8,
      f: 0.8,
    },
  });
}

/**
 * Generate Analogs using RDKit via Research Engine
 */
export async function generateAnalogs(
  parentSmiles: string,
  numAnalogs: number = 10
): Promise<PythonResult> {
  const url = `${RESEARCH_ENGINE_URL}/rdkit/generate-analogs`;
  return callPythonService(url, 'POST', {
    smiles: parentSmiles,
    num_analogs: numAnalogs,
  });
}

/**
 * Query PubChem Database
 */
export async function queryPubChem(compoundName: string): Promise<PythonResult> {
  const url = `${RESEARCH_ENGINE_URL}/pubchem/compound/${encodeURIComponent(compoundName)}`;
  return callPythonService(url, 'GET');
}

/**
 * Query ChEMBL Database
 */
export async function queryChEMBL(query: string): Promise<PythonResult> {
  const url = `${RESEARCH_ENGINE_URL}/chembl/search`;
  return callPythonService(url, 'POST', { query });
}

/**
 * Query External Database (Generic)
 */
export async function queryExternalDatabase(
  database: string,
  query: string
): Promise<PythonResult> {
  // Route to appropriate service based on database type
  switch (database.toLowerCase()) {
    case 'pubchem':
      return queryPubChem(query);
    case 'chembl':
      return queryChEMBL(query);
    default:
      return {
        success: false,
        error: `Unknown database: ${database}`,
      };
  }
}

/**
 * Get compound properties via RDKit
 */
export async function getCompoundProperties(smiles: string): Promise<PythonResult> {
  const url = `${RESEARCH_ENGINE_URL}/rdkit/properties`;
  return callPythonService(url, 'POST', { smiles });
}

/**
 * BioTransformer metabolite prediction
 */
export async function predictMetabolites(smiles: string): Promise<PythonResult> {
  const url = `${BIOTRANSFORMER_URL}/predict`;
  return callPythonService(url, 'POST', { smiles });
}
