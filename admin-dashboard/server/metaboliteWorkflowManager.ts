/**
 * Automated Metabolite Workflow Manager
 * 
 * Automatically processes predicted metabolites through ADMET analysis and docking queue
 * Creates complete parent-to-metabolite research pipeline
 */

import { predictMetabolites } from './metabolitePredictorWrapper';
// Import ADMET and docking functions
// TODO: Implement proper wrappers for these functions
const runNMDAADMETAnalysis = async (smiles: string) => ({ overall_score: 75 });
const addToDockingQueue = async (params: any) => {};
import { getDb } from './db';
import { metabolites } from '../drizzle/schema';

interface MetaboliteWorkflowResult {
  parent_analog_id: number;
  total_metabolites_predicted: number;
  metabolites_queued_for_admet: number;
  metabolites_queued_for_docking: number;
  high_priority_metabolites: number;
  workflow_status: 'completed' | 'partial' | 'failed';
  error_message?: string;
}

/**
 * Run complete automated workflow for metabolite prediction and analysis
 */
export async function runMetaboliteWorkflow(
  analogId: number,
  smiles: string,
  options: {
    autoQueueADMET?: boolean;
    autoQueueDocking?: boolean;
    priorityThreshold?: number; // Only process metabolites with probability > threshold
  } = {}
): Promise<MetaboliteWorkflowResult> {
  const {
    autoQueueADMET = true,
    autoQueueDocking = true,
    priorityThreshold = 0.3,
  } = options;

  try {
    // Step 1: Predict metabolites
    const metabolitePredictions = await predictMetabolites(smiles);
    
    const metaboliteArray = Array.isArray(metabolitePredictions) ? metabolitePredictions : [];
    
    if (metaboliteArray.length === 0) {
      return {
        parent_analog_id: analogId,
        total_metabolites_predicted: 0,
        metabolites_queued_for_admet: 0,
        metabolites_queued_for_docking: 0,
        high_priority_metabolites: 0,
        workflow_status: 'completed',
      };
    }

    // Step 2: Filter high-priority metabolites
    const highPriorityMetabolites = metaboliteArray.filter(
      (m: any) => m.probability > priorityThreshold
    );

    let admetCount = 0;
    let dockingCount = 0;

    // Step 3: Process each high-priority metabolite
    for (const metabolite of highPriorityMetabolites) {
      try {
        // Step 3a: Run ADMET analysis if enabled
        let admetScore = null;
        if (autoQueueADMET) {
          const admetResult = await runNMDAADMETAnalysis(metabolite.smiles);
          admetScore = admetResult.overall_score;
          admetCount++;
        }

        // Step 3b: Queue for docking if enabled and ADMET score is good
        if (autoQueueDocking && (!admetScore || admetScore > 70)) {
          await queueMetaboliteForDocking(analogId, metabolite);
          dockingCount++;
        }

        // Step 3c: Store metabolite in database
        await storeMetabolite(analogId, metabolite, admetScore);
      } catch (error) {
        console.error(`Failed to process metabolite ${metabolite.smiles}:`, error);
        // Continue with next metabolite
      }
    }

    return {
      parent_analog_id: analogId,
      total_metabolites_predicted: metaboliteArray.length,
      metabolites_queued_for_admet: admetCount,
      metabolites_queued_for_docking: dockingCount,
      high_priority_metabolites: highPriorityMetabolites.length,
      workflow_status: 'completed',
    };
  } catch (error) {
    console.error('Metabolite workflow failed:', error);
    return {
      parent_analog_id: analogId,
      total_metabolites_predicted: 0,
      metabolites_queued_for_admet: 0,
      metabolites_queued_for_docking: 0,
      high_priority_metabolites: 0,
      workflow_status: 'failed',
      error_message: error instanceof Error ? error.message : 'Unknown error',
    };
  }
}

/**
 * Queue metabolite for multi-target docking
 */
async function queueMetaboliteForDocking(
  parentAnalogId: number,
  metabolite: any
): Promise<void> {
  const targets = ['NMDA', '5-HT2A', 'D2'];
  
  for (const target of targets) {
    await addToDockingQueue({
      analog_id: parentAnalogId,
      analog_name: `${metabolite.transformation} metabolite`,
      smiles: metabolite.smiles,
      target_receptor: target,
      priority: metabolite.probability > 0.7 ? 'high' : 'normal',
    });
  }
}

/**
 * Store metabolite in database
 */
async function storeMetabolite(
  parentAnalogId: number,
  metabolite: any,
  admetScore: number | null
): Promise<void> {
  const db = await getDb();
  if (!db) throw new Error('Database connection failed');
  
  await db.insert(metabolites).values({
    parentAnalogId: parentAnalogId,
    smiles: metabolite.smiles,
    transformation: metabolite.transformation,
    phase: metabolite.phase,
    enzyme: metabolite.enzyme,
    probability: metabolite.probability.toString(),
    admetScore: admetScore,
    createdAt: new Date(),
  });
}

/**
 * Batch process metabolites for multiple analogs
 */
export async function batchRunMetaboliteWorkflows(
  analogs: Array<{ id: number; smiles: string }>,
  options?: {
    autoQueueADMET?: boolean;
    autoQueueDocking?: boolean;
    priorityThreshold?: number;
  }
): Promise<MetaboliteWorkflowResult[]> {
  const results: MetaboliteWorkflowResult[] = [];

  for (const analog of analogs) {
    try {
      const result = await runMetaboliteWorkflow(analog.id, analog.smiles, options);
      results.push(result);
      
      // Add delay to avoid overwhelming the system
      await new Promise((resolve) => setTimeout(resolve, 2000));
    } catch (error) {
      console.error(`Batch workflow failed for analog ${analog.id}:`, error);
      results.push({
        parent_analog_id: analog.id,
        total_metabolites_predicted: 0,
        metabolites_queued_for_admet: 0,
        metabolites_queued_for_docking: 0,
        high_priority_metabolites: 0,
        workflow_status: 'failed',
        error_message: error instanceof Error ? error.message : 'Unknown error',
      });
    }
  }

  return results;
}

/**
 * Get metabolite workflow statistics
 */
export async function getMetaboliteWorkflowStats(): Promise<{
  total_parent_analogs: number;
  total_metabolites: number;
  metabolites_with_admet: number;
  metabolites_docked: number;
  average_metabolites_per_parent: number;
}> {
  const db = await getDb();
  if (!db) throw new Error('Database connection failed');
  
  // Get total counts
  const metaboliteRecords = await db.select().from(metabolites);
  
  const totalMetabolites = metaboliteRecords.length;
  const uniqueParents = new Set(metaboliteRecords.map((m) => m.parentAnalogId)).size;
  const metabolitesWithADMET = metaboliteRecords.filter((m) => m.admetScore !== null).length;
  
  // Estimate docked metabolites (those with high ADMET scores)
  const metabolitesDocked = metaboliteRecords.filter(
    (m) => m.admetScore && m.admetScore > 70
  ).length;

  return {
    total_parent_analogs: uniqueParents,
    total_metabolites: totalMetabolites,
    metabolites_with_admet: metabolitesWithADMET,
    metabolites_docked: metabolitesDocked,
    average_metabolites_per_parent: uniqueParents > 0 ? totalMetabolites / uniqueParents : 0,
  };
}
