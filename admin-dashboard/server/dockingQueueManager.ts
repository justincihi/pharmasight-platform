import { getDb } from './db';
import { dockingQueue, analogDiscoveries } from '../drizzle/schema';
import { eq, and, or } from 'drizzle-orm';
import { runMolecularDocking } from './molecularDockingWrapper';

const DOCKING_TARGETS = ['NMDA', '5-HT2A', 'D2'] as const;
type DockingTarget = typeof DOCKING_TARGETS[number];

/**
 * Add analog to docking queue for all targets
 */
export async function queueAnalogForDocking(analogId: number, priority: number = 5) {
  const db = await getDb();
  if (!db) throw new Error('Database not available');
  
  // Add entry for each target
  const entries = DOCKING_TARGETS.map(target => ({
    analogId,
    target,
    priority,
    status: 'pending' as const,
  }));

  await db.insert(dockingQueue).values(entries);
  
  console.log(`[Docking Queue] Added analog ${analogId} to queue for ${DOCKING_TARGETS.length} targets`);
}

/**
 * Process next pending docking job
 */
export async function processNextDockingJob(): Promise<boolean> {
  const db = await getDb();
  if (!db) throw new Error('Database not available');
  
  // Find next pending job (highest priority first)
  const jobs = await db
    .select()
    .from(dockingQueue)
    .where(eq(dockingQueue.status, 'pending'))
    .orderBy(dockingQueue.priority)
    .limit(1);

  if (jobs.length === 0) {
    return false; // No jobs to process
  }

  const job = jobs[0];
  
  try {
    // Mark as running
    await db
      .update(dockingQueue)
      .set({
        status: 'running',
        startedAt: new Date(),
      })
      .where(eq(dockingQueue.id, job.id));

    // Get analog SMILES
    const analogs = await db
      .select()
      .from(analogDiscoveries)
      .where(eq(analogDiscoveries.id, job.analogId))
      .limit(1);

    if (analogs.length === 0) {
      throw new Error(`Analog ${job.analogId} not found`);
    }

    const analog = analogs[0];

    // Run docking
    console.log(`[Docking Queue] Processing job ${job.id}: Analog ${job.analogId} → ${job.target}`);
    const result = await runMolecularDocking(analog.smiles, analog.compoundId, job.target);

    if (!result.success) {
      throw new Error(result.error || 'Docking failed');
    }

    // Update job with results
    await db
      .update(dockingQueue)
      .set({
        status: 'completed',
        bindingAffinity: result.binding_affinity?.toString(),
        dockingScore: Math.round((10 - Math.abs(result.binding_affinity || 0)) * 10), // Normalize to 0-100
        completedAt: new Date(),
      })
      .where(eq(dockingQueue.id, job.id));

    // Update analog with best docking score if this is better
    if (result.binding_affinity) {
      const currentBest = analog.bindingAffinity 
        ? parseFloat(analog.bindingAffinity) 
        : 0;
      
      if (!currentBest || result.binding_affinity < currentBest) {
        await db
          .update(analogDiscoveries)
          .set({
            bindingAffinity: result.binding_affinity.toString(),
            dockingScore: Math.round((10 - Math.abs(result.binding_affinity)) * 10),
            dockingTarget: job.target,
          })
          .where(eq(analogDiscoveries.id, job.analogId));
      }
    }

    console.log(`[Docking Queue] Job ${job.id} completed: ${result.binding_affinity} kcal/mol`);
    return true;
  } catch (error: any) {
    console.error(`[Docking Queue] Job ${job.id} failed:`, error);
    
    // Mark as failed
    await db
      .update(dockingQueue)
      .set({
        status: 'failed',
        errorMessage: error.message,
        completedAt: new Date(),
      })
      .where(eq(dockingQueue.id, job.id));

    return true; // Job was processed (even though it failed)
  }
}

/**
 * Start docking queue worker (processes jobs continuously)
 */
export async function startDockingQueueWorker(intervalMs: number = 30000) {
  console.log('[Docking Queue] Worker started');
  
  const processJobs = async () => {
    try {
      let processed = true;
      while (processed) {
        processed = await processNextDockingJob();
        if (processed) {
          // Small delay between jobs
          await new Promise(resolve => setTimeout(resolve, 1000));
        }
      }
    } catch (error) {
      console.error('[Docking Queue] Worker error:', error);
    }
  };

  // Process immediately
  await processJobs();

  // Then check periodically
  setInterval(processJobs, intervalMs);
}

/**
 * Get queue status
 */
export async function getDockingQueueStatus() {
  const db = await getDb();
  if (!db) throw new Error('Database not available');
  
  const pending = await db
    .select()
    .from(dockingQueue)
    .where(eq(dockingQueue.status, 'pending'));

  const running = await db
    .select()
    .from(dockingQueue)
    .where(eq(dockingQueue.status, 'running'));

  const completed = await db
    .select()
    .from(dockingQueue)
    .where(eq(dockingQueue.status, 'completed'));

  const failed = await db
    .select()
    .from(dockingQueue)
    .where(eq(dockingQueue.status, 'failed'));

  return {
    pending: pending.length,
    running: running.length,
    completed: completed.length,
    failed: failed.length,
    total: pending.length + running.length + completed.length + failed.length,
  };
}

/**
 * Get docking results for an analog
 */
export async function getAnalogDockingResults(analogId: number) {
  const db = await getDb();
  if (!db) throw new Error('Database not available');
  
  return await db
    .select()
    .from(dockingQueue)
    .where(and(
      eq(dockingQueue.analogId, analogId),
      eq(dockingQueue.status, 'completed')
    ));
}
