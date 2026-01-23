import { notifyOwner } from "./_core/notification";

/**
 * Send notification when new analog is discovered
 */
export async function notifyAnalogDiscovery(analogData: {
  compoundId: string;
  parentCompound: string;
  confidence: number;
  similarity: number;
}) {
  const title = `🔬 New Analog Discovered: ${analogData.compoundId}`;
  const content = `
Parent Compound: ${analogData.parentCompound}
Confidence: ${analogData.confidence}%
Similarity: ${analogData.similarity}%

A new high-confidence analog has been discovered by the autonomous research engine.
View details in the PharmaSight™ Admin Dashboard.
  `.trim();

  return await notifyOwner({ title, content });
}

/**
 * Send notification when synthesis routes are generated
 */
export async function notifySynthesisRoutes(routeData: {
  compoundId: string;
  compoundName: string;
  routeCount: number;
  bestFeasibility: number;
  lowestCost: number;
}) {
  const title = `⚗️ Synthesis Routes Generated: ${routeData.compoundName}`;
  const content = `
Compound: ${routeData.compoundId}
Routes Generated: ${routeData.routeCount}
Best Feasibility: ${routeData.bestFeasibility}/100
Lowest Cost: $${routeData.lowestCost}

AI-powered retrosynthesis analysis complete. View detailed routes in the dashboard.
  `.trim();

  return await notifyOwner({ title, content });
}

/**
 * Send notification when batch analysis completes
 */
export async function notifyBatchComplete(batchData: {
  batchId: string;
  compoundCount: number;
  successCount: number;
  duration: string;
}) {
  const title = `📊 Batch Analysis Complete: ${batchData.batchId}`;
  const content = `
Compounds Analyzed: ${batchData.compoundCount}
Successful: ${batchData.successCount}
Duration: ${batchData.duration}

Batch cheminformatics analysis has completed. View results in the dashboard.
  `.trim();

  return await notifyOwner({ title, content });
}

/**
 * Send notification when scheduler runs
 */
export async function notifySchedulerRun(runData: {
  discoveryCount: number;
  timestamp: string;
  status: "success" | "error";
  message?: string;
}) {
  const title = runData.status === "success" 
    ? `✅ Scheduler Run Complete` 
    : `❌ Scheduler Run Failed`;
  
  const content = `
Timestamp: ${runData.timestamp}
Discoveries Imported: ${runData.discoveryCount}
Status: ${runData.status.toUpperCase()}
${runData.message ? `\nMessage: ${runData.message}` : ''}

Autonomous research engine has completed its scheduled run.
  `.trim();

  return await notifyOwner({ title, content });
}
