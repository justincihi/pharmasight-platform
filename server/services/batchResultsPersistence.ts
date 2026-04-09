import { getDb } from "../db";
import { batchDockingJobs, batchDockingResults, analogDiscoveries } from "../../drizzle/schema";
import { eq } from "drizzle-orm";
import { randomUUID } from "crypto";

export interface BatchJobConfig {
  userId: string;
  jobName: string;
  targetName: string;
  parametersId?: number;
  totalCompounds: number;
}

export interface BatchResultData {
  analogId: string;
  bindingAffinity: number;
  numPoses: number;
  dockingScore: number;
  topPoses: any[];
  status: "completed" | "failed";
  errorMessage?: string;
}

/**
 * Create a new batch docking job
 */
export async function createBatchJob(config: BatchJobConfig) {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const jobId = randomUUID();
    
    await db.insert(batchDockingJobs).values({
      id: jobId,
      userId: config.userId,
      jobName: config.jobName,
      status: "pending",
      totalCompounds: config.totalCompounds,
      completedCompounds: 0,
      failedCompounds: 0,
      targetName: config.targetName,
      parametersId: config.parametersId,
      createdAt: new Date(),
      updatedAt: new Date(),
    });

    console.log(`[Batch Persistence] Created job: ${jobId}`);
    return jobId;
  } catch (error) {
    console.error("[Batch Persistence] Error creating job:", error);
    throw error;
  }
}

/**
 * Update batch job status
 */
export async function updateBatchJobStatus(
  jobId: string,
  status: "pending" | "running" | "completed" | "failed" | "cancelled",
  completedCount?: number,
  failedCount?: number
) {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const updateData: any = {
      status,
      updatedAt: new Date(),
    };

    if (completedCount !== undefined) updateData.completedCompounds = completedCount;
    if (failedCount !== undefined) updateData.failedCompounds = failedCount;
    if (status === "running") updateData.startedAt = new Date();
    if (status === "completed" || status === "failed") updateData.completedAt = new Date();

    await db
      .update(batchDockingJobs)
      .set(updateData)
      .where(eq(batchDockingJobs.id, jobId));

    console.log(`[Batch Persistence] Updated job ${jobId} to ${status}`);
  } catch (error) {
    console.error("[Batch Persistence] Error updating job status:", error);
    throw error;
  }
}

/**
 * Store individual batch docking result
 */
export async function storeBatchResult(
  jobId: string,
  analogId: string,
  result: BatchResultData
) {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    // Find analog ID from compound ID
    const analog = await db
      .select()
      .from(analogDiscoveries)
      .where(eq(analogDiscoveries.compoundId, analogId))
      .limit(1);

    if (!analog || analog.length === 0) {
      throw new Error(`Analog not found: ${analogId}`);
    }

    await db.insert(batchDockingResults).values({
      jobId,
      analogId: analog[0].id,
      status: result.status,
      bindingAffinity: result.bindingAffinity.toString(),
      dockingScore: result.dockingScore,
      numPoses: result.numPoses,
      topPoses: JSON.stringify(result.topPoses),
      errorMessage: result.errorMessage,
      completedAt: result.status === "completed" ? new Date() : undefined,
      createdAt: new Date(),
      updatedAt: new Date(),
    });

    console.log(`[Batch Persistence] Stored result for analog: ${analogId}`);
  } catch (error) {
    console.error("[Batch Persistence] Error storing result:", error);
    throw error;
  }
}

/**
 * Get batch job details
 */
export async function getBatchJobDetails(jobId: string) {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, jobId))
      .limit(1);

    if (!job || job.length === 0) {
      throw new Error(`Job not found: ${jobId}`);
    }

    const results = await db
      .select()
      .from(batchDockingResults)
      .where(eq(batchDockingResults.jobId, jobId));

    return {
      job: job[0],
      results,
    };
  } catch (error) {
    console.error("[Batch Persistence] Error fetching job details:", error);
    throw error;
  }
}

/**
 * Calculate batch job statistics
 */
export function calculateBatchStatistics(results: any[]) {
  if (results.length === 0) {
    return {
      totalResults: 0,
      successCount: 0,
      failureCount: 0,
      averageAffinity: 0,
      bestAffinity: 0,
      worstAffinity: 0,
      averageScore: 0,
    };
  }

  const successful = results.filter((r) => r.status === "completed");
  const affinities = successful
    .map((r) => parseFloat(r.bindingAffinity || "0"))
    .filter((a) => !isNaN(a));

  return {
    totalResults: results.length,
    successCount: successful.length,
    failureCount: results.filter((r) => r.status === "failed").length,
    averageAffinity: affinities.length > 0 ? affinities.reduce((a, b) => a + b, 0) / affinities.length : 0,
    bestAffinity: affinities.length > 0 ? Math.min(...affinities) : 0,
    worstAffinity: affinities.length > 0 ? Math.max(...affinities) : 0,
    averageScore: successful.length > 0
      ? successful.reduce((sum, r) => sum + (r.dockingScore || 0), 0) / successful.length
      : 0,
  };
}

/**
 * Update batch job with summary statistics
 */
export async function updateBatchJobSummary(jobId: string, statistics: any) {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    await db
      .update(batchDockingJobs)
      .set({
        resultsSummary: JSON.stringify(statistics),
        updatedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, jobId));

    console.log(`[Batch Persistence] Updated summary for job: ${jobId}`);
  } catch (error) {
    console.error("[Batch Persistence] Error updating summary:", error);
    throw error;
  }
}
