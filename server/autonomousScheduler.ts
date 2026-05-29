import { CronJob } from "cron";
import { importAutonomousDiscoveries } from "./importDiscoveries";
import { runAutonomousResearch } from "./runAutonomousResearch";
import { getDb } from "./db";
import { analogDiscoveries, researchRuns } from "../drizzle/schema";
import { desc, gte, and, eq } from "drizzle-orm";
import { notifyOwner } from "./_core/notification";
import { randomUUID } from "crypto";

/**
 * Autonomous Research Scheduler
 * Runs daily to import new discoveries and notify admins of high-confidence analogs
 */

interface SchedulerConfig {
  cronSchedule: string;
  minConfidenceForNotification: number;
  enabled: boolean;
}

const defaultConfig: SchedulerConfig = {
  cronSchedule: process.env.SCHEDULER_CRON || "0 9 * * *",
  minConfidenceForNotification: 85,
  enabled: process.env.SCHEDULER_ENABLED !== "false",
};

let schedulerJob: CronJob | null = null;

type ProgressLevel = "info" | "success" | "warning" | "error";

/**
 * Main scheduler task that runs on schedule
 */
async function runScheduledTask(triggeredBy: "manual" | "scheduled" = "scheduled", userId?: number) {
  const startTime = Date.now();
  const runId = randomUUID();
  
  console.log(`[Scheduler] Starting research run ${runId} at ${new Date().toISOString()}`);

  const db = await getDb();
  const progressLog: Array<{ timestamp: string; message: string; level: ProgressLevel }> = [];

  const log = (message: string, level: ProgressLevel = "info") => {
    const entry = { timestamp: new Date().toISOString(), message, level };
    progressLog.push(entry);
    console.log(`[Scheduler] [${level.toUpperCase()}] ${message}`);
  };

  // Create run record
  let runDbId: number | null = null;
  if (db) {
    try {
      const result = await db.insert(researchRuns).values({
        runId,
        triggeredBy,
        triggeredByUserId: userId,
        status: "running",
        startedAt: new Date(),
      });
      runDbId = (result as any).insertId ?? null;
    } catch (e) {
      console.error("[Scheduler] Failed to create run record:", e);
    }
  }

  const updateRun = async (updates: Partial<{
    status: "running" | "completed" | "failed";
    discoveriesCount: number;
    highConfidenceCount: number;
    articlesScanned: number;
    goalsUsed: string[];
    topDiscoveries: any[];
    articlesLog: any[];
    progressLog: any[];
    errorMessage: string;
    completedAt: Date;
    durationMs: number;
  }>) => {
    if (!db || !runDbId) return;
    try {
      await db.update(researchRuns).set(updates).where(eq(researchRuns.id, runDbId));
    } catch (e) {
      console.error("[Scheduler] Failed to update run record:", e);
    }
  };

  try {
    log("Starting autonomous research engine...", "info");
    await updateRun({ progressLog: [...progressLog] });

    // Step 1: Run autonomous research engine
    const researchResult = await runAutonomousResearch();
    
    if (!researchResult.success) {
      log(`Research engine failed: ${researchResult.error}`, "error");
      await updateRun({
        status: "failed",
        errorMessage: researchResult.error,
        progressLog: [...progressLog],
        completedAt: new Date(),
        durationMs: Date.now() - startTime,
      });
      return;
    }
    
    const discoveryCount = researchResult.discoveries?.length || 0;
    const articleCount = researchResult.articlesScanned || 0;
    const goalsUsed = researchResult.goalsUsed || [];
    const articlesLog = researchResult.articles || [];

    log(`Research engine discovered ${discoveryCount} new compounds from ${articleCount} articles`, "success");
    await updateRun({
      discoveriesCount: discoveryCount,
      articlesScanned: articleCount,
      goalsUsed,
      articlesLog,
      progressLog: [...progressLog],
    });

    // Step 2: Import discoveries into database
    log("Importing discoveries into database...", "info");
    const importResult = await importAutonomousDiscoveries();
    
    if (!importResult.success) {
      log(`Failed to import discoveries: ${importResult.error}`, "warning");
    } else {
      log(`Successfully imported ${importResult.importedCount} discoveries`, "success");
    }

    await updateRun({ progressLog: [...progressLog] });

    // Step 3: Check for high-confidence discoveries from the last 24 hours
    if (!db) {
      log("Database not available for notification check", "warning");
      return;
    }

    const oneDayAgo = new Date(Date.now() - 24 * 60 * 60 * 1000);
    
    const highConfidenceDiscoveries = await db
      .select()
      .from(analogDiscoveries)
      .where(
        and(
          gte(analogDiscoveries.discoveredAt, oneDayAgo),
          gte(analogDiscoveries.confidenceScore, defaultConfig.minConfidenceForNotification)
        )
      )
      .orderBy(desc(analogDiscoveries.confidenceScore))
      .limit(10);

    const topDiscoveries = highConfidenceDiscoveries.map(d => ({
      compoundId: d.compoundId,
      compoundName: d.compoundName,
      parentCompound: d.parentCompound,
      confidenceScore: d.confidenceScore,
      safetyScore: d.safetyScore,
      efficacyScore: d.efficacyScore,
      patentStatus: d.patentStatus,
      smiles: d.smiles,
    }));

    log(`Found ${highConfidenceDiscoveries.length} high-confidence discoveries (≥${defaultConfig.minConfidenceForNotification}%)`, 
      highConfidenceDiscoveries.length > 0 ? "success" : "info");

    // Step 4: Send notification if high-confidence discoveries found
    if (highConfidenceDiscoveries.length > 0) {
      const notificationContent = generateNotificationContent(highConfidenceDiscoveries);
      
      const notified = await notifyOwner({
        title: `🔬 ${highConfidenceDiscoveries.length} High-Confidence Analog(s) Discovered`,
        content: notificationContent,
      });

      if (notified) {
        log(`Notification sent for ${highConfidenceDiscoveries.length} high-confidence discoveries`, "success");
      } else {
        log("Notification service temporarily unavailable", "warning");
      }

      await logNotifications(db, highConfidenceDiscoveries);
    }

    // Final update
    const duration = Date.now() - startTime;
    log(`Research run completed in ${(duration / 1000).toFixed(1)}s`, "success");
    
    await updateRun({
      status: "completed",
      highConfidenceCount: highConfidenceDiscoveries.length,
      topDiscoveries,
      progressLog: [...progressLog],
      completedAt: new Date(),
      durationMs: duration,
    });

  } catch (error) {
    const errMsg = error instanceof Error ? error.message : String(error);
    log(`Unexpected error: ${errMsg}`, "error");
    console.error("[Scheduler] Error during scheduled task:", error);
    
    await updateRun({
      status: "failed",
      errorMessage: errMsg,
      progressLog: [...progressLog],
      completedAt: new Date(),
      durationMs: Date.now() - startTime,
    });
  }

  return runId;
}

/**
 * Generate notification content from discoveries
 */
function generateNotificationContent(discoveries: any[]): string {
  let content = `PharmaSight™ has discovered ${discoveries.length} new high-confidence analog(s) in the last 24 hours:\n\n`;

  discoveries.forEach((discovery, index) => {
    const patentStatus = discovery.patentStatus === "patent_free" ? "✅ Patent-Free" : "⚠️ Patent Opportunity";
    
    content += `${index + 1}. **${discovery.compoundName}**\n`;
    content += `   - Confidence: ${discovery.confidenceScore}%\n`;
    content += `   - Parent: ${discovery.parentCompound}\n`;
    content += `   - Safety: ${discovery.safetyScore}/100\n`;
    content += `   - Efficacy: ${discovery.efficacyScore}/100\n`;
    content += `   - Patent Status: ${patentStatus}\n`;
    content += `   - Market Value: $${discovery.marketValue}\n\n`;
  });

  content += `\nView full details in the PharmaSight™ Admin Dashboard.`;
  
  return content;
}

/**
 * Log notifications to database for real-time display
 */
async function logNotifications(db: any, discoveries: any[]) {
  const { notifications } = await import("../drizzle/schema");
  
  for (const discovery of discoveries) {
    try {
      await db.insert(notifications).values({
        userId: 1,
        analogId: discovery.id,
        title: `🔬 High-Confidence Discovery: ${discovery.compoundName}`,
        message: `New high-confidence analog discovered: ${discovery.compoundName} (${discovery.confidenceScore}% confidence). Parent: ${discovery.parentCompound}. Safety: ${discovery.safetyScore}/100. Efficacy: ${discovery.efficacyScore}/100.`,
        notificationType: "high-confidence" as const,
        isRead: 0,
      });
    } catch (error) {
      console.error(`[Scheduler] Failed to log notification for ${discovery.compoundName}:`, error);
    }
  }
}

/**
 * Start the autonomous research scheduler
 */
export function startScheduler(config: Partial<SchedulerConfig> = {}) {
  const finalConfig = { ...defaultConfig, ...config };

  if (!finalConfig.enabled) {
    console.log("[Scheduler] Autonomous research scheduler is disabled");
    return;
  }

  if (schedulerJob) {
    console.log("[Scheduler] Scheduler already running");
    return;
  }

  try {
    schedulerJob = new CronJob(
      finalConfig.cronSchedule,
      () => { runScheduledTask("scheduled").catch(console.error); },
      null,
      true,
      "America/New_York"
    );

    console.log(`[Scheduler] Autonomous research scheduler started with cron: ${finalConfig.cronSchedule}`);
    console.log(`[Scheduler] Next run: ${schedulerJob.nextDate().toISO()}`);
  } catch (error) {
    console.error("[Scheduler] Failed to start scheduler:", error);
  }
}

/**
 * Stop the scheduler
 */
export function stopScheduler() {
  if (schedulerJob) {
    schedulerJob.stop();
    schedulerJob = null;
    console.log("[Scheduler] Scheduler stopped");
  }
}

/**
 * Run the scheduler task immediately (for testing / manual trigger)
 */
export async function runSchedulerNow(userId?: number) {
  console.log("[Scheduler] Running scheduler task immediately...");
  return await runScheduledTask("manual", userId);
}

/**
 * Get scheduler status
 */
export function getSchedulerStatus() {
  return {
    running: schedulerJob !== null,
    nextRun: schedulerJob?.nextDate().toISO() || null,
    config: defaultConfig,
  };
}
