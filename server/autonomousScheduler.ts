import { CronJob } from "cron";
import { importAutonomousDiscoveries } from "./importDiscoveries";
import { runAutonomousResearch } from "./runAutonomousResearch";
import { getDb } from "./db";
import { analogDiscoveries } from "../drizzle/schema";
import { desc, gte, and, eq } from "drizzle-orm";
import { notifyOwner } from "./_core/notification";

/**
 * Autonomous Research Scheduler
 * Runs daily to import new discoveries and notify admins of high-confidence analogs
 */

interface SchedulerConfig {
  cronSchedule: string; // Default: "0 9 * * *" (9 AM daily)
  minConfidenceForNotification: number; // Default: 85
  enabled: boolean;
}

const defaultConfig: SchedulerConfig = {
  cronSchedule: process.env.SCHEDULER_CRON || "0 9 * * *", // 9 AM daily
  minConfidenceForNotification: 85,
  enabled: process.env.SCHEDULER_ENABLED !== "false",
};

let schedulerJob: CronJob | null = null;

/**
 * Main scheduler task that runs on schedule
 */
async function runScheduledTask() {
  console.log(`[Scheduler] Running autonomous research import at ${new Date().toISOString()}`);

  try {
    // Step 1: Run autonomous research engine to discover new analogs
    console.log("[Scheduler] Running autonomous research engine...");
    const researchResult = await runAutonomousResearch();
    
    if (!researchResult.success) {
      console.error("[Scheduler] Research engine failed:", researchResult.error);
      return;
    }
    
    console.log(`[Scheduler] Research engine discovered ${researchResult.discoveries?.length || 0} new compounds`);
    
    // Step 2: Import discoveries into database
    const importResult = await importAutonomousDiscoveries();
    
    if (!importResult.success) {
      console.error("[Scheduler] Failed to import discoveries:", importResult.error);
      return;
    }

    console.log(`[Scheduler] Successfully imported ${importResult.importedCount} discoveries`);

    // Step 2: Check for high-confidence discoveries from the last 24 hours
    const db = await getDb();
    if (!db) {
      console.error("[Scheduler] Database not available");
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

    // Step 3: Send notification if high-confidence discoveries found
    if (highConfidenceDiscoveries.length > 0) {
      const notificationContent = generateNotificationContent(highConfidenceDiscoveries);
      
      const notified = await notifyOwner({
        title: `🔬 ${highConfidenceDiscoveries.length} High-Confidence Analog(s) Discovered`,
        content: notificationContent,
      });

      if (notified) {
        console.log(`[Scheduler] Notification sent for ${highConfidenceDiscoveries.length} high-confidence discoveries`);
      } else {
        console.warn("[Scheduler] Failed to send notification");
      }

      // Also log to database notifications table
      await logNotifications(db, highConfidenceDiscoveries);
    } else {
      console.log("[Scheduler] No new high-confidence discoveries in the last 24 hours");
    }

  } catch (error) {
    console.error("[Scheduler] Error during scheduled task:", error);
  }
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
 * Log notifications to database
 */
async function logNotifications(db: any, discoveries: any[]) {
  const { notifications } = await import("../drizzle/schema");
  
  for (const discovery of discoveries) {
    try {
      await db.insert(notifications).values({
        analogId: discovery.id,
        type: "high_confidence_discovery",
        message: `New high-confidence analog discovered: ${discovery.compoundName} (${discovery.confidenceScore}% confidence)`,
        sentAt: new Date(),
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
      runScheduledTask,
      null, // onComplete
      true, // start immediately
      "America/New_York" // timezone
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
 * Run the scheduler task immediately (for testing)
 */
export async function runSchedulerNow() {
  console.log("[Scheduler] Running scheduler task immediately...");
  await runScheduledTask();
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
