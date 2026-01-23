import { startScheduler, getSchedulerStatus } from "./autonomousScheduler";

/**
 * Initialize the autonomous research scheduler on server startup
 */
export function initializeScheduler() {
  console.log("[Init] Initializing autonomous research scheduler...");
  
  // Start the scheduler with default configuration
  startScheduler({
    enabled: process.env.SCHEDULER_ENABLED !== "false",
    cronSchedule: process.env.SCHEDULER_CRON || "0 9 * * *", // 9 AM daily
    minConfidenceForNotification: parseInt(process.env.MIN_CONFIDENCE || "85"),
  });

  const status = getSchedulerStatus();
  
  if (status.running) {
    console.log(`[Init] ✅ Scheduler running - Next run: ${status.nextRun}`);
  } else {
    console.log("[Init] ⚠️  Scheduler is disabled");
  }
}
