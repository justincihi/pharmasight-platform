import { startScheduler, getSchedulerStatus, setCronSchedule } from "./autonomousScheduler";
import { getSetting } from "./db";

const DEFAULT_CRON = process.env.SCHEDULER_CRON || "0 9 * * *";

/**
 * Initialize the autonomous research scheduler on server startup.
 * Loads the persisted cron schedule from the database (app_settings key "scheduler_cron")
 * so the schedule survives server restarts.
 */
export async function initializeScheduler() {
  console.log("[Init] Initializing autonomous research scheduler...");

  // Load persisted cron schedule from DB (falls back to env/default)
  let cronSchedule = DEFAULT_CRON;
  try {
    const saved = await getSetting("scheduler_cron");
    if (saved) {
      cronSchedule = saved;
      console.log(`[Init] Loaded persisted cron schedule from DB: ${cronSchedule}`);
    }
  } catch (err) {
    console.warn("[Init] Could not load cron schedule from DB, using default:", err);
  }

  // Apply the loaded schedule to the in-memory state
  setCronSchedule(cronSchedule);

  startScheduler({
    enabled: process.env.SCHEDULER_ENABLED !== "false",
    cronSchedule,
    minConfidenceForNotification: parseInt(process.env.MIN_CONFIDENCE || "85"),
  });

  const status = getSchedulerStatus();
  if (status.running) {
    console.log(`[Init] ✅ Scheduler running - Next run: ${status.nextRun}`);
  } else {
    console.log("[Init] ⚠️  Scheduler is disabled");
  }
}
