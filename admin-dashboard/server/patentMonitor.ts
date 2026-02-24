import { getDb } from './db';
import { analogDiscoveries, notifications } from '../drizzle/schema';
import { eq } from 'drizzle-orm';
import { searchPatents, checkExpiringPatents } from './patentSearchIntegration';

/**
 * Monitor high-value analogs for patent status changes
 */
export async function monitorPatentStatus() {
  const db = await getDb();
  if (!db) throw new Error('Database not available');

  console.log('[Patent Monitor] Starting patent status check...');

  // Get all high-confidence, patent-free analogs
  const analogs = await db
    .select()
    .from(analogDiscoveries)
    .where(eq(analogDiscoveries.patentStatus, 'patent-free'));

  let checked = 0;
  let alerts = 0;

  for (const analog of analogs) {
    try {
      // Search for patents
      const result = await searchPatents(analog.smiles, analog.compoundName);

      // Check for expiring patents that might create opportunities
      const expirationAlerts = await checkExpiringPatents(analog.id, result.patents);

      if (expirationAlerts.length > 0) {
        // Create notifications for expiring patents
        for (const alert of expirationAlerts) {
          await db.insert(notifications).values({
            userId: 1, // Admin user
            analogId: analog.id,
            title: `Patent Expiration Alert: ${analog.compoundName}`,
            message: alert,
            notificationType: 'patent-alert',
            isRead: 0,
          });
          alerts++;
        }
      }

      // Update patent status if it has changed
      if (result.found && result.patents.length > 0) {
        const patentNumbers = result.patents.map(p => p.patentNumber);
        
        await db
          .update(analogDiscoveries)
          .set({
            patentNumbers: JSON.stringify(patentNumbers),
            patentStatus: result.freedomToOperate ? 'patent-free' : 'patented',
          })
          .where(eq(analogDiscoveries.id, analog.id));
      }

      checked++;
    } catch (error) {
      console.error(`[Patent Monitor] Error checking analog ${analog.id}:`, error);
    }

    // Rate limiting - don't overwhelm APIs
    await new Promise(resolve => setTimeout(resolve, 2000));
  }

  console.log(`[Patent Monitor] Checked ${checked} analogs, created ${alerts} alerts`);
  
  return {
    checked,
    alerts,
  };
}

/**
 * Start patent monitoring scheduler (runs daily)
 */
export function startPatentMonitor() {
  console.log('[Patent Monitor] Scheduler started - will run daily at 2 AM');

  // Run immediately on startup
  setTimeout(() => {
    monitorPatentStatus().catch(console.error);
  }, 5000);

  // Then run daily at 2 AM
  const runDaily = () => {
    const now = new Date();
    const next2AM = new Date(now);
    next2AM.setHours(2, 0, 0, 0);
    
    if (next2AM <= now) {
      next2AM.setDate(next2AM.getDate() + 1);
    }

    const msUntil2AM = next2AM.getTime() - now.getTime();

    setTimeout(() => {
      monitorPatentStatus().catch(console.error);
      runDaily(); // Schedule next run
    }, msUntil2AM);
  };

  runDaily();
}
