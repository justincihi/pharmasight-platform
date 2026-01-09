import { eq, or, like, desc, gte } from "drizzle-orm";
import { drizzle } from "drizzle-orm/mysql2";
import { InsertUser, users, analogDiscoveries, testResults, notifications, chatMessages, bookmarks, InsertAnalogDiscovery, InsertTestResult, InsertNotification, InsertChatMessage, InsertBookmark } from "../drizzle/schema";
import { ENV } from './_core/env';

let _db: ReturnType<typeof drizzle> | null = null;

// Lazily create the drizzle instance so local tooling can run without a DB.
export async function getDb() {
  if (!_db && process.env.DATABASE_URL) {
    try {
      _db = drizzle(process.env.DATABASE_URL);
    } catch (error) {
      console.warn("[Database] Failed to connect:", error);
      _db = null;
    }
  }
  return _db;
}

export async function upsertUser(user: InsertUser): Promise<void> {
  if (!user.openId) {
    throw new Error("User openId is required for upsert");
  }

  const db = await getDb();
  if (!db) {
    console.warn("[Database] Cannot upsert user: database not available");
    return;
  }

  try {
    const values: InsertUser = {
      openId: user.openId,
    };
    const updateSet: Record<string, unknown> = {};

    const textFields = ["name", "email", "loginMethod"] as const;
    type TextField = (typeof textFields)[number];

    const assignNullable = (field: TextField) => {
      const value = user[field];
      if (value === undefined) return;
      const normalized = value ?? null;
      values[field] = normalized;
      updateSet[field] = normalized;
    };

    textFields.forEach(assignNullable);

    if (user.lastSignedIn !== undefined) {
      values.lastSignedIn = user.lastSignedIn;
      updateSet.lastSignedIn = user.lastSignedIn;
    }
    if (user.role !== undefined) {
      values.role = user.role;
      updateSet.role = user.role;
    } else if (user.openId === ENV.ownerOpenId) {
      values.role = 'admin';
      updateSet.role = 'admin';
    }

    if (!values.lastSignedIn) {
      values.lastSignedIn = new Date();
    }

    if (Object.keys(updateSet).length === 0) {
      updateSet.lastSignedIn = new Date();
    }

    await db.insert(users).values(values).onDuplicateKeyUpdate({
      set: updateSet,
    });
  } catch (error) {
    console.error("[Database] Failed to upsert user:", error);
    throw error;
  }
}

export async function getUserByOpenId(openId: string) {
  const db = await getDb();
  if (!db) {
    console.warn("[Database] Cannot get user: database not available");
    return undefined;
  }

  const result = await db.select().from(users).where(eq(users.openId, openId)).limit(1);

  return result.length > 0 ? result[0] : undefined;
}

/**
 * Analog Discovery Queries
 */

/**
 * Get test results for specific analog IDs
 */
export async function getTestResultsByAnalogIds(analogIds: number[]) {
  const db = await getDb();
  if (!db || analogIds.length === 0) return [];

  try {
    const results = await (db as any)
      .select()
      .from(testResults)
      .where((col: any) => analogIds.includes(col.analogId));
    return results;
  } catch (error) {
    console.error("Error fetching test results by analog IDs:", error);
    return [];
  }
}

export async function getAnalogDiscoveries(
  limit: number = 50,
  offset: number = 0,
  filters?: { patentStatus?: string; minConfidence?: number }
) {
  const db = await getDb();
  if (!db) return [];

  try {
    // Cast to any to bypass strict type checking
    const query = (db as any).select().from(analogDiscoveries);
    
    let finalQuery = query;
    if (filters?.minConfidence !== undefined) {
      const minConf = filters.minConfidence;
      finalQuery = finalQuery.where((col: any) => col.confidenceScore >= minConf);
    }
    if (filters?.patentStatus) {
      const status = filters.patentStatus;
      finalQuery = finalQuery.where((col: any) => col.patentStatus === status);
    }

    return await finalQuery
      .orderBy((col: any) => col.discoveredAt)
      .limit(limit)
      .offset(offset);
  } catch (error) {
    console.error("Error fetching analogs:", error);
    return [];
  }
}

export async function getAnalogById(id: number) {
  const db = await getDb();
  if (!db) return null;

  const result = await (db as any)
    .select()
    .from(analogDiscoveries)
    .where((col: any) => col.id === id)
    .limit(1);

  return result.length > 0 ? result[0] : null;
}

export async function searchAnalogs(
  searchQuery: string,
  limit: number = 50
) {
  const db = await getDb();
  if (!db) return [];

  // Search by compound name or SMILES
  return await (db as any)
    .select()
    .from(analogDiscoveries)
    .where(
      (col: any) =>
        col.compoundName.like(`%${searchQuery}%`) ||
        col.smiles.like(`%${searchQuery}%`)
    )
    .limit(limit);
}

export async function createAnalogDiscovery(
  data: InsertAnalogDiscovery
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db.insert(analogDiscoveries).values(data);
  
  // Return the created record
  const result = await (db as any)
    .select()
    .from(analogDiscoveries)
    .where((col: any) => col.compoundId === data.compoundId)
    .limit(1);

  return result[0];
}

export async function getAnalyticsStats() {
  const db = await getDb();
  if (!db) return null;

  const allAnalogs = await db.select().from(analogDiscoveries);
  
  const totalDiscovered = allAnalogs.length;
  const highConfidence = allAnalogs.filter(
    (a) => a.confidenceScore >= 85
  ).length;
  const patentFree = allAnalogs.filter(
    (a) => a.patentStatus === "patent-free"
  ).length;
  const patented = allAnalogs.filter(
    (a) => a.patentStatus === "patented"
  ).length;

  return {
    totalDiscovered,
    highConfidenceCount: highConfidence,
    highConfidencePercentage: totalDiscovered > 0 
      ? Math.round((highConfidence / totalDiscovered) * 100)
      : 0,
    patentFreeCount: patentFree,
    patentedCount: patented,
    patentFreePercentage: totalDiscovered > 0
      ? Math.round((patentFree / totalDiscovered) * 100)
      : 0,
    averageConfidence: totalDiscovered > 0
      ? Math.round(
          allAnalogs.reduce((sum, a) => sum + a.confidenceScore, 0) /
            totalDiscovered
        )
      : 0,
  };
}

export async function getDiscoveryTimeline(days: number = 30) {
  const db = await getDb();
  if (!db) return [];

  const cutoffDate = new Date();
  cutoffDate.setDate(cutoffDate.getDate() - days);

  return await (db as any)
    .select()
    .from(analogDiscoveries)
    .where((col: any) => col.discoveredAt >= cutoffDate)
    .orderBy((col: any) => col.discoveredAt);
}

/**
 * Test Results Queries
 */
export async function createTestResult(
  data: InsertTestResult
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db.insert(testResults).values(data);
  
  const result = await (db as any)
    .select()
    .from(testResults)
    .orderBy((col: any) => col.id)
    .limit(1);

  return result[0];
}

export async function getTestResults(analogId: number) {
  const db = await getDb();
  if (!db) return [];

  return await (db as any)
    .select()
    .from(testResults)
    .where((col: any) => col.analogId === analogId);
}

export async function getAllTestResults() {
  const db = await getDb();
  if (!db) return [];

  return await (db as any)
    .select()
    .from(testResults)
    .orderBy((col: any) => col.createdAt);
}

/**
 * Notification Queries
 */
export async function createNotification(
  data: InsertNotification
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db.insert(notifications).values(data);
}

export async function getAdminNotifications(
  userId: number,
  limit: number = 20
) {
  const db = await getDb();
  if (!db) return [];

  // Get all notifications for admin users (not filtered by userId)
  return await db
    .select()
    .from(notifications)
    .orderBy(desc(notifications.createdAt))
    .limit(limit);
}

export async function getUnreadNotificationCount(
  userId: number
): Promise<number> {
  const db = await getDb();
  if (!db) return 0;

  const result = await db
    .select()
    .from(notifications)
    .where(eq(notifications.isRead, 0));
  
  return result.length;
}

export async function markNotificationAsRead(
  notificationId: number
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db
    .update(notifications)
    .set({ isRead: 1 })
    .where(eq(notifications.id, notificationId));
}

export async function markAllNotificationsAsRead(
  userId: number
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db
    .update(notifications)
    .set({ isRead: 1 })
    .where(eq(notifications.isRead, 0));
}

export async function getNewNotifications(
  userId: number,
  since: string | null
) {
  const db = await getDb();
  if (!db) return { notifications: [], lastChecked: new Date().toISOString() };

  let query = db.select().from(notifications);
  
  if (since) {
    const sinceDate = new Date(since);
    query = query.where(gte(notifications.createdAt, sinceDate)) as any;
  }
  
  const results = await query.orderBy(desc(notifications.createdAt)).limit(50);
  
  return {
    notifications: results,
    lastChecked: new Date().toISOString(),
  };
}

/**
 * Chat Message Queries
 */
export async function saveChatMessage(
  data: InsertChatMessage
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db.insert(chatMessages).values(data);
}

export async function getChatHistory(
  userId: number,
  limit: number = 50
) {
  const db = await getDb();
  if (!db) return [];

  return await (db as any)
    .select()
    .from(chatMessages)
    .where((col: any) => col.userId === userId)
    .orderBy((col: any) => col.createdAt)
    .limit(limit);
}


/**
 * Bookmark Queries
 */
export async function createBookmark(data: InsertBookmark) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  const result = await db.insert(bookmarks).values(data);
  return { id: result[0].insertId };
}

export async function getUserBookmarks(userId: number, limit: number = 50) {
  const db = await getDb();
  if (!db) return [];

  return await db
    .select()
    .from(bookmarks)
    .where(eq(bookmarks.userId, userId))
    .orderBy(desc(bookmarks.createdAt))
    .limit(limit);
}

export async function getBookmarkById(bookmarkId: number) {
  const db = await getDb();
  if (!db) return null;

  const results = await db
    .select()
    .from(bookmarks)
    .where(eq(bookmarks.id, bookmarkId))
    .limit(1);

  return results[0] || null;
}

export async function updateBookmark(
  bookmarkId: number,
  data: Partial<{ title: string; notes: string; category: string }>
) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db
    .update(bookmarks)
    .set(data as any)
    .where(eq(bookmarks.id, bookmarkId));
}

export async function deleteBookmark(bookmarkId: number) {
  const db = await getDb();
  if (!db) throw new Error("Database not available");

  await db.delete(bookmarks).where(eq(bookmarks.id, bookmarkId));
}

export async function isAnalogBookmarked(userId: number, analogId: number): Promise<boolean> {
  const db = await getDb();
  if (!db) return false;

  const results = await db
    .select()
    .from(bookmarks)
    .where(eq(bookmarks.userId, userId))
    .limit(100);

  return results.some((b: any) => b.analogId === analogId);
}

export async function getBookmarkByAnalogId(userId: number, analogId: number) {
  const db = await getDb();
  if (!db) return null;

  const results = await db
    .select()
    .from(bookmarks)
    .where(eq(bookmarks.userId, userId))
    .limit(100);

  return results.find((b: any) => b.analogId === analogId) || null;
}

export async function getBookmarksByCategory(userId: number, category: string) {
  const db = await getDb();
  if (!db) return [];

  return await db
    .select()
    .from(bookmarks)
    .where(eq(bookmarks.userId, userId))
    .orderBy(desc(bookmarks.createdAt));
}
