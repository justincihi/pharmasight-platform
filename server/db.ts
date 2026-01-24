import { eq, or, like } from "drizzle-orm";
import { drizzle } from "drizzle-orm/mysql2";
import { InsertUser, users, analogDiscoveries, testResults, notifications, chatMessages, InsertAnalogDiscovery, InsertTestResult, InsertNotification, InsertChatMessage } from "../drizzle/schema";
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

  return await (db as any)
    .select()
    .from(notifications)
    .where((col: any) => col.userId === userId)
    .orderBy((col: any) => col.createdAt)
    .limit(limit);
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

