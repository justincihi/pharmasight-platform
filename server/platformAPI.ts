/**
 * REST API endpoints for PharmaSight Platform integration
 * These endpoints allow the Python research engine to communicate with the dashboard
 */

import type { Express, Request, Response } from "express";
import { getDb } from "./db";
import { analogDiscoveries } from "../drizzle/schema";
import { eq, desc } from "drizzle-orm";
import { syncToMasterFile } from "./masterFileSync";

/**
 * Register platform API routes
 */
export function registerPlatformAPI(app: Express) {
  // Health check endpoint
  app.get("/api/platform/health", (_req: Request, res: Response) => {
    res.json({ status: "ok", service: "pharmasight-admin-dashboard", timestamp: new Date().toISOString() });
  });

  // Import analog discoveries from research engine
  app.post("/api/platform/discoveries/import", async (req: Request, res: Response) => {
    try {
      const { discoveries, apiKey } = req.body;

      // Simple API key validation (in production, use proper auth)
      if (apiKey !== process.env.PLATFORM_API_KEY) {
        return res.status(401).json({ error: "Unauthorized" });
      }

      if (!Array.isArray(discoveries)) {
        return res.status(400).json({ error: "Invalid request: discoveries must be an array" });
      }

      const db = await getDb();
      if (!db) {
        return res.status(500).json({ error: "Database not available" });
      }

      const imported: any[] = [];
      const errors: any[] = [];

      for (const discovery of discoveries) {
        try {
          // Check if analog already exists
          const existing = await db
            .select()
            .from(analogDiscoveries)
            .where(eq(analogDiscoveries.compoundId, discovery.compoundId))
            .limit(1);

          if (existing.length === 0) {
            // Insert new analog
            await db.insert(analogDiscoveries).values({
              compoundId: discovery.compoundId,
              compoundName: discovery.compoundName || discovery.compoundId,
              smiles: discovery.smiles,
              parentCompound: discovery.parentCompound,
              mechanismOfAction: discovery.mechanismOfAction || "",
              keyDifferences: discovery.keyDifferences || "",
              confidenceScore: discovery.confidence || 0,
              similarityScore: typeof discovery.similarity === 'number' && discovery.similarity <= 1 ? Math.round(discovery.similarity * 100) : (discovery.similarity || 0),
              safetyScore: discovery.safetyScore || 0,
              efficacyScore: discovery.efficacyScore || 0,
              drugLikenessScore: discovery.drugLikenessScore || 100,
              patentStatus: discovery.patentStatus || "unknown",
              marketValue: discovery.marketValue || "0",
              discoveredBy: "platform-api",
              discoveryMethod: discovery.discoveryMethod || "external-import",
            });

            imported.push(discovery.compoundId);
          } else {
            console.log(`[Platform API] Analog ${discovery.compoundId} already exists, skipping`);
          }
        } catch (error) {
          console.error(`[Platform API] Error importing ${discovery.compoundId}:`, error);
          errors.push({ compoundId: discovery.compoundId, error: String(error) });
        }
      }

      // Sync to master file after import
      await syncToMasterFile();
      
      res.json({
        success: true,
        imported: imported.length,
        skipped: discoveries.length - imported.length - errors.length,
        errors: errors.length,
        details: { imported, errors },
      });
    } catch (error) {
      console.error("[Platform API] Error in import endpoint:", error);
      res.status(500).json({ error: "Internal server error", message: String(error) });
    }
  });

  // Get recent analog discoveries
  app.get("/api/platform/discoveries/recent", async (req: Request, res: Response) => {
    try {
      const { apiKey, limit = 10 } = req.query;

      if (apiKey !== process.env.PLATFORM_API_KEY) {
        return res.status(401).json({ error: "Unauthorized" });
      }

      const db = await getDb();
      if (!db) {
        return res.status(500).json({ error: "Database not available" });
      }

      const recentAnalogs = await db
        .select()
        .from(analogDiscoveries)
        .orderBy(desc(analogDiscoveries.discoveredAt))
        .limit(Number(limit));

      res.json({
        success: true,
        count: recentAnalogs.length,
        discoveries: recentAnalogs,
      });
    } catch (error) {
      console.error("[Platform API] Error in recent discoveries endpoint:", error);
      res.status(500).json({ error: "Internal server error", message: String(error) });
    }
  });

  // Get analog by compound ID
  app.get("/api/platform/analogs/:compoundId", async (req: Request, res: Response) => {
    try {
      const { apiKey } = req.query;
      const { compoundId } = req.params;

      if (apiKey !== process.env.PLATFORM_API_KEY) {
        return res.status(401).json({ error: "Unauthorized" });
      }

      const db = await getDb();
      if (!db) {
        return res.status(500).json({ error: "Database not available" });
      }

      const analog = await db
        .select()
        .from(analogDiscoveries)
        .where(eq(analogDiscoveries.compoundId, compoundId))
        .limit(1);

      if (analog.length === 0) {
        return res.status(404).json({ error: "Analog not found" });
      }

      res.json({
        success: true,
        analog: analog[0],
      });
    } catch (error) {
      console.error("[Platform API] Error in get analog endpoint:", error);
      res.status(500).json({ error: "Internal server error", message: String(error) });
    }
  });

  // Update analog data
  app.put("/api/platform/analogs/:compoundId", async (req: Request, res: Response) => {
    try {
      const { apiKey, ...updateData } = req.body;
      const { compoundId } = req.params;

      if (apiKey !== process.env.PLATFORM_API_KEY) {
        return res.status(401).json({ error: "Unauthorized" });
      }

      const db = await getDb();
      if (!db) {
        return res.status(500).json({ error: "Database not available" });
      }

      await db
        .update(analogDiscoveries)
        .set(updateData)
        .where(eq(analogDiscoveries.compoundId, compoundId));

      res.json({
        success: true,
        message: `Analog ${compoundId} updated successfully`,
      });
    } catch (error) {
      console.error("[Platform API] Error in update analog endpoint:", error);
      res.status(500).json({ error: "Internal server error", message: String(error) });
    }
  });

  console.log("[Platform API] REST API endpoints registered");
}
