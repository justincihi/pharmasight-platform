/**
 * Discovery Audit Router - tRPC procedures for discovery management
 */

import { router, protectedProcedure } from "../_core/trpc";
import { z } from "zod";

export const discoveryAuditRouter = router({
  /**
   * Get recent discoveries from the last 60 days
   */
  getRecentDiscoveries: protectedProcedure
    .input(
      z.object({
        days: z.number().default(60),
        limit: z.number().default(50),
        offset: z.number().default(0),
      })
    )
    .query(async ({ input, ctx }) => {
      try {
        // Only admins can view discoveries
        if (ctx.user?.role !== "admin") {
          throw new Error("Unauthorized: Admin access required");
        }

        const { getDb } = await import("../db");
        const { analogDiscoveries } = await import("../../drizzle/schema");
        const { gte } = await import("drizzle-orm");

        const db = await getDb();
        if (!db) throw new Error("Database connection failed");

        // Calculate date range
        const daysAgo = new Date();
        daysAgo.setDate(daysAgo.getDate() - input.days);

        // Build query
        const discoveries = await db
          .select()
          .from(analogDiscoveries)
          .where(gte(analogDiscoveries.createdAt, daysAgo))
          .limit(input.limit)
          .offset(input.offset);

        return {
          success: true,
          discoveries: discoveries.map((d) => ({
            id: d.id,
            compoundName: d.compoundName,
            smiles: d.smiles,
            confidenceScore: d.confidenceScore,
            patentStatus: d.patentStatus,
            createdAt: d.createdAt,
            source: "automated_research",
          })),
          total: discoveries.length,
        };
      } catch (error) {
        console.error("Error fetching recent discoveries:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
          discoveries: [],
          total: 0,
        };
      }
    }),

  /**
   * Get discovery statistics
   */
  getDiscoveryStats: protectedProcedure
    .input(z.object({ days: z.number().default(60) }))
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin") {
          throw new Error("Unauthorized: Admin access required");
        }

        const { getDb } = await import("../db");
        const { analogDiscoveries } = await import("../../drizzle/schema");
        const { gte } = await import("drizzle-orm");

        const db = await getDb();
        if (!db) throw new Error("Database connection failed");

        const daysAgo = new Date();
        daysAgo.setDate(daysAgo.getDate() - input.days);

        const discoveries = await db
          .select()
          .from(analogDiscoveries)
          .where(gte(analogDiscoveries.createdAt, daysAgo));

        const stats = {
          totalDiscoveries: discoveries.length,
          completeCompounds: discoveries.filter(
            (d) => d.confidenceScore >= 85
          ).length,
          partialScaffolds: discoveries.filter(
            (d) => d.confidenceScore < 85
          ).length,
          patentFree: discoveries.filter((d) => d.patentStatus === "patent-free")
            .length,
          avgConfidence:
            discoveries.length > 0
              ? Math.round(
                  discoveries.reduce((sum, d) => sum + d.confidenceScore, 0) /
                    discoveries.length
                )
              : 0,
        };

        return {
          success: true,
          stats,
        };
      } catch (error) {
        console.error("Error fetching discovery stats:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
          stats: null,
        };
      }
    }),

  /**
   * Save discovery to master list
   */
  saveToMasterList: protectedProcedure
    .input(
      z.object({
        discoveryId: z.number(),
        compoundName: z.string(),
        smiles: z.string(),
        confidenceScore: z.number(),
        patentStatus: z.enum(["patent-free", "patented", "unknown"]),
      })
    )
    .mutation(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin") {
          throw new Error("Unauthorized: Admin access required");
        }

        const { getDb } = await import("../db");
        const { analogDiscoveries } = await import("../../drizzle/schema");
        const { eq } = await import("drizzle-orm");

        const db = await getDb();
        if (!db) throw new Error("Database connection failed");

        // Update discovery status
        await db
          .update(analogDiscoveries)
          .set({
            approvalStatus: "approved",
            approvedBy: ctx.user.openId,
            approvedAt: new Date(),
          })
          .where(eq(analogDiscoveries.id, input.discoveryId));

        return {
          success: true,
          message: `${input.compoundName} added to master list`,
        };
      } catch (error) {
        console.error("Error saving to master list:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * Save discovery to scaffold library
   */
  saveToScaffoldLibrary: protectedProcedure
    .input(
      z.object({
        discoveryId: z.number(),
        scaffoldName: z.string(),
        smiles: z.string(),
        scaffoldType: z.enum([
          "fragment",
          "core",
          "side-chain",
          "linker",
          "other",
        ]),
        potentialApplications: z.string().optional(),
      })
    )
    .mutation(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin") {
          throw new Error("Unauthorized: Admin access required");
        }

        const { getDb } = await import("../db");
        const { analogDiscoveries } = await import("../../drizzle/schema");
        const { eq } = await import("drizzle-orm");

        const db = await getDb();
        if (!db) throw new Error("Database connection failed");

        // Update discovery with scaffold info
        await db
          .update(analogDiscoveries)
          .set({
            approvalStatus: "approved",
            approvedBy: ctx.user.openId,
            approvedAt: new Date(),
          })
          .where(eq(analogDiscoveries.id, input.discoveryId));

        return {
          success: true,
          message: `${input.scaffoldName} saved to scaffold library`,
        };
      } catch (error) {
        console.error("Error saving to scaffold library:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * Get scaffold library
   */
  getScaffoldLibrary: protectedProcedure
    .input(
      z.object({
        limit: z.number().default(50),
        offset: z.number().default(0),
      })
    )
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin") {
          throw new Error("Unauthorized: Admin access required");
        }

        const { getDb } = await import("../db");
        const { analogDiscoveries } = await import("../../drizzle/schema");
        const { eq } = await import("drizzle-orm");

        const db = await getDb();
        if (!db) throw new Error("Database connection failed");

        const scaffolds = await db
          .select()
          .from(analogDiscoveries)
          .limit(input.limit)
          .offset(input.offset);

        return {
          success: true,
          scaffolds: scaffolds.map((s) => ({
            id: s.id,
            scaffoldName: s.compoundName,
            smiles: s.smiles,
            createdAt: s.createdAt,
          })),
          total: scaffolds.length,
        };
      } catch (error) {
        console.error("Error fetching scaffold library:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
          scaffolds: [],
          total: 0,
        };
      }
    }),
});
