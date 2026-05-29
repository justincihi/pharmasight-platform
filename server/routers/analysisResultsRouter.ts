import { protectedProcedure, router } from '../_core/trpc';
import { z } from 'zod';

/**
 * Analysis Results Router
 * Handles persistent logging and retrieval of analysis results
 * Stores all docking, toxicity, ADMET, and PK/PD results with timestamps
 */
export const analysisResultsRouter = router({
  /**
   * Log an analysis result to the database
   */
  logResult: protectedProcedure
    .input(
      z.object({
        analogId: z.number(),
        analysisType: z.enum(['docking', 'toxicity', 'admet', 'pkpd']),
        smiles: z.string(),
        target: z.string().optional(),
        result: z.record(z.string(), z.unknown()),
        source: z.enum(['python', 'api', 'fallback']).optional(),
        executionTime: z.number().optional(), // milliseconds
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const { getDb } = await import('../db');
      const { analysisResults } = await import('../../drizzle/schema');
      const db = await getDb();

      if (!db) {
        throw new Error('Database connection failed');
      }

      try {
        await db.insert(analysisResults).values({
          analogId: input.analogId,
          analysisType: input.analysisType,
          smiles: input.smiles,
          target: input.target,
          result: input.result,
          source: input.source || 'python',
          executionTime: input.executionTime || 0,
          createdBy: ctx.user.id,
        });

        return {
          success: true,
          message: 'Result logged successfully',
        };
      } catch (error) {
        console.error('[Analysis Results] Failed to log result:', error);
        throw new Error(`Failed to log analysis result: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Get all results for an analog
   */
  getByAnalogId: protectedProcedure
    .input(
      z.object({
        analogId: z.number(),
        analysisType: z.enum(['docking', 'toxicity', 'admet', 'pkpd']).optional(),
        limit: z.number().optional().default(50),
        offset: z.number().optional().default(0),
      })
    )
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const { getDb } = await import('../db');
      const { analysisResults } = await import('../../drizzle/schema');
      const { eq, and } = await import('drizzle-orm');
      const db = await getDb();

      if (!db) {
        throw new Error('Database connection failed');
      }

      try {
        let query = db.select().from(analysisResults).where(eq(analysisResults.analogId, input.analogId));

        if (input.analysisType) {
          query = db
            .select()
            .from(analysisResults)
            .where(and(eq(analysisResults.analogId, input.analogId), eq(analysisResults.analysisType, input.analysisType)));
        }

        const results = await query.limit(input.limit).offset(input.offset);

        return {
          success: true,
          results: results.map((r) => ({
            ...r,
            result: typeof r.result === 'string' ? JSON.parse(r.result) : r.result,
          })),
        };
      } catch (error) {
        console.error('[Analysis Results] Failed to retrieve results:', error);
        throw new Error(`Failed to retrieve analysis results: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Get results by SMILES string
   */
  getBySMILES: protectedProcedure
    .input(
      z.object({
        smiles: z.string(),
        analysisType: z.enum(['docking', 'toxicity', 'admet', 'pkpd']).optional(),
        limit: z.number().optional().default(50),
      })
    )
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const { getDb } = await import('../db');
      const { analysisResults } = await import('../../drizzle/schema');
      const { eq, and } = await import('drizzle-orm');
      const db = await getDb();

      if (!db) {
        throw new Error('Database connection failed');
      }

      try {
        let query = db.select().from(analysisResults).where(eq(analysisResults.smiles, input.smiles));

        if (input.analysisType) {
          query = db
            .select()
            .from(analysisResults)
            .where(and(eq(analysisResults.smiles, input.smiles), eq(analysisResults.analysisType, input.analysisType)));
        }

        const results = await query.limit(input.limit);

        return {
          success: true,
          results: results.map((r) => ({
            ...r,
            result: typeof r.result === 'string' ? JSON.parse(r.result) : r.result,
          })),
        };
      } catch (error) {
        console.error('[Analysis Results] Failed to retrieve results by SMILES:', error);
        throw new Error(`Failed to retrieve analysis results: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Get analysis statistics
   */
  getStatistics: protectedProcedure
    .input(
      z.object({
        startDate: z.date().optional(),
        endDate: z.date().optional(),
        analysisType: z.enum(['docking', 'toxicity', 'admet', 'pkpd']).optional(),
      })
    )
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const { getDb } = await import('../db');
      const { analysisResults } = await import('../../drizzle/schema');
      const { sql } = await import('drizzle-orm');
      const db = await getDb();

      if (!db) {
        throw new Error('Database connection failed');
      }

      try {
        // Get total count
        const totalCount = await db.select({ count: sql`COUNT(*)` }).from(analysisResults);

        // Get count by analysis type
        const countByType = await db
          .select({
            type: analysisResults.analysisType,
            count: sql`COUNT(*)`,
          })
          .from(analysisResults)
          .groupBy(analysisResults.analysisType);

        // Get count by source
        const countBySource = await db
          .select({
            source: analysisResults.source,
            count: sql`COUNT(*)`,
          })
          .from(analysisResults)
          .groupBy(analysisResults.source);

        // Get average execution time
        const avgExecutionTime = await db
          .select({
            avgTime: sql`AVG(${analysisResults.executionTime})`,
          })
          .from(analysisResults);

        return {
          success: true,
          statistics: {
            totalResults: (totalCount[0]?.count as number) || 0,
            byType: countByType.map((row) => ({
              type: row.type,
              count: row.count as number,
            })),
            bySource: countBySource.map((row) => ({
              source: row.source,
              count: row.count as number,
            })),
            averageExecutionTime: (avgExecutionTime[0]?.avgTime as number) || 0,
          },
        };
      } catch (error) {
        console.error('[Analysis Results] Failed to get statistics:', error);
        throw new Error(`Failed to get analysis statistics: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Delete old results (cleanup)
   */
  deleteOlderThan: protectedProcedure
    .input(
      z.object({
        days: z.number().min(1).max(365),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const { getDb } = await import('../db');
      const { analysisResults } = await import('../../drizzle/schema');
      const { lt } = await import('drizzle-orm');
      const db = await getDb();

      if (!db) {
        throw new Error('Database connection failed');
      }

      try {
        const cutoffDate = new Date();
        cutoffDate.setDate(cutoffDate.getDate() - input.days);

        await db.delete(analysisResults).where(lt(analysisResults.createdAt, cutoffDate));

        return {
          success: true,
          message: `Deleted results older than ${input.days} days`,
        };
      } catch (error) {
        console.error('[Analysis Results] Failed to delete old results:', error);
        throw new Error(`Failed to delete old results: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),
});
