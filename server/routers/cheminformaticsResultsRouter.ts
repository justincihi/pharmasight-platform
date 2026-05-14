import { z } from 'zod';
import { protectedProcedure, router } from '../_core/trpc';
import { getDb } from '../db';
import { cheminformaticsResults } from '../../drizzle/schema';
import { eq, desc, and } from 'drizzle-orm';
import { TRPCError } from '@trpc/server';

/**
 * Cheminformatics Results Router
 * Handles persistence and retrieval of pipeline execution results
 */
export const cheminformaticsResultsRouter = router({
  /**
   * Save cheminformatics pipeline results
   */
  saveResults: protectedProcedure
    .input(
      z.object({
        inputSmiles: z.string().min(1),
        canonicalSmiles: z.string().optional(),
        workflow: z.enum(['similarity', 'brics', 'validate', 'full_pipeline']),
        threshold: z.number().min(0).max(1).optional().default(0.70),
        maxHits: z.number().min(1).max(100).optional().default(25),
        results: z.object({
          success: z.boolean(),
          error: z.string().optional(),
          data: z.any().optional(),
        }),
        executionTime: z.number().optional(),
        notes: z.string().optional(),
      })
    )
    .mutation(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      await db.insert(cheminformaticsResults).values({
        userId: ctx.user.id,
        inputSmiles: input.inputSmiles,
        canonicalSmiles: input.canonicalSmiles,
        workflow: input.workflow,
        threshold: input.threshold.toString() as any,
        maxHits: input.maxHits,
        results: input.results as any,
        executionTime: input.executionTime,
        status: input.results.success ? 'completed' : 'failed',
        notes: input.notes,
      });

      return {
        success: true,
        message: 'Results saved successfully',
      };
    }),

  /**
   * Get a specific result by ID
   */
  getResult: protectedProcedure
    .input(z.object({ id: z.number() }))
    .query(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      const result = await db
        .select()
        .from(cheminformaticsResults)
        .where(
          and(
            eq(cheminformaticsResults.id, input.id),
            eq(cheminformaticsResults.userId, ctx.user.id)
          )
        )
        .limit(1);

      if (!result.length) {
        throw new TRPCError({
          code: 'NOT_FOUND',
          message: 'Result not found',
        });
      }

      return result[0];
    }),

  /**
   * List user's cheminformatics results with pagination
   */
  listResults: protectedProcedure
    .input(
      z.object({
        limit: z.number().min(1).max(100).optional().default(20),
        offset: z.number().min(0).optional().default(0),
        workflow: z.enum(['similarity', 'brics', 'validate', 'full_pipeline']).optional(),
        status: z.enum(['pending', 'running', 'completed', 'failed']).optional(),
      })
    )
    .query(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      const whereConditions = [eq(cheminformaticsResults.userId, ctx.user.id)];

      if (input.workflow) {
        whereConditions.push(eq(cheminformaticsResults.workflow, input.workflow));
      }

      if (input.status) {
        whereConditions.push(eq(cheminformaticsResults.status, input.status));
      }

      const results = await db
        .select()
        .from(cheminformaticsResults)
        .where(and(...whereConditions))
        .orderBy(desc(cheminformaticsResults.createdAt))
        .limit(input.limit)
        .offset(input.offset);

      const total = await db
        .select()
        .from(cheminformaticsResults)
        .where(and(...whereConditions));

      return {
        results,
        total: total.length,
        limit: input.limit,
        offset: input.offset,
      };
    }),

  /**
   * Search results by SMILES or workflow
   */
  searchResults: protectedProcedure
    .input(
      z.object({
        query: z.string().optional(),
        workflow: z.enum(['similarity', 'brics', 'validate', 'full_pipeline']).optional(),
        limit: z.number().min(1).max(100).optional().default(20),
      })
    )
    .query(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      const results = await db
        .select()
        .from(cheminformaticsResults)
        .where(eq(cheminformaticsResults.userId, ctx.user.id))
        .orderBy(desc(cheminformaticsResults.createdAt))
        .limit(input.limit);

      let filtered = results;

      if (input.workflow) {
        filtered = results.filter((r: any) => r.workflow === input.workflow);
      }

      if (input.query) {
        const queryLower = input.query.toLowerCase();
        filtered = filtered.filter(
          (r: any) =>
            r.inputSmiles.toLowerCase().includes(queryLower) ||
            r.canonicalSmiles?.toLowerCase().includes(queryLower)
        );
      }

      return filtered;
    }),

  /**
   * Delete a result
   */
  deleteResult: protectedProcedure
    .input(z.object({ id: z.number() }))
    .mutation(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      const result = await db
        .select()
        .from(cheminformaticsResults)
        .where(
          and(
            eq(cheminformaticsResults.id, input.id),
            eq(cheminformaticsResults.userId, ctx.user.id)
          )
        )
        .limit(1);

      if (!result.length) {
        throw new TRPCError({
          code: 'NOT_FOUND',
          message: 'Result not found',
        });
      }

      await db
        .delete(cheminformaticsResults)
        .where(eq(cheminformaticsResults.id, input.id));

      return {
        success: true,
        message: 'Result deleted successfully',
      };
    }),

  /**
   * Get statistics about user's pipeline runs
   */
  getStatistics: protectedProcedure.query(async ({ ctx }) => {
    if (!ctx.user) {
      throw new TRPCError({ code: 'UNAUTHORIZED' });
    }

    const db = await getDb();
    if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

    const allResults = await db
      .select()
      .from(cheminformaticsResults)
      .where(eq(cheminformaticsResults.userId, ctx.user.id));

    const byWorkflow = {
      similarity: allResults.filter((r: any) => r.workflow === 'similarity').length,
      brics: allResults.filter((r: any) => r.workflow === 'brics').length,
      validate: allResults.filter((r: any) => r.workflow === 'validate').length,
      full_pipeline: allResults.filter((r: any) => r.workflow === 'full_pipeline').length,
    };

    const byStatus = {
      completed: allResults.filter((r: any) => r.status === 'completed').length,
      failed: allResults.filter((r: any) => r.status === 'failed').length,
      pending: allResults.filter((r: any) => r.status === 'pending').length,
      running: allResults.filter((r: any) => r.status === 'running').length,
    };

    const avgExecutionTime =
      allResults.reduce((sum: number, r: any) => sum + (r.executionTime || 0), 0) /
      Math.max(allResults.length, 1);

    return {
      totalRuns: allResults.length,
      byWorkflow,
      byStatus,
      avgExecutionTime,
      lastRun: allResults.length > 0 ? allResults[0].createdAt : null,
    };
  }),

  /**
   * Export results as JSON
   */
  exportResults: protectedProcedure
    .input(
      z.object({
        ids: z.array(z.number()).optional(),
        workflow: z.enum(['similarity', 'brics', 'validate', 'full_pipeline']).optional(),
      })
    )
    .query(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      const allResults = await db
        .select()
        .from(cheminformaticsResults)
        .where(eq(cheminformaticsResults.userId, ctx.user.id));

      let results = allResults;

      if (input.ids && input.ids.length > 0) {
        results = allResults.filter((r: any) => input.ids!.includes(r.id));
      }

      if (input.workflow) {
        results = results.filter((r: any) => r.workflow === input.workflow);
      }

      return {
        success: true,
        data: results,
        count: results.length,
      };
    }),
});
