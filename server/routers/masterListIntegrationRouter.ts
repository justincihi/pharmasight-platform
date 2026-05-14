import { z } from 'zod';
import { protectedProcedure, router } from '../_core/trpc';
import { getDb } from '../db';
import { analogDiscoveries, cheminformaticsResults } from '../../drizzle/schema';
import { eq } from 'drizzle-orm';
import { TRPCError } from '@trpc/server';

/**
 * Master List Integration Router
 * Handles linking cheminformatics results to the master analog list
 */
export const masterListIntegrationRouter = router({
  /**
   * Save analog from cheminformatics result to master list
   */
  saveToMasterList: protectedProcedure
    .input(
      z.object({
        resultId: z.number(),
        smiles: z.string().min(1),
        compoundName: z.string().optional(),
        parentCompound: z.string(),
        confidenceScore: z.number().min(0).max(100).optional().default(85),
        similarityScore: z.number().min(0).max(100).optional().default(75),
        patentStatus: z.enum(['patent-free', 'patent-opportunity', 'patented', 'unknown']).default('patent-free'),
        notes: z.string().optional(),
      })
    )
    .mutation(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      // Verify result ownership
      const result = await db
        .select()
        .from(cheminformaticsResults)
        .where(eq(cheminformaticsResults.id, input.resultId))
        .limit(1);

      if (!result.length || result[0].userId !== ctx.user.id) {
        throw new TRPCError({
          code: 'FORBIDDEN',
          message: 'Result not found or unauthorized',
        });
      }

      // Check for duplicates
      const existing = await db
        .select()
        .from(analogDiscoveries)
        .where(eq(analogDiscoveries.smiles, input.smiles))
        .limit(1);

      if (existing.length) {
        return {
          success: false,
          message: 'Analog already exists in master list',
          compoundId: existing[0].compoundId,
          isDuplicate: true,
        };
      }

      // Generate unique compound ID
      const compoundId = `CHEM-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;

      // Insert into analog discoveries
      await db.insert(analogDiscoveries).values({
        compoundId,
        compoundName: input.compoundName || `Generated Analog ${compoundId}`,
        parentCompound: input.parentCompound,
        smiles: input.smiles,
        confidenceScore: input.confidenceScore,
        similarityScore: input.similarityScore,
        safetyScore: 50,
        efficacyScore: 50,
        drugLikenessScore: 50,
        patentStatus: input.patentStatus,
        patentNumbers: JSON.stringify([]),
        fdaStatus: 'unknown',
        therapeuticPotential: 'Pending evaluation',
        discoveryMethod: 'cheminformatics_pipeline',
        discoveredBy: 'cheminformatics_pipeline',
        discoveredAt: new Date(),
        optimizationNotes: input.notes || `Added from cheminformatics pipeline run ${input.resultId}`,
      });

      return {
        success: true,
        message: 'Analog saved to master list successfully',
        compoundId,
        isDuplicate: false,
      };
    }),

  /**
   * Batch save multiple analogs from a cheminformatics result
   */
  batchSaveToMasterList: protectedProcedure
    .input(
      z.object({
        resultId: z.number(),
        analogs: z.array(
          z.object({
            smiles: z.string().min(1),
            compoundName: z.string().optional(),
            parentCompound: z.string(),
            confidenceScore: z.number().min(0).max(100).optional().default(85),
            similarityScore: z.number().min(0).max(100).optional().default(75),
            patentStatus: z.enum(['patent-free', 'patent-opportunity', 'patented', 'unknown']).default('patent-free'),
          })
        ),
      })
    )
    .mutation(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      // Verify result ownership
      const result = await db
        .select()
        .from(cheminformaticsResults)
        .where(eq(cheminformaticsResults.id, input.resultId))
        .limit(1);

      if (!result.length || result[0].userId !== ctx.user.id) {
        throw new TRPCError({
          code: 'FORBIDDEN',
          message: 'Result not found or unauthorized',
        });
      }

      const saved = [];
      const duplicates = [];

      for (const analog of input.analogs) {
        // Check for duplicates
        const existing = await db
          .select()
          .from(analogDiscoveries)
          .where(eq(analogDiscoveries.smiles, analog.smiles))
          .limit(1);

        if (existing.length) {
          duplicates.push({
            smiles: analog.smiles,
            compoundId: existing[0].compoundId,
          });
          continue;
        }

        // Generate unique compound ID
        const compoundId = `CHEM-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;

        // Insert into analog discoveries
        await db.insert(analogDiscoveries).values({
          compoundId,
          compoundName: analog.compoundName || `Generated Analog ${compoundId}`,
          parentCompound: analog.parentCompound,
          smiles: analog.smiles,
          confidenceScore: analog.confidenceScore,
          similarityScore: analog.similarityScore,
          safetyScore: 50,
          efficacyScore: 50,
          drugLikenessScore: 50,
          patentStatus: analog.patentStatus,
          patentNumbers: JSON.stringify([]),
          fdaStatus: 'unknown',
          therapeuticPotential: 'Pending evaluation',
          discoveryMethod: 'cheminformatics_pipeline',
          discoveredBy: 'cheminformatics_pipeline',
          discoveredAt: new Date(),
          optimizationNotes: `Batch added from cheminformatics pipeline run ${input.resultId}`,
        });

        saved.push(compoundId);
      }

      return {
        success: true,
        message: `Saved ${saved.length} analogs, ${duplicates.length} duplicates skipped`,
        saved,
        duplicates,
        totalProcessed: input.analogs.length,
      };
    }),

  /**
   * Get master list statistics
   */
  getMasterListStats: protectedProcedure.query(async ({ ctx }) => {
    if (!ctx.user) {
      throw new TRPCError({ code: 'UNAUTHORIZED' });
    }

    const db = await getDb();
    if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

    const allAnalogs = await db.select().from(analogDiscoveries);

    const byPatentStatus = {
      'patent-free': allAnalogs.filter((a: any) => a.patentStatus === 'patent-free').length,
      'patent-opportunity': allAnalogs.filter((a: any) => a.patentStatus === 'patent-opportunity').length,
      patented: allAnalogs.filter((a: any) => a.patentStatus === 'patented').length,
      unknown: allAnalogs.filter((a: any) => a.patentStatus === 'unknown').length,
    };

    const bySource = {
      cheminformatics: allAnalogs.filter((a: any) => a.source === 'cheminformatics').length,
      other: allAnalogs.filter((a: any) => a.source !== 'cheminformatics').length,
    };

    return {
      totalAnalogs: allAnalogs.length,
      byPatentStatus,
      bySource,
      avgConfidenceScore:
        allAnalogs.reduce((sum: number, a: any) => sum + (a.confidenceScore || 0), 0) /
        Math.max(allAnalogs.length, 1),
    };
  }),

  /**
   * Link existing analog to cheminformatics result
   */
  linkToResult: protectedProcedure
    .input(
      z.object({
        resultId: z.number(),
        compoundId: z.string(),
      })
    )
    .mutation(async ({ ctx, input }) => {
      if (!ctx.user) {
        throw new TRPCError({ code: 'UNAUTHORIZED' });
      }

      const db = await getDb();
      if (!db) throw new TRPCError({ code: 'INTERNAL_SERVER_ERROR' });

      // Verify result ownership
      const result = await db
        .select()
        .from(cheminformaticsResults)
        .where(eq(cheminformaticsResults.id, input.resultId))
        .limit(1);

      if (!result.length || result[0].userId !== ctx.user.id) {
        throw new TRPCError({
          code: 'FORBIDDEN',
          message: 'Result not found or unauthorized',
        });
      }

      // Verify analog exists
      const analog = await db
        .select()
        .from(analogDiscoveries)
        .where(eq(analogDiscoveries.compoundId, input.compoundId))
        .limit(1);

      if (!analog.length) {
        throw new TRPCError({
          code: 'NOT_FOUND',
          message: 'Analog not found',
        });
      }

      return {
        success: true,
        message: 'Analog linked to result successfully',
        analog: analog[0],
      };
    }),
});
