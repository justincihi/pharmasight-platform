import { router, publicProcedure, protectedProcedure } from '../_core/trpc';
import { z } from 'zod';
import {
  confirmAndFetchSimilars,
  checkPatentStatus,
  screenAndFlagForMasterlist,
  generateBricsAnalogs,
  enumerateSubstituentAnalogs,
  fullAnalogPipeline,
  validateSmiles,
} from '../_core/cheminformatics';

/**
 * Cheminformatics Router
 * Exposes PubChem similarity screening, patent checking, and analog generation workflows
 */

export const cheminformaticsRouter = router({
  /**
   * Validate a SMILES string
   */
  validateSmiles: publicProcedure
    .input(
      z.object({
        smiles: z.string().min(1, 'SMILES cannot be empty'),
      })
    )
    .query(async ({ input }) => {
      const result = await validateSmiles(input.smiles);
      return result;
    }),

  /**
   * Confirm compound in PubChem and fetch similar compounds
   */
  confirmAndFetchSimilars: publicProcedure
    .input(
      z.object({
        nameOrSmiles: z.string().min(1, 'Compound name or SMILES required'),
        threshold: z.number().min(0).max(1).default(0.70),
        maxHits: z.number().min(1).max(100).default(25),
      })
    )
    .query(async ({ input }) => {
      const result = await confirmAndFetchSimilars(
        input.nameOrSmiles,
        input.threshold,
        input.maxHits
      );
      return result;
    }),

  /**
   * Check patent status for a compound by CID
   */
  checkPatentStatus: publicProcedure
    .input(
      z.object({
        cid: z.number().int().positive('CID must be a positive integer'),
      })
    )
    .query(async ({ input }) => {
      const result = await checkPatentStatus(input.cid);
      return result;
    }),

  /**
   * Screen compounds and flag patent-free candidates
   */
  screenAndFlagForMasterlist: publicProcedure
    .input(
      z.object({
        hits: z.array(
          z.object({
            cid: z.number().int(),
            name: z.string(),
            smiles: z.string(),
            mw: z.number(),
            tanimoto: z.number(),
          })
        ),
        maxHits: z.number().min(1).max(100).default(25),
      })
    )
    .query(async ({ input }) => {
      const result = await screenAndFlagForMasterlist(input.hits, input.maxHits);
      return result;
    }),

  /**
   * Generate BRICS analogs from a parent SMILES
   */
  generateBricsAnalogs: publicProcedure
    .input(
      z.object({
        smiles: z.string().min(1, 'SMILES required'),
        n: z.number().min(1).max(100).default(25),
      })
    )
    .query(async ({ input }) => {
      const result = await generateBricsAnalogs(input.smiles, input.n);
      return result;
    }),

  /**
   * Enumerate substituent analogs with R-group substitution
   */
  enumerateSubstituentAnalogs: publicProcedure
    .input(
      z.object({
        baseSmiles: z.string().min(1, 'Base SMILES required'),
        attachmentIdx: z.number().int().min(0, 'Attachment index must be >= 0'),
        n: z.number().min(1).max(100).default(25),
      })
    )
    .query(async ({ input }) => {
      const result = await enumerateSubstituentAnalogs(
        input.baseSmiles,
        input.attachmentIdx,
        input.n
      );
      return result;
    }),

  /**
   * Full integrated analog pipeline
   * Orchestrates all workflows: canonicalize → generate → screen → patent check
   */
  fullAnalogPipeline: publicProcedure
    .input(
      z.object({
        inputSmiles: z.string().min(1, 'Input SMILES required'),
        threshold: z.number().min(0).max(1).default(0.70),
        maxHits: z.number().min(1).max(100).default(25),
      })
    )
    .query(async ({ input }) => {
      const result = await fullAnalogPipeline(
        input.inputSmiles,
        input.threshold,
        input.maxHits
      );
      return result;
    }),
});
