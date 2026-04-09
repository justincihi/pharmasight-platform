import { router, protectedProcedure } from "../_core/trpc";
import { z } from "zod";
import {
  runBatchKetamineDocking,
  calculateSelectivityScores,
  exportBatchResultsToCSV,
  exportBatchResultsToJSON,
  KetamineBatchConfig,
} from "../services/batchKetamineService";

export const batchKetamineRouter = router({
  /**
   * Run batch docking for all ketamine analogs
   */
  runBatchDocking: protectedProcedure
    .input(
      z.object({
        compoundFilter: z.string().optional(),
        receptors: z.array(z.string()),
        dockingParams: z
          .object({
            exhaustiveness: z.number().optional(),
            numPoses: z.number().optional(),
            centerX: z.number().optional(),
            centerY: z.number().optional(),
            centerZ: z.number().optional(),
            sizeX: z.number().optional(),
            sizeY: z.number().optional(),
            sizeZ: z.number().optional(),
          })
          .optional(),
      })
    )
    .mutation(async ({ input }) => {
      try {
        const config: KetamineBatchConfig = {
          compoundFilter: input.compoundFilter,
          receptors: input.receptors,
          dockingParams: input.dockingParams,
        };

        const result = await runBatchKetamineDocking(config);
        return result;
      } catch (error) {
        const errorMsg = error instanceof Error ? error.message : String(error);
        console.error("[Batch Ketamine Router] Error:", errorMsg);
        throw new Error(`Batch docking failed: ${errorMsg}`);
      }
    }),

  /**
   * Calculate selectivity scores for batch results
   */
  calculateSelectivity: protectedProcedure
    .input(
      z.object({
        results: z.array(
          z.object({
            analogId: z.string(),
            compoundName: z.string(),
            smiles: z.string(),
            receptor: z.string(),
            bindingAffinity: z.number(),
            rmsd: z.number(),
            numPoses: z.number(),
            timestamp: z.date(),
          })
        ),
        primaryReceptor: z.string(),
      })
    )
    .query(({ input }) => {
      try {
        const selectivityScores = calculateSelectivityScores(
          input.results,
          input.primaryReceptor
        );
        return selectivityScores;
      } catch (error) {
        const errorMsg = error instanceof Error ? error.message : String(error);
        console.error("[Batch Ketamine Router] Selectivity calculation error:", errorMsg);
        throw new Error(`Selectivity calculation failed: ${errorMsg}`);
      }
    }),

  /**
   * Export batch results to CSV
   */
  exportToCSV: protectedProcedure
    .input(
      z.object({
        results: z.array(
          z.object({
            analogId: z.string(),
            compoundName: z.string(),
            smiles: z.string(),
            receptor: z.string(),
            bindingAffinity: z.number(),
            rmsd: z.number(),
            numPoses: z.number(),
            timestamp: z.date(),
          })
        ),
      })
    )
    .query(({ input }) => {
      try {
        const csv = exportBatchResultsToCSV(input.results);
        return { success: true, csv };
      } catch (error) {
        const errorMsg = error instanceof Error ? error.message : String(error);
        console.error("[Batch Ketamine Router] CSV export error:", errorMsg);
        throw new Error(`CSV export failed: ${errorMsg}`);
      }
    }),

  /**
   * Export batch results to JSON
   */
  exportToJSON: protectedProcedure
    .input(
      z.object({
        results: z.array(
          z.object({
            analogId: z.string(),
            compoundName: z.string(),
            smiles: z.string(),
            receptor: z.string(),
            bindingAffinity: z.number(),
            rmsd: z.number(),
            numPoses: z.number(),
            timestamp: z.date(),
          })
        ),
      })
    )
    .query(({ input }) => {
      try {
        const json = exportBatchResultsToJSON(input.results);
        return { success: true, json };
      } catch (error) {
        const errorMsg = error instanceof Error ? error.message : String(error);
        console.error("[Batch Ketamine Router] JSON export error:", errorMsg);
        throw new Error(`JSON export failed: ${errorMsg}`);
      }
    }),
});
