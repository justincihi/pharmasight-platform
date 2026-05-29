import { protectedProcedure, router } from '../_core/trpc';
import { z } from 'zod';

/**
 * Lead Optimization Router
 * Handles metabolite prediction, ADMET property prediction, SAR analysis, and lead optimization
 */
export const leadOptimizationRouter = router({
  /**
   * Predict metabolites using biotransformer
   */
  predictMetabolites: protectedProcedure
    .input(
      z.object({
        smiles: z.string(),
        parentName: z.string().optional(),
        maxPhase1: z.number().optional().default(10),
        maxPhase2: z.number().optional().default(10),
        maxPhase3: z.number().optional().default(5),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');

        const result = (await executePythonScriptSafe('metabolite_predictor.py', 'predict_all_metabolites', [
          input.smiles,
          input.parentName || 'Unknown',
        ])) as any;

        if (result?.error) {
          throw new Error(result.error);
        }

        return {
          success: true,
          data: result?.data || result,
        };
      } catch (error) {
        console.error('[Lead Optimization] Metabolite prediction failed:', error);
        throw new Error(`Metabolite prediction failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Predict ADMET properties using ChemProp ADMET-AI
   */
  predictADMET: protectedProcedure
    .input(
      z.object({
        smiles: z.string(),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');

        const result = (await executePythonScriptSafe('chemprop_admet.py', 'predict_all_admet', [
          input.smiles,
        ])) as any;

        if (result?.error) {
          throw new Error(result.error);
        }

        return {
          success: true,
          data: result?.data || result,
        };
      } catch (error) {
        console.error('[Lead Optimization] ADMET prediction failed:', error);
        throw new Error(`ADMET prediction failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Analyze structure-activity relationships (SAR)
   */
  analyzeSAR: protectedProcedure
    .input(
      z.object({
        smiles: z.string(),
        compoundSeries: z
          .array(
            z.object({
              smiles: z.string(),
              activity: z.number().optional(),
            })
          )
          .optional(),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');

        const result = (await executePythonScriptSafe('sar_analyzer.py', 'analyze_sar', [
          input.smiles,
        ])) as any;

        if (result?.error) {
          throw new Error(result.error);
        }

        return {
          success: true,
          data: result?.data || result,
        };
      } catch (error) {
        console.error('[Lead Optimization] SAR analysis failed:', error);
        throw new Error(`SAR analysis failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Optimize leads using multi-objective optimization
   */
  optimizeLeads: protectedProcedure
    .input(
      z.object({
        parentSmiles: z.string(),
        targetSmiles: z.string().optional(),
        numAnalogs: z.number().optional().default(20),
        topN: z.number().optional().default(10),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');

        const result = (await executePythonScriptSafe('lead_optimizer.py', 'optimize_leads', [
          input.parentSmiles,
          input.targetSmiles || '',
          input.numAnalogs,
          input.topN,
        ])) as any;

        if (result?.error) {
          throw new Error(result.error);
        }

        return {
          success: true,
          data: result?.data || result,
        };
      } catch (error) {
        console.error('[Lead Optimization] Lead optimization failed:', error);
        throw new Error(`Lead optimization failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Run complete lead optimization pipeline
   */
  runCompletePipeline: protectedProcedure
    .input(
      z.object({
        parentSmiles: z.string(),
        targetSmiles: z.string().optional(),
        numAnalogs: z.number().optional().default(20),
        topN: z.number().optional().default(5),
        includeMetabolites: z.boolean().optional().default(true),
        includeSAR: z.boolean().optional().default(true),
      })
    )
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const results: Record<string, any> = {};
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');

        // Step 1: Predict ADMET for parent compound
        const admetResult = (await executePythonScriptSafe('chemprop_admet.py', 'predict_all_admet', [
          input.parentSmiles,
        ])) as any;
        if (!admetResult?.error) {
          results.parentADMET = admetResult?.data || admetResult;
        }

        // Step 2: Predict metabolites if requested
        if (input.includeMetabolites) {
          const metResult = (await executePythonScriptSafe('metabolite_predictor.py', 'predict_all_metabolites', [
            input.parentSmiles,
            'Parent',
          ])) as any;
          if (!metResult?.error) {
            results.metabolites = metResult?.data || metResult;
          }
        }

        // Step 3: Analyze SAR if requested
        if (input.includeSAR) {
          const sarResult = (await executePythonScriptSafe('sar_analyzer.py', 'analyze_sar', [
            input.parentSmiles,
          ])) as any;
          if (!sarResult?.error) {
            results.sar = sarResult?.data || sarResult;
          }
        }

        // Step 4: Optimize leads
        const optimResult = (await executePythonScriptSafe('lead_optimizer.py', 'optimize_leads', [
          input.parentSmiles,
          input.targetSmiles || '',
          input.numAnalogs,
          input.topN,
        ])) as any;
        if (!optimResult?.error) {
          results.optimizedLeads = optimResult?.data || optimResult;
        }

        return {
          success: true,
          data: results,
          pipelineSteps: {
            admet: !!results.parentADMET,
            metabolites: !!results.metabolites,
            sar: !!results.sar,
            optimization: !!results.optimizedLeads,
          },
        };
      } catch (error) {
        console.error('[Lead Optimization] Complete pipeline failed:', error);
        throw new Error(`Pipeline failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),

  /**
   * Compare multiple lead compounds
   */
  compareLeads: protectedProcedure
    .input(
      z.object({
        smilesList: z.array(z.string()).min(2).max(10),
        metrics: z.array(z.string()).optional(),
      })
    )
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      try {
        const { executePythonScriptSafe } = await import('../_core/pythonBridgeSafe');
        const comparisons = [];

        for (const smiles of input.smilesList) {
          const result = (await executePythonScriptSafe('chemprop_admet.py', 'predict_all_admet', [
            smiles,
          ])) as any;
          if (!result?.error) {
            comparisons.push({
              smiles,
              admet: result?.data || result,
            });
          }
        }

        return {
          success: true,
          comparisons,
          count: comparisons.length,
        };
      } catch (error) {
        console.error('[Lead Optimization] Lead comparison failed:', error);
        throw new Error(`Lead comparison failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
      }
    }),
});
