/**
 * Batch Testing Router - tRPC procedures for analog batch testing
 */

import { router, protectedProcedure } from "../_core/trpc";
import { z } from "zod";
import {
  createBatchTestingJob,
  processBatchTestingJob,
  generateBatchSummary,
  type AnalogTestingJob,
} from "../services/analogBatchTestingService";

// In-memory job storage (in production, use database)
const batchJobs = new Map<string, AnalogTestingJob>();

export const batchTestingRouter = router({
  /**
   * Start a batch testing job for multiple analogs
   */
  startBatchTesting: protectedProcedure
    .input(
      z.object({
        analogIds: z.array(z.number()).min(1),
        receptorPath: z.string().optional(),
      })
    )
    .mutation(async ({ input, ctx }) => {
      try {
        // Create new batch job
        const job = await createBatchTestingJob(input.analogIds);

        // Store in memory
        batchJobs.set(job.id, job);

        // TODO: In production, save to database with user attribution
        console.log(
          `Batch testing job created: ${job.id} for user ${ctx.user?.id}`
        );

        return {
          success: true,
          jobId: job.id,
          message: `Batch testing started for ${input.analogIds.length} analogs`,
        };
      } catch (error) {
        console.error("Error starting batch testing:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * Get batch job status and results
   */
  getBatchJobStatus: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .query(({ input }) => {
      const job = batchJobs.get(input.jobId);

      if (!job) {
        return {
          success: false,
          error: "Job not found",
        };
      }

      // Calculate progress
      const completedCount = job.results.filter(
        (r) => r.status === "completed"
      ).length;
      const failedCount = job.results.filter((r) => r.status === "failed").length;
      const progress = Math.round(
        ((completedCount + failedCount) / job.results.length) * 100
      );

      return {
        success: true,
        job,
        progress,
        completedCount,
        failedCount,
        totalCount: job.results.length,
      };
    }),

  /**
   * Get batch job results and summary
   */
  getBatchJobResults: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .query(({ input }) => {
      const job = batchJobs.get(input.jobId);

      if (!job) {
        return {
          success: false,
          error: "Job not found",
        };
      }

      if (job.status !== "completed") {
        return {
          success: false,
          error: "Job is still running or failed",
          status: job.status,
        };
      }

      const summary = generateBatchSummary(job);

      return {
        success: true,
        job,
        summary,
        results: job.results,
      };
    }),

  /**
   * List all batch jobs for current user
   */
  listBatchJobs: protectedProcedure.query(({ ctx }) => {
    // TODO: In production, filter by user ID from database
    const jobs = Array.from(batchJobs.values());

    return {
      success: true,
      jobs: jobs.map((job) => ({
        id: job.id,
        status: job.status,
        createdAt: job.createdAt,
        completedAt: job.completedAt,
        totalAnalogs: job.results.length,
        completedTests: job.results.filter((r) => r.status === "completed").length,
        failedTests: job.results.filter((r) => r.status === "failed").length,
      })),
    };
  }),

  /**
   * Cancel a batch testing job
   */
  cancelBatchJob: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .mutation(({ input }) => {
      const job = batchJobs.get(input.jobId);

      if (!job) {
        return {
          success: false,
          error: "Job not found",
        };
      }

      if (job.status === "completed" || job.status === "failed") {
        return {
          success: false,
          error: "Cannot cancel a completed or failed job",
        };
      }

      job.status = "failed";
      job.error = "Job cancelled by user";
      job.completedAt = new Date();

      return {
        success: true,
        message: "Batch job cancelled",
      };
    }),

  /**
   * Export batch results as CSV
   */
  exportBatchResultsCSV: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .query(({ input }) => {
      const job = batchJobs.get(input.jobId);

      if (!job) {
        return {
          success: false,
          error: "Job not found",
        };
      }

      // Generate CSV
      const headers = [
        "Compound Name",
        "SMILES",
        "Status",
        "Binding Affinity",
        "RMSD",
        "hERG",
        "Hepatotoxicity",
        "Mutagenicity",
        "Carcinogenicity",
        "Tmax",
        "Half-life",
        "Error",
      ];

      const rows = job.results.map((result) => [
        result.compoundName,
        result.smiles,
        result.status,
        result.dockingResult?.bindingAffinity ?? "N/A",
        result.dockingResult?.rmsd ?? "N/A",
        result.admetResult?.herg ?? "N/A",
        result.admetResult?.hepatotoxicity ?? "N/A",
        result.admetResult?.mutagenicity ?? "N/A",
        result.admetResult?.carcinogenicity ?? "N/A",
        result.pkpdResult?.tmax ?? "N/A",
        result.pkpdResult?.halfLife ?? "N/A",
        result.error ?? "",
      ]);

      const csv = [headers, ...rows].map((row) => row.join(",")).join("\n");

      return {
        success: true,
        csv,
        filename: `batch-results-${input.jobId}.csv`,
      };
    }),
});
