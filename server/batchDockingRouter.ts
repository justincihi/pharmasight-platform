import { router, protectedProcedure } from "./_core/trpc";
import { z } from "zod";
import { getDb } from "./db";
import {
  batchDockingJobs,
  batchDockingResults,
  analogDiscoveries,
} from "../drizzle/schema";
import { eq, and, desc } from "drizzle-orm";
import { TRPCError } from "@trpc/server";
import { runMolecularDocking } from "./molecularDockingWrapper";
import crypto from "crypto";

const generateJobId = () => `job_${crypto.randomBytes(8).toString("hex")}`;

let db: any = null;

// Initialize db on first use
async function initDb() {
  if (!db) {
    db = await getDb();
    if (!db) {
      throw new TRPCError({
        code: "INTERNAL_SERVER_ERROR",
        message: "Database connection failed",
      });
    }
  }
  return db;
}

export const batchDockingRouter = router({
  /**
   * Submit a new batch docking job
   */  submitJob: protectedProcedure
    .input(
      z.object({
        jobName: z.string().min(1, "Job name required"),
        compoundIds: z.array(z.number()).min(1, "At least one compound required"),
        targetName: z.string().min(1, "Target name required"),
        parametersId: z.number().optional(),
      })
    )
    .mutation(async ({ ctx, input }: any) => {
      const jobId = generateJobId();
      const database = await initDb();

      // Create batch job record
      const result = await database
        .insert(batchDockingJobs)
        .values({
          id: jobId,
          userId: ctx.user.id.toString(),
          jobName: input.jobName,
          status: "pending",
          totalCompounds: input.compoundIds.length,
          completedCompounds: 0,
          failedCompounds: 0,
          targetName: input.targetName,
          parametersId: input.parametersId,
        })
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to create batch job: ${err.message}`,
          });
        });

      // Create result records for each compound
      for (const analogId of input.compoundIds) {
        await db
          .insert(batchDockingResults)
          .values({
            jobId,
            analogId,
            status: "pending",
          })
          .catch((err: any) => {
            console.error(`Failed to create result record: ${err.message}`);
          });
      }

      // Start batch processing in background (non-blocking)
      processBatchJob(jobId, input.compoundIds, input.targetName, ctx.user.id.toString()).catch(
        (err) => console.error(`Batch processing error: ${err.message}`)
      );

      return {
        success: true,
        jobId,
        message: `Batch job submitted with ${input.compoundIds.length} compounds`,
      };
    }),

  /**
   * List all batch jobs for the current user
   */
  listJobs: protectedProcedure    .input(
      z.object({
        limit: z.number().default(20),
        offset: z.number().default(0),
      })
    )
    .query(async ({ ctx, input }: any) => {
      const database = await initDb();
      const jobs = await database
        .select()
        .from(batchDockingJobs)
        .where(eq(batchDockingJobs.userId, ctx.user.id.toString()))
        .orderBy(desc(batchDockingJobs.createdAt))
        .limit(input.limit)
        .offset(input.offset)
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch job: ${err.message}`,
          });
        });

      return jobs;
    }),

  /**
   * Get detailed information about a specific batch job
   */  getJobDetails: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .query(async ({ ctx, input }: any) => {
      const database = await initDb();
      const job = await database
        .select()
        .from(batchDockingJobs)
        .where(
          and(
            eq(batchDockingJobs.id, input.jobId),
            eq(batchDockingJobs.userId, ctx.user.id.toString())
          )
        )
        .then((rows: any) => rows[0])
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch job: ${err.message}`,
          });
        });

      if (!job) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: "Job not found",
        });
      }

      // Get results for this job
      const results = await database
        .select()
        .from(batchDockingResults)
        .where(eq(batchDockingResults.jobId, input.jobId))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch results: ${err.message}`,
          });
        });

      return {
        job,
        results,
        progressPercentage: Math.round(
          ((job.completedCompounds + job.failedCompounds) / job.totalCompounds) * 100
        ),
      };
    }),

  /**
   * Cancel a batch job
   */  cancelJob: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .mutation(async ({ ctx, input }: any) => {
      const database = await initDb();
      const job = await database
        .select()
        .from(batchDockingJobs)
        .where(
          and(
            eq(batchDockingJobs.id, input.jobId),
            eq(batchDockingJobs.userId, ctx.user.id.toString())
          )
        )
        .then((rows: any) => rows[0]);

      if (!job) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: "Job not found",
        });
      }

      if ((job as any).status === "completed" || (job as any).status === "failed") {
        throw new TRPCError({
          code: "BAD_REQUEST",
          message: "Cannot cancel a completed or failed job",
        });
      }

      await database
        .update(batchDockingJobs)
        .set({ status: "cancelled", updatedAt: new Date() })
        .where(eq(batchDockingJobs.id, input.jobId))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to cancel job: ${err.message}`,
          });
        });

      return { success: true, message: "Job cancelled successfully" };
    }),

  /**
   * Export batch results as CSV or JSON
   */  exportResults: protectedProcedure
    .input(
      z.object({
        jobId: z.string(),
        format: z.enum(["csv", "json"]),
      })
    )
    .query(async ({ ctx, input }: any) => {
      const database = await initDb();
      const job = await database
        .select()
        .from(batchDockingJobs)
        .where(
          and(
            eq(batchDockingJobs.id, input.jobId),
            eq(batchDockingJobs.userId, ctx.user.id.toString())
          )
        )
        .then((rows: any) => rows[0]);

      if (!job) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: "Job not found",
        });
      }

      // Get results with analog details
      const results = await database
        .select({
          result: batchDockingResults,
          analog: analogDiscoveries,
        })
        .from(batchDockingResults)
        .innerJoin(
          analogDiscoveries,
          eq(batchDockingResults.analogId, analogDiscoveries.id)
        )
        .where(eq(batchDockingResults.jobId, input.jobId))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch results: ${err.message}`,
          });
        });

      if (input.format === "csv") {
        return generateCSVExport(job, results);
      } else {
        return generateJSONExport(job, results);
      }
    }),

  /**
   * Get batch job statistics
   */  getStatistics: protectedProcedure
    .input(z.object({ jobId: z.string() }))
    .query(async ({ ctx, input }: any) => {
      const database = await initDb();
      const results = await database
        .select()
        .from(batchDockingResults)
        .where(eq(batchDockingResults.jobId, input.jobId))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch results: ${err.message}`,
          });
        });

      const completedResults = results.filter((r: any) => r.status === "completed");
      const affinities = completedResults
        .map((r: any) => {
          if (!r.bindingAffinity) return null;
          const val = parseFloat(r.bindingAffinity);
          return isNaN(val) ? null : val;
        })
        .filter((v: any) => v !== null) as number[];

      const stats = {
        totalCompounds: results.length,
        completedCompounds: completedResults.length,
        failedCompounds: results.filter((r: any) => r.status === "failed").length,
        pendingCompounds: results.filter((r: any) => r.status === "pending").length,
        successRate: results.length > 0 ? (completedResults.length / results.length) * 100 : 0,
        averageAffinity: affinities.length > 0 ? affinities.reduce((a, b) => a + b) / affinities.length : 0,
        bestAffinity: affinities.length > 0 ? Math.min(...affinities) : null,
        worstAffinity: affinities.length > 0 ? Math.max(...affinities) : null,
      };

      return stats;
    }),
});

/**
 * Background job processor - runs docking for all compounds in batch
 */
async function processBatchJob(
  jobId: string,
  compoundIds: number[],
  targetName: string,
  userId: string
) {
  try {
    const database = await initDb();

    // Update job status to running
    await database
      .update(batchDockingJobs)
      .set({
        status: "running",
        startedAt: new Date(),
        updatedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, jobId));

    let completedCount = 0;
    let failedCount = 0;

    // Process each compound
    for (const analogId of compoundIds) {
      try {
        // Get analog details
        const analog = await db
          .select()
          .from(analogDiscoveries)
        .where(eq(analogDiscoveries.id, analogId))
        .then((rows: any) => rows[0]);

        if (!analog) {
          await db
            .update(batchDockingResults)
            .set({
              status: "failed",
              errorMessage: "Analog not found",
              updatedAt: new Date(),
            })
            .where(
              and(
                eq(batchDockingResults.jobId, jobId),
                eq(batchDockingResults.analogId, analogId)
              )
            );
          failedCount++;
          continue;
        }

        // Run docking
        const dockingResult = await runMolecularDocking(
          analog.smiles,
          analog.id.toString(),
          targetName
        );

        // Update result record
        await database
          .update(batchDockingResults)
          .set({
            status: "completed",
            bindingAffinity: dockingResult.binding_affinity?.toString(),
            dockingScore: Math.round((dockingResult.binding_affinity || 0) * 10),
            numPoses: dockingResult.num_poses,
            topPoses: JSON.stringify(dockingResult.top_poses || []),
            completedAt: new Date(),
            updatedAt: new Date(),
          })
          .where(
            and(
              eq(batchDockingResults.jobId, jobId),
              eq(batchDockingResults.analogId, analogId)
            )
          );

        completedCount++;
      } catch (error: any) {
        console.error(`Docking failed for analog ${analogId}: ${error.message}`);
        await database
          .update(batchDockingResults)
          .set({
            status: "failed",
            errorMessage: error.message,
            updatedAt: new Date(),
          })
          .where(
            and(
              eq(batchDockingResults.jobId, jobId),
              eq(batchDockingResults.analogId, analogId)
            )
          );
        failedCount++;
      }

      // Update job progress
      await database
        .update(batchDockingJobs)
        .set({
          completedCompounds: completedCount,
          failedCompounds: failedCount,
          updatedAt: new Date(),
        })
        .where(eq(batchDockingJobs.id, jobId));
    }

    // Mark job as completed
    await db
      .update(batchDockingJobs)
      .set({
        status: "completed",
        completedAt: new Date(),
        updatedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, jobId));
  } catch (error: any) {
    console.error(`Batch job processing failed: ${error.message}`);
    await db
      .update(batchDockingJobs)
      .set({
        status: "failed",
        updatedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, jobId));
  }
}

/**
 * Generate CSV export
 */
function generateCSVExport(job: any, results: any[]): string {
  const headers = [
    "Compound Name",
    "SMILES",
    "Status",
    "Binding Affinity (kcal/mol)",
    "Docking Score",
    "Number of Poses",
    "Error Message",
  ];

  const rows = results.map((item) => [
    item.analog.compoundName,
    item.analog.smiles,
    item.result.status,
    item.result.bindingAffinity || "N/A",
    item.result.dockingScore || "N/A",
    item.result.numPoses || "N/A",
    item.result.errorMessage || "",
  ]);

  const csvContent = [
    `Batch Docking Results - ${job.jobName}`,
    `Target: ${job.targetName}`,
    `Date: ${new Date().toISOString()}`,
    `Status: ${job.status}`,
    `Completed: ${job.completedCompounds}/${job.totalCompounds}`,
    "",
    headers.join(","),
    ...rows.map((row: any) =>
      row.map((cell: any) => `"${cell}"`).join(",")
    ),
  ].join("\n");

  return csvContent;
}

/**
 * Generate JSON export
 */
function generateJSONExport(job: any, results: any[]): string {
  const exportData = {
    job: {
      id: job.id,
      name: job.jobName,
      target: job.targetName,
      status: job.status,
      createdAt: job.createdAt,
      completedAt: job.completedAt,
      totalCompounds: job.totalCompounds,
      completedCompounds: job.completedCompounds,
      failedCompounds: job.failedCompounds,
    },
    results: results.map((item) => ({
      compound: {
        id: item.analog.id,
        name: item.analog.compoundName,
        smiles: item.analog.smiles,
        parentCompound: item.analog.parentCompound,
      },
      docking: {
        status: item.result.status,
        bindingAffinity: item.result.bindingAffinity,
        dockingScore: item.result.dockingScore,
        numPoses: item.result.numPoses,
        topPoses: item.result.topPoses ? JSON.parse(item.result.topPoses) : [],
        errorMessage: item.result.errorMessage,
        completedAt: item.result.completedAt,
      },
    })),
  };

  return JSON.stringify(exportData, null, 2);
}
