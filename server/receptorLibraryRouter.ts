import { router, protectedProcedure, adminProcedure } from "./_core/trpc";
import { z } from "zod";
import { getDb } from "./db";
import { receptorLibrary } from "../drizzle/schema";
import { eq, like, desc } from "drizzle-orm";
import { TRPCError } from "@trpc/server";
import { storagePut } from "./storage";
import fs from "fs";
import path from "path";

let db: any = null;

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

export const receptorLibraryRouter = router({
  /**
   * List all available receptors in the library
   */
  listReceptors: protectedProcedure
    .input(
      z.object({
        category: z.string().optional(),
        search: z.string().optional(),
        limit: z.number().default(50),
        offset: z.number().default(0),
      })
    )
    .query(async ({ input }: any) => {
      const database = await initDb();

      let query = database
        .select()
        .from(receptorLibrary)
        .where(eq(receptorLibrary.isActive, true));

      if (input.category) {
        query = query.where(eq(receptorLibrary.category, input.category));
      }

      if (input.search) {
        query = query.where(
          like(receptorLibrary.targetName, `%${input.search}%`)
        );
      }

      const receptors = await query
        .orderBy(desc(receptorLibrary.createdAt))
        .limit(input.limit)
        .offset(input.offset)
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch receptors: ${err.message}`,
          });
        });

      return receptors;
    }),

  /**
   * Get a specific receptor by target name
   */
  getReceptor: protectedProcedure
    .input(z.object({ targetName: z.string() }))
    .query(async ({ input }: any) => {
      const database = await initDb();

      const receptor = await database
        .select()
        .from(receptorLibrary)
        .where(eq(receptorLibrary.targetName, input.targetName))
        .then((rows: any) => rows[0])
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to fetch receptor: ${err.message}`,
          });
        });

      if (!receptor) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: `Receptor ${input.targetName} not found`,
        });
      }

      return receptor;
    }),

  /**
   * Get download URL for a receptor PDBQT file
   */
  getDownloadUrl: protectedProcedure
    .input(z.object({ targetName: z.string(), expiresIn: z.number().default(3600) }))
    .query(async ({ input }: any) => {
      const database = await initDb();

      const receptor = await database
        .select()
        .from(receptorLibrary)
        .where(eq(receptorLibrary.targetName, input.targetName))
        .then((rows: any) => rows[0]);

      if (!receptor) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: `Receptor ${input.targetName} not found`,
        });
      }

      // Return the S3 URL (already public)
      return {
        targetName: receptor.targetName,
        pdbqtUrl: receptor.pdbqtUrl,
        pdbUrl: receptor.pdbUrl,
        fileName: receptor.pdbqtFileName,
      };
    }),

  /**
   * Admin: Upload a new receptor to the library
   */
  uploadReceptor: adminProcedure
    .input(
      z.object({
        targetName: z.string().min(1),
        description: z.string().min(1),
        pdbId: z.string().min(1),
        category: z.string().optional(),
        tags: z.array(z.string()).optional(),
        notes: z.string().optional(),
        defaultBoxCenterX: z.number().optional(),
        defaultBoxCenterY: z.number().optional(),
        defaultBoxCenterZ: z.number().optional(),
        defaultBoxSizeX: z.number().optional(),
        defaultBoxSizeY: z.number().optional(),
        defaultBoxSizeZ: z.number().optional(),
        defaultExhaustiveness: z.number().optional(),
        defaultNumPoses: z.number().optional(),
      })
    )
    .mutation(async ({ ctx, input }: any) => {
      const database = await initDb();

      // Check if receptor already exists
      const existing = await database
        .select()
        .from(receptorLibrary)
        .where(eq(receptorLibrary.targetName, input.targetName))
        .then((rows: any) => rows[0]);

      if (existing) {
        throw new TRPCError({
          code: "BAD_REQUEST",
          message: `Receptor ${input.targetName} already exists`,
        });
      }

      // Upload PDBQT file to S3
      const receptorsDir = path.join(process.cwd(), "receptors");
      const pdbqtFileName = `${input.targetName}_${input.pdbId}.pdbqt`;
      const pdbqtPath = path.join(receptorsDir, pdbqtFileName);

      if (!fs.existsSync(pdbqtPath)) {
        throw new TRPCError({
          code: "BAD_REQUEST",
          message: `PDBQT file not found: ${pdbqtFileName}`,
        });
      }

      const fileContent = fs.readFileSync(pdbqtPath);
      const fileSize = fileContent.length;

      // Upload to S3
      const s3Key = `receptors/${input.targetName}/${pdbqtFileName}`;
      const { url: pdbqtUrl } = await storagePut(s3Key, fileContent, "text/plain");

      // Create database record
      const result = await database
        .insert(receptorLibrary)
        .values({
          targetName: input.targetName,
          description: input.description,
          pdbId: input.pdbId,
          pdbFileName: `${input.targetName}_${input.pdbId}.pdb`,
          pdbqtFileName: pdbqtFileName,
          pdbqtUrl: pdbqtUrl,
          fileSize: fileSize,
          category: input.category || "psychiatric",
          tags: input.tags ? JSON.stringify(input.tags) : null,
          notes: input.notes,
          uploadedBy: ctx.user.id.toString(),
          defaultBoxCenterX: input.defaultBoxCenterX,
          defaultBoxCenterY: input.defaultBoxCenterY,
          defaultBoxCenterZ: input.defaultBoxCenterZ,
          defaultBoxSizeX: input.defaultBoxSizeX,
          defaultBoxSizeY: input.defaultBoxSizeY,
          defaultBoxSizeZ: input.defaultBoxSizeZ,
          defaultExhaustiveness: input.defaultExhaustiveness || 8,
          defaultNumPoses: input.defaultNumPoses || 5,
        })
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to create receptor: ${err.message}`,
          });
        });

      return {
        success: true,
        targetName: input.targetName,
        pdbqtUrl: pdbqtUrl,
        fileSize: fileSize,
      };
    }),

  /**
   * Admin: Delete a receptor from the library
   */
  deleteReceptor: adminProcedure
    .input(z.object({ targetName: z.string() }))
    .mutation(async ({ input }: any) => {
      const database = await initDb();

      const receptor = await database
        .select()
        .from(receptorLibrary)
        .where(eq(receptorLibrary.targetName, input.targetName))
        .then((rows: any) => rows[0]);

      if (!receptor) {
        throw new TRPCError({
          code: "NOT_FOUND",
          message: `Receptor ${input.targetName} not found`,
        });
      }

      // Mark as inactive instead of deleting
      await database
        .update(receptorLibrary)
        .set({
          isActive: false,
          updatedAt: new Date(),
        })
        .where(eq(receptorLibrary.targetName, input.targetName))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to delete receptor: ${err.message}`,
          });
        });

      return {
        success: true,
        targetName: input.targetName,
      };
    }),

  /**
   * Admin: Update receptor metadata
   */
  updateReceptor: adminProcedure
    .input(
      z.object({
        targetName: z.string(),
        description: z.string().optional(),
        notes: z.string().optional(),
        tags: z.array(z.string()).optional(),
        defaultBoxCenterX: z.number().optional(),
        defaultBoxCenterY: z.number().optional(),
        defaultBoxCenterZ: z.number().optional(),
        defaultBoxSizeX: z.number().optional(),
        defaultBoxSizeY: z.number().optional(),
        defaultBoxSizeZ: z.number().optional(),
        defaultExhaustiveness: z.number().optional(),
        defaultNumPoses: z.number().optional(),
      })
    )
    .mutation(async ({ input }: any) => {
      const database = await initDb();

      const updates: any = {
        updatedAt: new Date(),
      };

      if (input.description) updates.description = input.description;
      if (input.notes) updates.notes = input.notes;
      if (input.tags) updates.tags = JSON.stringify(input.tags);
      if (input.defaultBoxCenterX !== undefined) updates.defaultBoxCenterX = input.defaultBoxCenterX;
      if (input.defaultBoxCenterY !== undefined) updates.defaultBoxCenterY = input.defaultBoxCenterY;
      if (input.defaultBoxCenterZ !== undefined) updates.defaultBoxCenterZ = input.defaultBoxCenterZ;
      if (input.defaultBoxSizeX !== undefined) updates.defaultBoxSizeX = input.defaultBoxSizeX;
      if (input.defaultBoxSizeY !== undefined) updates.defaultBoxSizeY = input.defaultBoxSizeY;
      if (input.defaultBoxSizeZ !== undefined) updates.defaultBoxSizeZ = input.defaultBoxSizeZ;
      if (input.defaultExhaustiveness !== undefined) updates.defaultExhaustiveness = input.defaultExhaustiveness;
      if (input.defaultNumPoses !== undefined) updates.defaultNumPoses = input.defaultNumPoses;

      await database
        .update(receptorLibrary)
        .set(updates)
        .where(eq(receptorLibrary.targetName, input.targetName))
        .catch((err: any) => {
          throw new TRPCError({
            code: "INTERNAL_SERVER_ERROR",
            message: `Failed to update receptor: ${err.message}`,
          });
        });

      return {
        success: true,
        targetName: input.targetName,
      };
    }),

  /**
   * Get receptor statistics
   */
  getStatistics: protectedProcedure.query(async (): Promise<any> => {
    const database = await initDb();

    const receptors = await database
      .select()
      .from(receptorLibrary)
      .where(eq(receptorLibrary.isActive, true))
      .catch((err: any) => {
        throw new TRPCError({
          code: "INTERNAL_SERVER_ERROR",
          message: `Failed to fetch statistics: ${err.message}`,
        });
      });

    const categories = new Set(receptors.map((r: any) => r.category));
    const totalSize = receptors.reduce((sum: number, r: any) => sum + (r.fileSize || 0), 0);

    return {
      totalReceptors: receptors.length,
      categories: Array.from(categories),
      totalFileSize: totalSize,
      averageFileSize: receptors.length > 0 ? totalSize / receptors.length : 0,
    };
  }),
});
