import { protectedProcedure, router } from './_core/trpc';
import { uploadPDBFile, getUserPDBFiles, getPDBFileById, deletePDBFile, getPDBFilesByTarget } from './pdbStorage';
import { z } from 'zod';

export const pdbRouter = router({
  /**
   * Upload a new PDB receptor file
   */
  upload: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) {
        return { fileContent: '', fileName: '', targetName: '', description: '' };
      }
      const obj = val as Record<string, unknown>;
      return {
        fileContent: typeof obj.fileContent === 'string' ? obj.fileContent : '',
        fileName: typeof obj.fileName === 'string' ? obj.fileName : '',
        targetName: typeof obj.targetName === 'string' ? obj.targetName : '',
        description: typeof obj.description === 'string' ? obj.description : '',
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.fileContent || !input.fileName || !input.targetName) {
        throw new Error('File content, filename, and target name are required');
      }

      // Validate PDB format
      if (!input.fileContent.includes('ATOM') && !input.fileContent.includes('HETATM')) {
        throw new Error('Invalid PDB file: must contain ATOM or HETATM records');
      }

      return uploadPDBFile(
        input.fileContent,
        input.fileName,
        ctx.user.id.toString(),
        input.targetName,
        input.description || undefined
      );
    }),

  /**
   * Get all PDB files for the current user
   */
  list: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { targetName: '' };
      const obj = val as Record<string, unknown>;
      return {
        targetName: typeof obj.targetName === 'string' ? obj.targetName : '',
      };
    })
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (input.targetName) {
        return getPDBFilesByTarget(input.targetName);
      }

      return getUserPDBFiles(ctx.user.id.toString());
    }),

  /**
   * Get a specific PDB file by ID
   */
  getById: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { fileId: '' };
      const obj = val as Record<string, unknown>;
      return {
        fileId: typeof obj.fileId === 'string' ? obj.fileId : '',
      };
    })
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.fileId) {
        throw new Error('File ID is required');
      }

      return getPDBFileById(input.fileId, ctx.user.id.toString());
    }),

  /**
   * Delete a PDB file
   */
  delete: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { fileId: '' };
      const obj = val as Record<string, unknown>;
      return {
        fileId: typeof obj.fileId === 'string' ? obj.fileId : '',
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.fileId) {
        throw new Error('File ID is required');
      }

      return deletePDBFile(input.fileId, ctx.user.id.toString());
    }),
});
