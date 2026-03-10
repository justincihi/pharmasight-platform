import { protectedProcedure, router } from './_core/trpc';
import { getDb } from './db';
import { dockingParameters } from '../drizzle/schema';
import { eq, and } from 'drizzle-orm';

export const dockingParametersRouter = router({
  /**
   * Create a new docking parameter configuration
   */
  create: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) {
        return {
          name: '',
          targetName: '',
          boxCenterX: '0',
          boxCenterY: '0',
          boxCenterZ: '0',
          boxSizeX: '20',
          boxSizeY: '20',
          boxSizeZ: '20',
          exhaustiveness: 8,
          numPoses: 9,
          isDefault: false,
        };
      }
      const obj = val as Record<string, unknown>;
      return {
        name: typeof obj.name === 'string' ? obj.name : '',
        targetName: typeof obj.targetName === 'string' ? obj.targetName : '',
        boxCenterX: typeof obj.boxCenterX === 'string' ? obj.boxCenterX : '0',
        boxCenterY: typeof obj.boxCenterY === 'string' ? obj.boxCenterY : '0',
        boxCenterZ: typeof obj.boxCenterZ === 'string' ? obj.boxCenterZ : '0',
        boxSizeX: typeof obj.boxSizeX === 'string' ? obj.boxSizeX : '20',
        boxSizeY: typeof obj.boxSizeY === 'string' ? obj.boxSizeY : '20',
        boxSizeZ: typeof obj.boxSizeZ === 'string' ? obj.boxSizeZ : '20',
        exhaustiveness: typeof obj.exhaustiveness === 'number' ? obj.exhaustiveness : 8,
        numPoses: typeof obj.numPoses === 'number' ? obj.numPoses : 9,
        isDefault: typeof obj.isDefault === 'boolean' ? obj.isDefault : false,
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.name || !input.targetName) {
        throw new Error('Name and target name are required');
      }

      // Validate exhaustiveness range
      if (input.exhaustiveness < 1 || input.exhaustiveness > 32) {
        throw new Error('Exhaustiveness must be between 1 and 32');
      }

      // Validate num poses
      if (input.numPoses < 1 || input.numPoses > 20) {
        throw new Error('Number of poses must be between 1 and 20');
      }

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      // If setting as default, unset other defaults for this target
      if (input.isDefault) {
        await db
          .update(dockingParameters)
          .set({ isDefault: 0 })
          .where(eq(dockingParameters.targetName, input.targetName));
      }

      await db.insert(dockingParameters).values({
        userId: ctx.user.id.toString(),
        name: input.name,
        targetName: input.targetName,
        boxCenterX: input.boxCenterX,
        boxCenterY: input.boxCenterY,
        boxCenterZ: input.boxCenterZ,
        boxSizeX: input.boxSizeX,
        boxSizeY: input.boxSizeY,
        boxSizeZ: input.boxSizeZ,
        exhaustiveness: input.exhaustiveness,
        numPoses: input.numPoses,
        isDefault: input.isDefault ? 1 : 0,
      });

      return { success: true };
    }),

  /**
   * Get all docking parameters for the current user
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

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      let query = db
        .select()
        .from(dockingParameters)
        .where(eq(dockingParameters.userId, ctx.user.id.toString()));

      if (input.targetName) {
        query = db
          .select()
          .from(dockingParameters)
          .where(
            and(
              eq(dockingParameters.userId, ctx.user.id.toString()),
              eq(dockingParameters.targetName, input.targetName)
            )
          );
      }

      return query;
    }),

  /**
   * Get default parameters for a target
   */
  getDefault: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { targetName: '' };
      const obj = val as Record<string, unknown>;
      return {
        targetName: typeof obj.targetName === 'string' ? obj.targetName : '',
      };
    })
    .query(async ({ input, ctx }) => {
      if (!input.targetName) {
        throw new Error('Target name is required');
      }

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      const result = await db
        .select()
        .from(dockingParameters)
        .where(
          and(
            eq(dockingParameters.targetName, input.targetName),
            eq(dockingParameters.isDefault, 1)
          )
        )
        .limit(1);

      return result[0] || null;
    }),

  /**
   * Update docking parameters
   */
  update: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { id: 0 };
      const obj = val as Record<string, unknown>;
      return {
        id: typeof obj.id === 'number' ? obj.id : 0,
        name: typeof obj.name === 'string' ? obj.name : '',
        boxCenterX: typeof obj.boxCenterX === 'string' ? obj.boxCenterX : '0',
        boxCenterY: typeof obj.boxCenterY === 'string' ? obj.boxCenterY : '0',
        boxCenterZ: typeof obj.boxCenterZ === 'string' ? obj.boxCenterZ : '0',
        boxSizeX: typeof obj.boxSizeX === 'string' ? obj.boxSizeX : '20',
        boxSizeY: typeof obj.boxSizeY === 'string' ? obj.boxSizeY : '20',
        boxSizeZ: typeof obj.boxSizeZ === 'string' ? obj.boxSizeZ : '20',
        exhaustiveness: typeof obj.exhaustiveness === 'number' ? obj.exhaustiveness : 8,
        numPoses: typeof obj.numPoses === 'number' ? obj.numPoses : 9,
        isDefault: typeof obj.isDefault === 'boolean' ? obj.isDefault : false,
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.id) {
        throw new Error('Parameter ID is required');
      }

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      // Verify ownership
      const existing = await db
        .select()
        .from(dockingParameters)
        .where(eq(dockingParameters.id, input.id))
        .limit(1);

      if (!existing || existing.length === 0 || existing[0].userId !== ctx.user.id.toString()) {
        throw new Error('Parameter not found or unauthorized');
      }

      await db
        .update(dockingParameters)
        .set({
          name: input.name,
          boxCenterX: input.boxCenterX,
          boxCenterY: input.boxCenterY,
          boxCenterZ: input.boxCenterZ,
          boxSizeX: input.boxSizeX,
          boxSizeY: input.boxSizeY,
          boxSizeZ: input.boxSizeZ,
          exhaustiveness: input.exhaustiveness,
          numPoses: input.numPoses,
          isDefault: input.isDefault ? 1 : 0,
        })
        .where(eq(dockingParameters.id, input.id));

      return { success: true };
    }),

  /**
   * Delete docking parameters
   */
  delete: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { id: 0 };
      const obj = val as Record<string, unknown>;
      return {
        id: typeof obj.id === 'number' ? obj.id : 0,
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      if (!input.id) {
        throw new Error('Parameter ID is required');
      }

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      // Verify ownership
      const existing = await db
        .select()
        .from(dockingParameters)
        .where(eq(dockingParameters.id, input.id))
        .limit(1);

      if (!existing || existing.length === 0 || existing[0].userId !== ctx.user.id.toString()) {
        throw new Error('Parameter not found or unauthorized');
      }

      await db.delete(dockingParameters).where(eq(dockingParameters.id, input.id));

      return { success: true };
    }),
});
