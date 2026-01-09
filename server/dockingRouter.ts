import { protectedProcedure, router } from './_core/trpc';
import { getDockingQueueStatus, getAnalogDockingResults, queueAnalogForDocking } from './dockingQueueManager';
import { getDb } from './db';
import { dockingQueue } from '../drizzle/schema';
import { eq, and } from 'drizzle-orm';

export const dockingRouter = router({
  enqueue: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { analogId: 0, priority: 5 };
      const obj = val as Record<string, unknown>;
      return {
        analogId: typeof obj.analogId === 'number' ? obj.analogId : 0,
        priority: typeof obj.priority === 'number' ? obj.priority : 5,
      };
    })
    .mutation(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      await queueAnalogForDocking(input.analogId, input.priority);
      return { success: true, message: 'Analog queued for docking' };
    }),
  getQueueStatus: protectedProcedure.query(async ({ ctx }) => {
    if (ctx.user?.role !== 'admin') {
      throw new Error('Unauthorized: Admin access required');
    }
    return getDockingQueueStatus();
  }),

  getQueueJobs: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return {};
      const obj = val as Record<string, unknown>;
      return {
        target: typeof obj.target === 'string' ? obj.target : undefined,
        status: typeof obj.status === 'string' ? obj.status : undefined,
      };
    })
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }

      const db = await getDb();
      if (!db) throw new Error('Database not available');

      let query = db.select().from(dockingQueue);

      const conditions = [];
      if (input.target) {
        conditions.push(eq(dockingQueue.target, input.target));
      }
      if (input.status) {
        conditions.push(eq(dockingQueue.status, input.status as any));
      }

      if (conditions.length > 0) {
        query = query.where(and(...conditions)) as any;
      }

      const jobs = await query.limit(100);
      return jobs;
    }),

  getAnalogResults: protectedProcedure
    .input((val: unknown) => {
      if (typeof val !== 'object' || val === null) return { analogId: 0 };
      const obj = val as Record<string, unknown>;
      return { analogId: typeof obj.analogId === 'number' ? obj.analogId : 0 };
    })
    .query(async ({ input, ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      return getAnalogDockingResults(input.analogId);
    }),

  clearCompleted: protectedProcedure.mutation(async ({ ctx }) => {
    if (ctx.user?.role !== 'admin') {
      throw new Error('Unauthorized: Admin access required');
    }

    const db = await getDb();
    if (!db) throw new Error('Database not available');

    await db.delete(dockingQueue).where(eq(dockingQueue.status, 'completed'));

    return { success: true };
  }),

  retryFailed: protectedProcedure.mutation(async ({ ctx }) => {
    if (ctx.user?.role !== 'admin') {
      throw new Error('Unauthorized: Admin access required');
    }

    const db = await getDb();
    if (!db) throw new Error('Database not available');

    await db
      .update(dockingQueue)
      .set({ status: 'pending', errorMessage: null })
      .where(eq(dockingQueue.status, 'failed'));

    return { success: true };
  }),
});
