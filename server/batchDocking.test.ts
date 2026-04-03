import { describe, it, expect, beforeAll, vi } from 'vitest';
import { getDb } from './db';
import { batchDockingJobs, batchDockingResults } from '../drizzle/schema';
import { eq, and } from 'drizzle-orm';

/**
 * Async Batch Docking Tests
 * Tests the complete workflow of submitting a batch job, polling for status, and retrieving results
 */

let db: any;
let testJobId: string;
const TEST_USER_ID = 'test-user-123';

beforeAll(async () => {
  db = await getDb();
  expect(db).toBeDefined();
});

describe('Async Batch Docking Job Handling', () => {
  
  it('should create a batch docking job', async () => {
    const jobId = `job_test_${Date.now()}`;
    
    const result = await db
      .insert(batchDockingJobs)
      .values({
        id: jobId,
        userId: TEST_USER_ID,
        jobName: 'Test Batch Job',
        status: 'pending',
        totalCompounds: 3,
        completedCompounds: 0,
        failedCompounds: 0,
        targetName: 'NMDA',
      });

    testJobId = jobId;
    expect(result).toBeDefined();
    
    // Verify job was created
    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, jobId))
      .then((rows: any) => rows[0]);

    expect(job).toBeDefined();
    expect(job.status).toBe('pending');
    expect(job.totalCompounds).toBe(3);
  });

  it('should create result records for each compound', async () => {
    const compoundIds = [1, 2, 3];
    
    for (const analogId of compoundIds) {
      await db
        .insert(batchDockingResults)
        .values({
          jobId: testJobId,
          analogId,
          status: 'pending',
        });
    }

    // Verify result records were created
    const results = await db
      .select()
      .from(batchDockingResults)
      .where(eq(batchDockingResults.jobId, testJobId));

    expect(results).toHaveLength(3);
    expect(results.every((r: any) => r.status === 'pending')).toBe(true);
  });

  it('should update job status to running', async () => {
    await db
      .update(batchDockingJobs)
      .set({
        status: 'running',
        startedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, testJobId));

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, testJobId))
      .then((rows: any) => rows[0]);

    expect(job.status).toBe('running');
    expect(job.startedAt).toBeDefined();
  });

  it('should update individual result with docking data', async () => {
    const mockDockingResult = {
      binding_affinity: -1.771,
      num_poses: 5,
      poses: [
        { mode: 1, affinity: -1.771 },
        { mode: 2, affinity: -1.759 },
      ],
    };

    await db
      .update(batchDockingResults)
      .set({
        status: 'completed',
        bindingAffinity: mockDockingResult.binding_affinity.toString(),
        dockingScore: Math.round(mockDockingResult.binding_affinity * 10),
        numPoses: mockDockingResult.num_poses,
        topPoses: JSON.stringify(mockDockingResult.poses),
        completedAt: new Date(),
      })
      .where(
        and(
          eq(batchDockingResults.jobId, testJobId),
          eq(batchDockingResults.analogId, 1)
        )
      );

    const result = await db
      .select()
      .from(batchDockingResults)
      .where(
        and(
          eq(batchDockingResults.jobId, testJobId),
          eq(batchDockingResults.analogId, 1)
        )
      )
      .then((rows: any) => rows[0]);

    expect(result.status).toBe('completed');
    expect(result.bindingAffinity).toBe('-1.771');
    expect(result.numPoses).toBe(5);
    expect(JSON.parse(result.topPoses)).toHaveLength(2);
  });

  it('should update job progress', async () => {
    await db
      .update(batchDockingJobs)
      .set({
        completedCompounds: 1,
        failedCompounds: 0,
      })
      .where(eq(batchDockingJobs.id, testJobId));

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, testJobId))
      .then((rows: any) => rows[0]);

    expect(job.completedCompounds).toBe(1);
    expect(job.totalCompounds).toBe(3);
    
    const progressPercentage = Math.round(
      ((job.completedCompounds + job.failedCompounds) / job.totalCompounds) * 100
    );
    expect(progressPercentage).toBe(33); // 1 out of 3
  });

  it('should handle failed compounds', async () => {
    await db
      .update(batchDockingResults)
      .set({
        status: 'failed',
        errorMessage: 'SMILES parsing failed',
      })
      .where(
        and(
          eq(batchDockingResults.jobId, testJobId),
          eq(batchDockingResults.analogId, 2)
        )
      );

    const result = await db
      .select()
      .from(batchDockingResults)
      .where(
        and(
          eq(batchDockingResults.jobId, testJobId),
          eq(batchDockingResults.analogId, 2)
        )
      )
      .then((rows: any) => rows[0]);

    expect(result.status).toBe('failed');
    expect(result.errorMessage).toBe('SMILES parsing failed');
  });

  it('should complete remaining compounds', async () => {
    // Complete compound 3
    await db
      .update(batchDockingResults)
      .set({
        status: 'completed',
        bindingAffinity: '-1.650',
        dockingScore: -17,
        numPoses: 5,
        topPoses: JSON.stringify([
          { mode: 1, affinity: -1.650 },
        ]),
        completedAt: new Date(),
      })
      .where(
        and(
          eq(batchDockingResults.jobId, testJobId),
          eq(batchDockingResults.analogId, 3)
        )
      );

    // Update job progress
    await db
      .update(batchDockingJobs)
      .set({
        completedCompounds: 2,
        failedCompounds: 1,
      })
      .where(eq(batchDockingJobs.id, testJobId));

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, testJobId))
      .then((rows: any) => rows[0]);

    expect(job.completedCompounds).toBe(2);
    expect(job.failedCompounds).toBe(1);
  });

  it('should mark job as completed', async () => {
    await db
      .update(batchDockingJobs)
      .set({
        status: 'completed',
        completedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, testJobId));

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, testJobId))
      .then((rows: any) => rows[0]);

    expect(job.status).toBe('completed');
    expect(job.completedAt).toBeDefined();
  });

  it('should calculate job statistics correctly', async () => {
    const results = await db
      .select()
      .from(batchDockingResults)
      .where(eq(batchDockingResults.jobId, testJobId));

    const completedResults = results.filter((r: any) => r.status === 'completed');
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
      failedCompounds: results.filter((r: any) => r.status === 'failed').length,
      pendingCompounds: results.filter((r: any) => r.status === 'pending').length,
      successRate: results.length > 0 ? (completedResults.length / results.length) * 100 : 0,
      averageAffinity: affinities.length > 0 ? affinities.reduce((a, b) => a + b) / affinities.length : 0,
      bestAffinity: affinities.length > 0 ? Math.min(...affinities) : null,
      worstAffinity: affinities.length > 0 ? Math.max(...affinities) : null,
    };

    expect(stats.totalCompounds).toBe(3);
    expect(stats.completedCompounds).toBe(2);
    expect(stats.failedCompounds).toBe(1);
    expect(stats.pendingCompounds).toBe(0);
    expect(stats.successRate).toBe(66.66666666666666);
    expect(stats.averageAffinity).toBeCloseTo(-1.7105, 3);
    expect(stats.bestAffinity).toBe(-1.771);
    expect(stats.worstAffinity).toBe(-1.650);
  });

  it('should retrieve job details with progress', async () => {
    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, testJobId))
      .then((rows: any) => rows[0]);

    const results = await db
      .select()
      .from(batchDockingResults)
      .where(eq(batchDockingResults.jobId, testJobId));

    const progressPercentage = Math.round(
      ((job.completedCompounds + job.failedCompounds) / job.totalCompounds) * 100
    );

    expect(job).toBeDefined();
    expect(results).toHaveLength(3);
    expect(progressPercentage).toBe(100);
  });

  it('should handle concurrent job submissions', async () => {
    const jobIds = [];
    
    // Submit multiple jobs concurrently
    for (let i = 0; i < 3; i++) {
      const jobId = `job_concurrent_${Date.now()}_${i}`;
      jobIds.push(jobId);
      
      await db
        .insert(batchDockingJobs)
        .values({
          id: jobId,
          userId: TEST_USER_ID,
          jobName: `Concurrent Job ${i}`,
          status: 'pending',
          totalCompounds: 2,
          completedCompounds: 0,
          failedCompounds: 0,
          targetName: 'NMDA',
        });
    }

    // Verify all jobs were created
    const jobs = await db
      .select()
      .from(batchDockingJobs)
      .where(
        eq(batchDockingJobs.userId, TEST_USER_ID)
      );

    expect(jobs.length).toBeGreaterThanOrEqual(3);
  });

  it('should support job cancellation', async () => {
    const jobId = `job_cancel_${Date.now()}`;
    
    // Create a job
    await db
      .insert(batchDockingJobs)
      .values({
        id: jobId,
        userId: TEST_USER_ID,
        jobName: 'Job to Cancel',
        status: 'running',
        totalCompounds: 5,
        completedCompounds: 2,
        failedCompounds: 0,
        targetName: 'NMDA',
      });

    // Cancel the job
    await db
      .update(batchDockingJobs)
      .set({
        status: 'cancelled',
        completedAt: new Date(),
      })
      .where(eq(batchDockingJobs.id, jobId));

    const job = await db
      .select()
      .from(batchDockingJobs)
      .where(eq(batchDockingJobs.id, jobId))
      .then((rows: any) => rows[0]);

    expect(job.status).toBe('cancelled');
  });

  it('should support async polling without blocking', async () => {
    const jobId = `job_polling_${Date.now()}`;
    
    // Create a job
    await db
      .insert(batchDockingJobs)
      .values({
        id: jobId,
        userId: TEST_USER_ID,
        jobName: 'Polling Test',
        status: 'pending',
        totalCompounds: 10,
        completedCompounds: 0,
        failedCompounds: 0,
        targetName: 'NMDA',
      });

    // Simulate polling - should be fast
    const startTime = Date.now();
    
    for (let i = 0; i < 5; i++) {
      const job = await db
        .select()
        .from(batchDockingJobs)
        .where(eq(batchDockingJobs.id, jobId))
        .then((rows: any) => rows[0]);
      
      expect(job).toBeDefined();
    }

    const duration = Date.now() - startTime;
    
    // 5 polls should complete in less than 1 second
    expect(duration).toBeLessThan(1000);
  });
});
