import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { createTRPCMsw } from 'msw-trpc';
import { appRouter } from './routers';

describe('Docking Parameters Router', () => {
  const mockUser = {
    id: 1,
    email: 'test@example.com',
    role: 'admin' as const,
  };

  const mockCtx = {
    user: mockUser,
    req: {} as any,
    res: {} as any,
  };

  describe('create', () => {
    it('should create a new docking parameter configuration', async () => {
      const caller = appRouter.createCaller(mockCtx);

      const result = await caller.dockingParams.create({
        name: 'NMDA-Test',
        targetName: 'NMDA Receptor',
        boxCenterX: '0.0',
        boxCenterY: '0.0',
        boxCenterZ: '0.0',
        boxSizeX: '25',
        boxSizeY: '25',
        boxSizeZ: '25',
        exhaustiveness: 16,
        numPoses: 9,
        isDefault: false,
      });

      expect(result).toBeDefined();
      expect(result.success).toBe(true);
    });

    it('should validate exhaustiveness range', async () => {
      const caller = appRouter.createCaller(mockCtx);

      try {
        await caller.dockingParams.create({
          name: 'Invalid-Exhaustiveness',
          targetName: 'NMDA Receptor',
          boxCenterX: '0',
          boxCenterY: '0',
          boxCenterZ: '0',
          boxSizeX: '20',
          boxSizeY: '20',
          boxSizeZ: '20',
          exhaustiveness: 50, // Invalid: > 32
          numPoses: 9,
          isDefault: false,
        });
        expect.fail('Should have thrown error for invalid exhaustiveness');
      } catch (error: any) {
        expect(error.message).toContain('between 1 and 32');
      }
    });

    it('should validate number of poses range', async () => {
      const caller = appRouter.createCaller(mockCtx);

      try {
        await caller.dockingParams.create({
          name: 'Invalid-Poses',
          targetName: 'NMDA Receptor',
          boxCenterX: '0',
          boxCenterY: '0',
          boxCenterZ: '0',
          boxSizeX: '20',
          boxSizeY: '20',
          boxSizeZ: '20',
          exhaustiveness: 8,
          numPoses: 25, // Invalid: > 20
          isDefault: false,
        });
        expect.fail('Should have thrown error for invalid poses');
      } catch (error: any) {
        expect(error.message).toContain('between 1 and 20');
      }
    });

    it('should require name and target name', async () => {
      const caller = appRouter.createCaller(mockCtx);

      try {
        await caller.dockingParams.create({
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
        });
        expect.fail('Should have thrown error for missing name/target');
      } catch (error: any) {
        expect(error.message).toContain('required');
      }
    });
  });

  describe('list', () => {
    it('should list parameters for a specific target', async () => {
      const caller = appRouter.createCaller(mockCtx);

      // Create a parameter first
      await caller.dockingParams.create({
        name: 'Test-List',
        targetName: '5HT2A Receptor',
        boxCenterX: '5',
        boxCenterY: '5',
        boxCenterZ: '5',
        boxSizeX: '20',
        boxSizeY: '20',
        boxSizeZ: '20',
        exhaustiveness: 12,
        numPoses: 9,
        isDefault: false,
      });

      // List parameters
      const result = await caller.dockingParams.list({ targetName: '5HT2A Receptor' });
      expect(Array.isArray(result)).toBe(true);
    });
  });

  describe('getDefault', () => {
    it('should get default parameters for a target', async () => {
      const caller = appRouter.createCaller(mockCtx);

      // Create a default parameter
      await caller.dockingParams.create({
        name: 'Default-Config',
        targetName: 'Dopamine D2 Receptor',
        boxCenterX: '0',
        boxCenterY: '0',
        boxCenterZ: '0',
        boxSizeX: '22',
        boxSizeY: '22',
        boxSizeZ: '22',
        exhaustiveness: 14,
        numPoses: 9,
        isDefault: true,
      });

      // Get default
      const result = await caller.dockingParams.getDefault({
        targetName: 'Dopamine D2 Receptor',
      });

      expect(result).toBeDefined();
      if (result) {
        expect(result.isDefault).toBe(1);
        expect(result.targetName).toBe('Dopamine D2 Receptor');
      }
    });
  });

  describe('authorization', () => {
    it('should reject non-admin users', async () => {
      const nonAdminCtx = {
        user: { ...mockUser, role: 'user' as const },
        req: {} as any,
        res: {} as any,
      };

      const caller = appRouter.createCaller(nonAdminCtx);

      try {
        await caller.dockingParams.create({
          name: 'Unauthorized',
          targetName: 'NMDA Receptor',
          boxCenterX: '0',
          boxCenterY: '0',
          boxCenterZ: '0',
          boxSizeX: '20',
          boxSizeY: '20',
          boxSizeZ: '20',
          exhaustiveness: 8,
          numPoses: 9,
          isDefault: false,
        });
        expect.fail('Should have thrown authorization error');
      } catch (error: any) {
        expect(error.message).toContain('Admin access required');
      }
    });
  });
});
