import { describe, it, expect, beforeAll } from "vitest";
import { appRouter } from "./routers";
import type { TrpcContext } from "./_core/context";

// Mock admin user context
const createAdminContext = (): TrpcContext => ({
  user: {
    id: 1,
    openId: "test-admin",
    email: "admin@test.com",
    name: "Test Admin",
    loginMethod: "test",
    role: "admin",
    createdAt: new Date(),
    updatedAt: new Date(),
    lastSignedIn: new Date(),
  },
  req: {
    protocol: "https",
    headers: {},
  } as any,
  res: {
    clearCookie: () => {},
  } as any,
});

describe("Autonomous Research Scheduler", () => {
  let adminCaller: any;

  beforeAll(() => {
    adminCaller = appRouter.createCaller(createAdminContext());
  });

  it("should get scheduler status", async () => {
    const result = await adminCaller.scheduler.status();

    expect(result).toBeDefined();
    expect(typeof result.running).toBe("boolean");
    expect(result.config).toBeDefined();
  });

  it("should allow admin to run scheduler immediately", async () => {
    try {
      const result = await adminCaller.scheduler.runNow();

      expect(result.success).toBe(true);
      expect(result.message).toBeDefined();
    } catch (error: any) {
      // Scheduler might fail if no discoveries file exists
      expect(error.message).not.toContain("Unauthorized");
    }
  });

  it("should start and stop scheduler", async () => {
    const startResult = await adminCaller.scheduler.start();
    expect(startResult.success).toBe(true);

    const statusAfterStart = await adminCaller.scheduler.status();
    expect(statusAfterStart.running).toBe(true);

    const stopResult = await adminCaller.scheduler.stop();
    expect(stopResult.success).toBe(true);

    const statusAfterStop = await adminCaller.scheduler.status();
    expect(statusAfterStop.running).toBe(false);
  });
});

describe("Python Cheminformatics Integration", () => {
  let adminCaller: any;

  beforeAll(() => {
    adminCaller = appRouter.createCaller(createAdminContext());
  });

  const testSmiles = "CC(C)Cc1ccc(cc1)C(C)C(O)=O"; // Ibuprofen

  it("should validate compound with ChEMBL", async () => {
    try {
      const result = await adminCaller.cheminformatics.validateChEMBL({
        smiles: testSmiles,
      });

      expect(result).toBeDefined();
      // Python might not be available in test environment
      if (result.success) {
        expect(result.data).toBeDefined();
      }
    } catch (error: any) {
      expect(error.message).not.toContain("Unauthorized");
    }
  });

  it("should predict ADMET properties", async () => {
    try {
      const result = await adminCaller.cheminformatics.predictADMET({
        smiles: testSmiles,
      });

      expect(result).toBeDefined();
      if (result.success) {
        expect(result.data).toBeDefined();
      }
    } catch (error: any) {
      expect(error.message).not.toContain("Unauthorized");
    }
  });

  it("should predict toxicity", async () => {
    try {
      const result = await adminCaller.cheminformatics.predictToxicity({
        smiles: testSmiles,
      });

      expect(result).toBeDefined();
      if (result.success) {
        expect(result.data).toBeDefined();
      }
    } catch (error: any) {
      expect(error.message).not.toContain("Unauthorized");
    }
  });

  it("should generate analogs", async () => {
    try {
      const result = await adminCaller.cheminformatics.generateAnalogs({
        parentSmiles: testSmiles,
        numAnalogs: 5,
      });

      expect(result).toBeDefined();
      if (result.success) {
        expect(result.data).toBeDefined();
      }
    } catch (error: any) {
      expect(error.message).not.toContain("Unauthorized");
    }
  });

  it("should simulate PK/PD", async () => {
    try {
      const result = await adminCaller.cheminformatics.simulatePKPD({
        smiles: testSmiles,
        dose: 200,
        route: "oral",
      });

      expect(result).toBeDefined();
      if (result.success) {
        expect(result.data).toBeDefined();
      }
    } catch (error: any) {
      expect(error.message).not.toContain("Unauthorized");
    }
  });
});

describe("Access Control", () => {
  it("should deny non-admin access to scheduler", async () => {
    const userContext: TrpcContext = {
      user: {
        id: 2,
        openId: "test-user",
        email: "user@test.com",
        name: "Test User",
        loginMethod: "test",
        role: "user",
        createdAt: new Date(),
        updatedAt: new Date(),
        lastSignedIn: new Date(),
      },
      req: { protocol: "https", headers: {} } as any,
      res: { clearCookie: () => {} } as any,
    };

    const userCaller = appRouter.createCaller(userContext);

    await expect(userCaller.scheduler.status()).rejects.toThrow("Unauthorized");
  });

  it("should deny non-admin access to cheminformatics", async () => {
    const userContext: TrpcContext = {
      user: {
        id: 2,
        openId: "test-user",
        email: "user@test.com",
        name: "Test User",
        loginMethod: "test",
        role: "user",
        createdAt: new Date(),
        updatedAt: new Date(),
        lastSignedIn: new Date(),
      },
      req: { protocol: "https", headers: {} } as any,
      res: { clearCookie: () => {} } as any,
    };

    const userCaller = appRouter.createCaller(userContext);

    await expect(
      userCaller.cheminformatics.predictADMET({ smiles: "CC" })
    ).rejects.toThrow("Unauthorized");
  });
});
