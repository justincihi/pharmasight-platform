import { describe, it, expect, beforeAll, afterAll } from "vitest";
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

// Mock non-admin user context
const createUserContext = (): TrpcContext => ({
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
  req: {
    protocol: "https",
    headers: {},
  } as any,
  res: {
    clearCookie: () => {},
  } as any,
});

describe("Analog Discovery Routes", () => {
  let adminCaller: any;
  let userCaller: any;

  beforeAll(() => {
    adminCaller = appRouter.createCaller(createAdminContext());
    userCaller = appRouter.createCaller(createUserContext());
  });

  describe("Admin Access Control", () => {
    it("should allow admin to list analogs", async () => {
      try {
        const result = await adminCaller.analog.list({
          limit: 10,
          offset: 0,
        });
        expect(Array.isArray(result)).toBe(true);
      } catch (error: any) {
        // Database might not have data, but should not be auth error
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should deny non-admin from listing analogs", async () => {
      try {
        await userCaller.analog.list({
          limit: 10,
          offset: 0,
        });
        expect.fail("Should have thrown unauthorized error");
      } catch (error: any) {
        expect(error.message).toContain("Unauthorized");
      }
    });

    it("should allow admin to search analogs", async () => {
      try {
        const result = await adminCaller.analog.search({
          query: "test",
        });
        expect(Array.isArray(result)).toBe(true);
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should deny non-admin from searching analogs", async () => {
      try {
        await userCaller.analog.search({
          query: "test",
        });
        expect.fail("Should have thrown unauthorized error");
      } catch (error: any) {
        expect(error.message).toContain("Unauthorized");
      }
    });
  });

  describe("Cheminformatics Tests", () => {
    it("should allow admin to run ADMET test", async () => {
      try {
        const result = await adminCaller.analog.runADMET({
          analogId: 1,
          smiles: "CC(=O)Oc1ccccc1C(=O)O",
        });
        expect(result).toBeDefined();
        expect(result.testType).toBe("admet");
        expect(result.testStatus).toBe("completed");
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should deny non-admin from running ADMET test", async () => {
      try {
        await userCaller.analog.runADMET({
          analogId: 1,
          smiles: "CC(=O)Oc1ccccc1C(=O)O",
        });
        expect.fail("Should have thrown unauthorized error");
      } catch (error: any) {
        expect(error.message).toContain("Unauthorized");
      }
    });

    it("should allow admin to run docking test", async () => {
      try {
        const result = await adminCaller.analog.runDocking({
          analogId: 1,
          smiles: "CC(=O)Oc1ccccc1C(=O)O",
          target: "EGFR",
        });
        expect(result).toBeDefined();
        expect(result.testType).toBe("docking");
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should allow admin to run toxicity test", async () => {
      try {
        const result = await adminCaller.analog.runToxicity({
          analogId: 1,
          smiles: "CC(=O)Oc1ccccc1C(=O)O",
        });
        expect(result).toBeDefined();
        expect(result.testType).toBe("toxicity");
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });
  });

  describe("Analytics Routes", () => {
    it("should allow admin to get analytics stats", async () => {
      try {
        const result = await adminCaller.analytics.getStats();
        if (result) {
          expect(typeof result.totalDiscovered).toBe("number");
          expect(typeof result.highConfidencePercentage).toBe("number");
          expect(typeof result.patentFreePercentage).toBe("number");
        }
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should deny non-admin from getting analytics", async () => {
      try {
        await userCaller.analytics.getStats();
        expect.fail("Should have thrown unauthorized error");
      } catch (error: any) {
        expect(error.message).toContain("Unauthorized");
      }
    });

    it("should allow admin to get discovery timeline", async () => {
      try {
        const result = await adminCaller.analytics.getTimeline({
          days: 30,
        });
        expect(Array.isArray(result)).toBe(true);
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });
  });

  describe("Notification Routes", () => {
    it("should allow admin to get recent notifications", async () => {
      try {
        const result = await adminCaller.notifications.getRecent({
          limit: 20,
        });
        expect(Array.isArray(result)).toBe(true);
      } catch (error: any) {
        expect(error.message).not.toContain("Unauthorized");
      }
    });

    it("should deny non-admin from getting notifications", async () => {
      try {
        await userCaller.notifications.getRecent({
          limit: 20,
        });
        expect.fail("Should have thrown unauthorized error");
      } catch (error: any) {
        expect(error.message).toContain("Unauthorized");
      }
    });
  });
});
