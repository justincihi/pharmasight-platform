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

describe("BioTransformer Integration", () => {
  let adminCaller: any;

  beforeAll(() => {
    const adminCtx = createAdminContext();
    adminCaller = appRouter.createCaller(adminCtx);
  });

  it("should get available metabolism types", async () => {
    const result = await adminCaller.biotransformer.getMetabolismTypes();

    expect(result).toBeDefined();
    expect(result).toHaveProperty('human');
    expect(result).toHaveProperty('cyp450');
    expect(result).toHaveProperty('phase2');
    expect(result).toHaveProperty('gut');
  });

  it("should predict metabolites for ketamine", async () => {
    const ketamineSmiles = "CCN(C1CCCCC1=O)c2cccc(F)c2Cl";

    const result = await adminCaller.biotransformer.predictMetabolites({
      smiles: ketamineSmiles,
      metabolismType: 'human',
      steps: 1,
    });

    expect(result).toBeDefined();
    expect(result.success).toBe(true);
    expect(result.parent_smiles).toBe(ketamineSmiles);
    expect(result.metabolites).toBeDefined();
    expect(Array.isArray(result.metabolites)).toBe(true);
    expect(result.metabolites.length).toBeGreaterThan(0);

    // Check first metabolite structure
    const firstMetabolite = result.metabolites[0];
    expect(firstMetabolite).toHaveProperty('smiles');
    expect(firstMetabolite).toHaveProperty('reaction');
    expect(firstMetabolite).toHaveProperty('enzyme');
  });

  it("should predict Phase II metabolites with steps=2", async () => {
    const ketamineSmiles = "CCN(C1CCCCC1=O)c2cccc(F)c2Cl";

    const result = await adminCaller.biotransformer.predictMetabolites({
      smiles: ketamineSmiles,
      metabolismType: 'human',
      steps: 2,
    });

    expect(result).toBeDefined();
    expect(result.success).toBe(true);
    expect(result.steps).toBe(2);
    expect(result.metabolites).toBeDefined();

    // Should have both Phase I and Phase II metabolites
    const hasPhaseI = result.metabolites.some((m: any) =>
      m.reaction === 'Hydroxylation' || m.reaction === 'N-oxidation'
    );
    const hasPhaseII = result.metabolites.some((m: any) =>
      m.reaction === 'Glucuronidation' || m.reaction === 'Sulfation'
    );

    expect(hasPhaseI).toBe(true);
    expect(hasPhaseII).toBe(true);
  });

  it("should handle batch predictions", async () => {
    const compounds = [
      { smiles: "CCN(C1CCCCC1=O)c2cccc(F)c2Cl", id: "ketamine" },
      { smiles: "CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3", id: "atropine" },
    ];

    const result = await adminCaller.biotransformer.batchPredict({
      compounds,
      metabolismType: 'human',
      steps: 1,
    });

    expect(result).toBeDefined();
    expect(Array.isArray(result)).toBe(true);
    expect(result.length).toBe(2);

    result.forEach((r: any) => {
      expect(r.success).toBe(true);
      expect(r).toHaveProperty('id');
      expect(r.metabolites).toBeDefined();
    });
  });

  it("should predict gut microbiome metabolism", async () => {
    const testSmiles = "CCN(C1CCCCC1=O)c2cccc(F)c2Cl";

    const result = await adminCaller.biotransformer.predictMetabolites({
      smiles: testSmiles,
      metabolismType: 'gut',
      steps: 1,
    });

    expect(result).toBeDefined();
    expect(result.success).toBe(true);
    expect(result.metabolism_type).toBe('gut');

    // Check for gut-specific transformation
    const hasGutReaction = result.metabolites.some((m: any) =>
      m.reaction.includes('dehalogenation') || m.enzyme.includes('microbiome')
    );
    expect(hasGutReaction).toBe(true);
  });

  it("should indicate mock mode when BioTransformer JAR is not available", async () => {
    const result = await adminCaller.biotransformer.predictMetabolites({
      smiles: "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",
      metabolismType: 'human',
      steps: 1,
    });

    expect(result).toBeDefined();
    expect(result).toHaveProperty('mock_mode');

    if (result.mock_mode === true) {
      expect(result).toHaveProperty('note');
      expect(result.note).toContain('mock data');
    }
  });
});
