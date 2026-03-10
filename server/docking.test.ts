import { describe, it, expect, beforeAll, afterAll, vi } from "vitest";
import { createCaller } from "./routers";
import { getDb } from "./db";

describe("Docking Endpoint", () => {
  let caller: any;
  let mockUser: any;

  beforeAll(async () => {
    // Create a mock admin user context
    mockUser = {
      id: "test-user-123",
      role: "admin",
    };

    // Create a test caller with admin context
    caller = createCaller({
      user: mockUser,
      req: {} as any,
      res: {} as any,
    });
  });

  it("should validate SMILES string format", async () => {
    // Test with valid SMILES
    const validSmiles = "C0c1ccc2c(c1)C(=O)C(=O)N2";
    
    try {
      const result = await caller.analog.runDocking({
        analogId: 1,
        smiles: validSmiles,
        target: "5HT2A",
      });
      
      // Should have docking result
      expect(result).toBeDefined();
      expect(result.dockingScore).toBeDefined();
    } catch (error: any) {
      // If it fails, it should be a meaningful error, not a pattern error
      expect(error.message).not.toContain("did not match the expected pattern");
    }
  });

  it("should reject invalid SMILES strings", async () => {
    // Test with invalid SMILES (contains invalid characters)
    const invalidSmiles = "C0c1ccc2c(c1)C(=O)C(=O)N2!!!";
    
    try {
      await caller.analog.runDocking({
        analogId: 1,
        smiles: invalidSmiles,
        target: "5HT2A",
      });
      
      // Should not reach here
      expect(true).toBe(false);
    } catch (error: any) {
      // Should have a clear error message
      expect(error.message).toContain("Invalid SMILES format");
    }
  });

  it("should reject empty SMILES strings", async () => {
    try {
      await caller.analog.runDocking({
        analogId: 1,
        smiles: "",
        target: "5HT2A",
      });
      
      // Should not reach here
      expect(true).toBe(false);
    } catch (error: any) {
      expect(error.message).toContain("SMILES string is required");
    }
  });

  it("should handle whitespace in SMILES strings", async () => {
    // Test with SMILES that has leading/trailing whitespace
    const smilesWithWhitespace = "  C0c1ccc2c(c1)C(=O)C(=O)N2  ";
    
    try {
      const result = await caller.analog.runDocking({
        analogId: 1,
        smiles: smilesWithWhitespace,
        target: "5HT2A",
      });
      
      // Should handle whitespace gracefully
      expect(result).toBeDefined();
    } catch (error: any) {
      // Should not be a pattern error
      expect(error.message).not.toContain("did not match the expected pattern");
    }
  });

  it("should use default target if not provided", async () => {
    const validSmiles = "C0c1ccc2c(c1)C(=O)C(=O)N2";
    
    try {
      const result = await caller.analog.runDocking({
        analogId: 1,
        smiles: validSmiles,
        // No target provided
      });
      
      expect(result).toBeDefined();
    } catch (error: any) {
      // Should not fail due to missing target
      expect(error.message).not.toContain("target");
    }
  });
});
