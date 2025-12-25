/**
 * Test Platform API authentication and endpoints
 */

import { describe, it, expect } from "vitest";
// Using native fetch (Node 18+)

const API_BASE_URL = "http://localhost:3000";
const API_KEY = process.env.PLATFORM_API_KEY || "";

describe("Platform API", () => {
  it("should have PLATFORM_API_KEY environment variable set", () => {
    expect(API_KEY).toBeTruthy();
    expect(API_KEY.length).toBeGreaterThan(10);
  });

  it("should return 401 for requests without API key", async () => {
    const response = await fetch(`${API_BASE_URL}/api/platform/discoveries/recent`);
    const data = await response.json();
    
    expect(response.status).toBe(401);
    expect(data.error).toBe("Unauthorized");
  });

  it("should accept requests with valid API key", async () => {
    const response = await fetch(
      `${API_BASE_URL}/api/platform/discoveries/recent?apiKey=${API_KEY}&limit=5`
    );
    const data = await response.json();
    
    expect(response.status).toBe(200);
    expect(data.success).toBe(true);
    expect(data).toHaveProperty("discoveries");
  });

  it("should successfully import a test discovery", async () => {
    const testDiscovery = {
      compoundId: `TEST-${Date.now()}`,
      compoundName: "Test Compound",
      smiles: "CC(C)C",
      parentCompound: "Test Parent",
      mechanismOfAction: "Test mechanism",
      keyDifferences: "Test differences",
      confidence: 85,
      similarity: 90,
      safetyScore: 80,
      efficacyScore: 85,
      drugLikenessScore: 90,
      patentStatus: "patent-free",
      marketValue: "$10M",
      discoveryMethod: "test"
    };

    const response = await fetch(`${API_BASE_URL}/api/platform/discoveries/import`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        apiKey: API_KEY,
        discoveries: [testDiscovery]
      })
    });

    const data = await response.json();
    
    expect(response.status).toBe(200);
    expect(data.success).toBe(true);
    expect(data.imported).toBeGreaterThanOrEqual(0);
  });

  it("should have working health check endpoint", async () => {
    const response = await fetch(`${API_BASE_URL}/api/platform/health`);
    const data = await response.json();
    
    expect(response.status).toBe(200);
    expect(data.status).toBe("ok");
    expect(data.service).toBe("pharmasight-admin-dashboard");
  });
});
