/**
 * Unit tests for multiLLM provider configuration and routing logic.
 * These tests run without real API keys — they verify the provider detection
 * and routing logic, not the actual LLM calls.
 */

import { describe, it, expect, afterEach, vi } from "vitest";
import { getAvailableProviders, getProvidersStatus } from "./multiLLM";

// Helper to cleanly set/restore env vars for each test
const setEnv = (vars: Record<string, string | undefined>) => {
  const originals: Record<string, string | undefined> = {};
  for (const [key, value] of Object.entries(vars)) {
    originals[key] = process.env[key];
    if (value === undefined) {
      delete process.env[key];
    } else {
      process.env[key] = value;
    }
  }
  return () => {
    for (const [key, original] of Object.entries(originals)) {
      if (original === undefined) {
        delete process.env[key];
      } else {
        process.env[key] = original;
      }
    }
  };
};

describe("getAvailableProviders", () => {
  afterEach(() => {
    // Clean up any keys the test may have set
    delete process.env.GEMINI_API_KEY;
    delete process.env.ANTHROPIC_API_KEY;
    delete process.env.PERPLEXITY_API_KEY;
    delete process.env.SONAR_API_KEY;
    delete process.env.XAI_API_KEY;
  });

  it("always includes openai", () => {
    const providers = getAvailableProviders();
    expect(providers).toContain("openai");
  });

  it("includes gemini when GEMINI_API_KEY is set", () => {
    const restore = setEnv({ GEMINI_API_KEY: "test-key" });
    try {
      expect(getAvailableProviders()).toContain("gemini");
    } finally {
      restore();
    }
  });

  it("does not include gemini when GEMINI_API_KEY is absent", () => {
    const restore = setEnv({ GEMINI_API_KEY: undefined });
    try {
      expect(getAvailableProviders()).not.toContain("gemini");
    } finally {
      restore();
    }
  });

  it("includes claude when ANTHROPIC_API_KEY is set", () => {
    const restore = setEnv({ ANTHROPIC_API_KEY: "test-key" });
    try {
      expect(getAvailableProviders()).toContain("claude");
    } finally {
      restore();
    }
  });

  it("includes perplexity when PERPLEXITY_API_KEY is set", () => {
    const restore = setEnv({ PERPLEXITY_API_KEY: "test-key", SONAR_API_KEY: undefined });
    try {
      expect(getAvailableProviders()).toContain("perplexity");
    } finally {
      restore();
    }
  });

  it("includes perplexity when legacy SONAR_API_KEY is set (backwards compat)", () => {
    const restore = setEnv({ PERPLEXITY_API_KEY: undefined, SONAR_API_KEY: "legacy-key" });
    try {
      expect(getAvailableProviders()).toContain("perplexity");
    } finally {
      restore();
    }
  });

  it("does not include perplexity when neither PERPLEXITY_API_KEY nor SONAR_API_KEY is set", () => {
    const restore = setEnv({ PERPLEXITY_API_KEY: undefined, SONAR_API_KEY: undefined });
    try {
      expect(getAvailableProviders()).not.toContain("perplexity");
    } finally {
      restore();
    }
  });

  it("includes xai when XAI_API_KEY is set", () => {
    const restore = setEnv({ XAI_API_KEY: "test-key" });
    try {
      expect(getAvailableProviders()).toContain("xai");
    } finally {
      restore();
    }
  });

  it("does not include xai when XAI_API_KEY is absent", () => {
    const restore = setEnv({ XAI_API_KEY: undefined });
    try {
      expect(getAvailableProviders()).not.toContain("xai");
    } finally {
      restore();
    }
  });

  it("returns all five providers when every key is set", () => {
    const restore = setEnv({
      GEMINI_API_KEY: "g",
      ANTHROPIC_API_KEY: "a",
      PERPLEXITY_API_KEY: "p",
      XAI_API_KEY: "x",
    });
    try {
      const providers = getAvailableProviders();
      expect(providers).toContain("openai");
      expect(providers).toContain("gemini");
      expect(providers).toContain("claude");
      expect(providers).toContain("perplexity");
      expect(providers).toContain("xai");
      expect(providers).toHaveLength(5);
    } finally {
      restore();
    }
  });
});

describe("getProvidersStatus", () => {
  it("returns a status entry for every supported provider", () => {
    const status = getProvidersStatus();
    expect(status).toHaveProperty("openai");
    expect(status).toHaveProperty("gemini");
    expect(status).toHaveProperty("claude");
    expect(status).toHaveProperty("perplexity");
    expect(status).toHaveProperty("xai");
  });

  it("each status entry has configured boolean and keyVar string", () => {
    const status = getProvidersStatus();
    for (const entry of Object.values(status)) {
      expect(typeof entry.configured).toBe("boolean");
      expect(typeof entry.keyVar).toBe("string");
      expect(entry.keyVar.length).toBeGreaterThan(0);
    }
  });

  it("reports perplexity as configured when PERPLEXITY_API_KEY is set", () => {
    const restore = setEnv({ PERPLEXITY_API_KEY: "p-key", SONAR_API_KEY: undefined });
    try {
      expect(getProvidersStatus().perplexity.configured).toBe(true);
    } finally {
      restore();
    }
  });

  it("reports perplexity as configured when legacy SONAR_API_KEY is set", () => {
    const restore = setEnv({ PERPLEXITY_API_KEY: undefined, SONAR_API_KEY: "s-key" });
    try {
      expect(getProvidersStatus().perplexity.configured).toBe(true);
    } finally {
      restore();
    }
  });

  it("reports xai as not configured when XAI_API_KEY is absent", () => {
    const restore = setEnv({ XAI_API_KEY: undefined });
    try {
      expect(getProvidersStatus().xai.configured).toBe(false);
    } finally {
      restore();
    }
  });

  it("reports xai as configured when XAI_API_KEY is set", () => {
    const restore = setEnv({ XAI_API_KEY: "grok-key" });
    try {
      expect(getProvidersStatus().xai.configured).toBe(true);
    } finally {
      restore();
    }
  });

  it("never leaks actual key values", () => {
    const restore = setEnv({ OPENAI_API_KEY: "secret-openai-key" });
    try {
      const status = getProvidersStatus();
      const statusJson = JSON.stringify(status);
      expect(statusJson).not.toContain("secret-openai-key");
    } finally {
      restore();
    }
  });
});
