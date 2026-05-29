import { describe, it, expect } from 'vitest';

describe('LLM API Keys Validation', () => {
  it('SONAR_API_KEY is set in environment', () => {
    expect(process.env.SONAR_API_KEY).toBeDefined();
    expect(process.env.SONAR_API_KEY!.length).toBeGreaterThan(10);
  });

  it('GEMINI_API_KEY is set in environment', () => {
    expect(process.env.GEMINI_API_KEY).toBeDefined();
    expect(process.env.GEMINI_API_KEY!.length).toBeGreaterThan(10);
  });

  it('can reach Perplexity API endpoint', async () => {
    const apiKey = process.env.SONAR_API_KEY;
    if (!apiKey) return; // skip if not set

    const response = await fetch('https://api.perplexity.ai/chat/completions', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${apiKey}`,
      },
      body: JSON.stringify({
        model: 'sonar',
        messages: [{ role: 'user', content: 'Say "ok" in one word.' }],
        max_tokens: 5,
      }),
    });

    // 200 = success, 401 = invalid key, 429 = rate limited (key is valid)
    expect([200, 429]).toContain(response.status);
  }, 15000);
});
