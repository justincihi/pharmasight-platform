import { readFile, writeFile, existsSync } from "fs";
import { promisify } from "util";

const readFileAsync = promisify(readFile);
const writeFileAsync = promisify(writeFile);

const TRENDS_FILE = "/home/ubuntu/pharmasight-admin-dashboard/data/medical_trends.json";
const SONAR_API_KEY = process.env.SONAR_API_KEY;
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;

interface MedicalTrend {
  title: string;
  description: string;
  priority: "High" | "Medium" | "Low";
  keywords: string[];
  source?: string;
}

interface MedicalTrendsData {
  trends: MedicalTrend[];
  lastUpdated: string;
}

/**
 * Get cached medical trends
 */
export async function getMedicalTrends(): Promise<MedicalTrendsData> {
  try {
    if (!existsSync(TRENDS_FILE)) {
      return {
        trends: [],
        lastUpdated: new Date().toISOString(),
      };
    }

    const data = await readFileAsync(TRENDS_FILE, "utf-8");
    return JSON.parse(data);
  } catch (error) {
    console.error("[Medical Trends] Error reading trends:", error);
    return {
      trends: [],
      lastUpdated: new Date().toISOString(),
    };
  }
}

const TRENDS_PROMPT = `Analyze recent medical breakthroughs and identify the top 5 most promising research areas for drug discovery. For each trend, provide:
1. A concise title (max 10 words)
2. A brief description (max 50 words) explaining why it's valuable
3. Priority level (High/Medium/Low) based on commercial potential and unmet need
4. 3-5 relevant keywords for research targeting

Format your response as a JSON array with this structure:
[
  {
    "title": "...",
    "description": "...",
    "priority": "High",
    "keywords": ["keyword1", "keyword2", "keyword3"]
  }
]

Focus on areas like: novel therapeutic targets, emerging drug classes, breakthrough mechanisms, underserved conditions, and patent-free opportunities.`;

/**
 * Refresh medical trends using Perplexity or Gemini API
 */
export async function refreshMedicalTrends(): Promise<void> {
  try {
    if (!SONAR_API_KEY && !GEMINI_API_KEY) {
      console.error("[Medical Trends] No API keys configured (need SONAR_API_KEY or GEMINI_API_KEY)");
      throw new Error("No LLM API key configured for medical trends");
    }

    let content: string | undefined;
    let source = "Unknown";

    // Try Perplexity first
    if (SONAR_API_KEY) {
      console.log("[Medical Trends] Trying Perplexity API...");
      try {
        const response = await fetch("https://api.perplexity.ai/chat/completions", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            "Authorization": `Bearer ${SONAR_API_KEY}`,
          },
          body: JSON.stringify({
            model: "sonar-pro",
            messages: [
              {
                role: "system",
                content: "You are a pharmaceutical research analyst. Identify the most promising and groundbreaking medical research trends from the past 6 months that could lead to valuable drug discovery opportunities."
              },
              {
                role: "user",
                content: TRENDS_PROMPT
              }
            ],
            temperature: 0.3,
            max_tokens: 2000,
          }),
        });

        if (response.ok) {
          const result = await response.json();
          content = result.choices?.[0]?.message?.content;
          source = "Perplexity Sonar Pro";
          console.log("[Medical Trends] Successfully fetched from Perplexity");
        } else {
          console.warn(`[Medical Trends] Perplexity failed: ${response.status} ${response.statusText}`);
        }
      } catch (error) {
        console.warn("[Medical Trends] Perplexity error:", error);
      }
    }

    // Fallback to Gemini if Perplexity failed
    if (!content && GEMINI_API_KEY) {
      console.log("[Medical Trends] Falling back to Gemini API...");
      try {
        const response = await fetch(
          `https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash-exp:generateContent?key=${GEMINI_API_KEY}`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              contents: [{
                parts: [{
                  text: `You are a pharmaceutical research analyst. ${TRENDS_PROMPT}`
                }]
              }],
              generationConfig: {
                temperature: 0.3,
                maxOutputTokens: 2000,
              }
            }),
          }
        );

        if (response.ok) {
          const result = await response.json();
          content = result.candidates?.[0]?.content?.parts?.[0]?.text;
          source = "Google Gemini";
          console.log("[Medical Trends] Successfully fetched from Gemini");
        } else {
          console.error(`[Medical Trends] Gemini failed: ${response.status} ${response.statusText}`);
        }
      } catch (error) {
        console.error("[Medical Trends] Gemini error:", error);
      }
    }

    if (!content) {
      throw new Error("Failed to fetch trends from any LLM provider");
    }

    // Extract JSON from response (handle markdown code blocks)
    let trendsData: MedicalTrend[];
    try {
      const jsonMatch = content.match(/```json\n([\s\S]*?)\n```/) || content.match(/```\n([\s\S]*?)\n```/);
      const jsonStr = jsonMatch ? jsonMatch[1] : content;
      trendsData = JSON.parse(jsonStr);
    } catch (parseError) {
      console.error("[Medical Trends] Failed to parse JSON:", content.substring(0, 200));
      throw new Error("Failed to parse trends from AI response");
    }

    // Validate and save trends
    if (!Array.isArray(trendsData)) {
      throw new Error("Invalid trends format");
    }

    const validatedTrends: MedicalTrend[] = trendsData
      .filter(trend => trend.title && trend.description && trend.keywords)
      .map(trend => ({
        title: trend.title,
        description: trend.description,
        priority: trend.priority || "Medium",
        keywords: Array.isArray(trend.keywords) ? trend.keywords : [],
        source,
      }));

    const data: MedicalTrendsData = {
      trends: validatedTrends,
      lastUpdated: new Date().toISOString(),
    };

    await writeFileAsync(TRENDS_FILE, JSON.stringify(data, null, 2));
    console.log(`[Medical Trends] Saved ${validatedTrends.length} trends from ${source}`);
  } catch (error) {
    console.error("[Medical Trends] Error refreshing trends:", error);
    throw error;
  }
}

/**
 * Alias for refreshMedicalTrends for backward compatibility
 */
export const analyzeMedicalTrends = refreshMedicalTrends;
