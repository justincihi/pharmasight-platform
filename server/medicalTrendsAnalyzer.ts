import { readFile, writeFile, existsSync } from "fs";
import { promisify } from "util";

const readFileAsync = promisify(readFile);
const writeFileAsync = promisify(writeFile);

const TRENDS_FILE = "/home/ubuntu/pharmasight-admin-dashboard/data/medical_trends.json";
const SONAR_API_KEY = process.env.SONAR_API_KEY;

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

/**
 * Refresh medical trends using Perplexity API
 */
export async function refreshMedicalTrends(): Promise<void> {
  try {
    if (!SONAR_API_KEY) {
      console.error("[Medical Trends] SONAR_API_KEY not configured");
      throw new Error("Perplexity API key not configured");
    }

    console.log("[Medical Trends] Fetching latest breakthroughs from Perplexity...");

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
            content: "You are a pharmaceutical research analyst. Identify the most promising and groundbreaking medical research trends from the past 6 months that could lead to valuable drug discovery opportunities. Focus on areas with high commercial potential and unmet medical needs."
          },
          {
            role: "user",
            content: `Analyze recent medical breakthroughs and identify the top 5 most promising research areas for drug discovery. For each trend, provide:
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

Focus on areas like: novel therapeutic targets, emerging drug classes, breakthrough mechanisms, underserved conditions, and patent-free opportunities.`
          }
        ],
        temperature: 0.3,
        max_tokens: 2000,
      }),
    });

    if (!response.ok) {
      throw new Error(`Perplexity API error: ${response.statusText}`);
    }

    const result = await response.json();
    const content = result.choices?.[0]?.message?.content;

    if (!content) {
      throw new Error("No content in Perplexity response");
    }

    // Extract JSON from response (handle markdown code blocks)
    let trendsData: MedicalTrend[];
    try {
      // Try to extract JSON from markdown code block
      const jsonMatch = content.match(/```json\n([\s\S]*?)\n```/) || content.match(/```\n([\s\S]*?)\n```/);
      const jsonStr = jsonMatch ? jsonMatch[1] : content;
      trendsData = JSON.parse(jsonStr);
    } catch (parseError) {
      console.error("[Medical Trends] Failed to parse JSON:", content);
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
        source: "Perplexity Sonar Pro",
      }));

    const data: MedicalTrendsData = {
      trends: validatedTrends,
      lastUpdated: new Date().toISOString(),
    };

    await writeFileAsync(TRENDS_FILE, JSON.stringify(data, null, 2));
    console.log(`[Medical Trends] Saved ${validatedTrends.length} trends`);
  } catch (error) {
    console.error("[Medical Trends] Error refreshing trends:", error);
    throw error;
  }
}

/**
 * Alias for refreshMedicalTrends for backward compatibility
 */
export const analyzeMedicalTrends = refreshMedicalTrends;
