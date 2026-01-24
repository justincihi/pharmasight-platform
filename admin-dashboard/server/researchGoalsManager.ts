import { readFile, writeFile, existsSync } from "fs";
import { promisify } from "util";

const readFileAsync = promisify(readFile);
const writeFileAsync = promisify(writeFile);

const GOALS_FILE = "/home/ubuntu/pharmasight-admin-dashboard/data/research_goals.json";

interface ResearchGoals {
  goals: string[];
  lastUpdated: string;
}

/**
 * Get saved research goals
 */
export async function getResearchGoals(): Promise<ResearchGoals> {
  try {
    if (!existsSync(GOALS_FILE)) {
      // Return default goals if file doesn't exist
      return {
        goals: ['psychedelics', 'nootropics', 'anxiolytics'],
        lastUpdated: new Date().toISOString(),
      };
    }

    const data = await readFileAsync(GOALS_FILE, "utf-8");
    return JSON.parse(data);
  } catch (error) {
    console.error("[Research Goals] Error reading goals:", error);
    return {
      goals: ['psychedelics', 'nootropics', 'anxiolytics'],
      lastUpdated: new Date().toISOString(),
    };
  }
}

/**
 * Save research goals
 */
export async function saveResearchGoals(goals: string[]): Promise<void> {
  try {
    const data: ResearchGoals = {
      goals,
      lastUpdated: new Date().toISOString(),
    };

    await writeFileAsync(GOALS_FILE, JSON.stringify(data, null, 2));
    console.log(`[Research Goals] Saved ${goals.length} goals`);
  } catch (error) {
    console.error("[Research Goals] Error saving goals:", error);
    throw error;
  }
}

/**
 * Alias for getResearchGoals for backward compatibility
 */
export const loadResearchGoals = getResearchGoals;
