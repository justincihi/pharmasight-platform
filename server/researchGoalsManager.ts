import { readFile, writeFile, existsSync, mkdirSync } from "fs";
import { promisify } from "util";

const readFileAsync = promisify(readFile);
const writeFileAsync = promisify(writeFile);

// Use /tmp for production compatibility (Cloud Run has ephemeral /tmp)
// Fall back to local data dir in dev
const DATA_DIR = process.env.DATA_DIR || "/tmp/pharmasight-data";
const GOALS_FILE = `${DATA_DIR}/research_goals.json`;

// Also check legacy path for backward compatibility
const LEGACY_GOALS_FILE = "/home/ubuntu/pharmasight-admin-dashboard/data/research_goals.json";

function ensureDataDir() {
  try {
    if (!existsSync(DATA_DIR)) {
      mkdirSync(DATA_DIR, { recursive: true });
    }
  } catch {
    // Ignore
  }
}

interface ResearchGoals {
  goals: string[];
  lastUpdated: string;
}

const DEFAULT_GOALS: ResearchGoals = {
  goals: ["psychedelics", "nootropics", "anxiolytics", "neuroprotection", "kava-analogs"],
  lastUpdated: new Date().toISOString(),
};

/**
 * Get saved research goals
 */
export async function getResearchGoals(): Promise<ResearchGoals> {
  ensureDataDir();
  
  // Try primary location
  for (const file of [GOALS_FILE, LEGACY_GOALS_FILE]) {
    if (existsSync(file)) {
      try {
        const data = await readFileAsync(file, "utf-8");
        return JSON.parse(data);
      } catch {
        continue;
      }
    }
  }

  // Return defaults if no file found
  return DEFAULT_GOALS;
}

/**
 * Save research goals
 */
export async function saveResearchGoals(goals: string[]): Promise<void> {
  ensureDataDir();
  
  const data: ResearchGoals = {
    goals,
    lastUpdated: new Date().toISOString(),
  };

  try {
    await writeFileAsync(GOALS_FILE, JSON.stringify(data, null, 2));
    console.log(`[Research Goals] Saved ${goals.length} goals to ${GOALS_FILE}`);
  } catch (error: any) {
    // Try legacy path as fallback
    try {
      await writeFileAsync(LEGACY_GOALS_FILE, JSON.stringify(data, null, 2));
      console.log(`[Research Goals] Saved ${goals.length} goals to legacy path`);
    } catch (e2) {
      console.error("[Research Goals] Error saving goals:", error);
      throw error;
    }
  }
}

/**
 * Alias for getResearchGoals for backward compatibility
 */
export const loadResearchGoals = getResearchGoals;
