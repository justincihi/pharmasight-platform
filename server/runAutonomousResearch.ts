import { spawn } from "child_process";
import { join, dirname } from "path";
import { fileURLToPath } from "url";
import { writeFileSync, readFileSync, existsSync } from "fs";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const PYTHON_MODULES_PATH = "/home/ubuntu/pharmasight-admin-dashboard/server/python_modules";
const OUTPUT_FILE = "/home/ubuntu/pharmasight-admin-dashboard/data/autonomous_research_output.json";

interface ResearchResult {
  success: boolean;
  discoveries?: any[];
  error?: string;
}

/**
 * Run the autonomous research engine to discover new analogs
 * This calls the daily_discovery_engine.py module
 */
export async function runAutonomousResearch(): Promise<ResearchResult> {
  return new Promise((resolve) => {
    const pythonScript = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    # Import the autonomous research engine
    from daily_discovery_engine import DailyDiscoveryEngine
    
    # Initialize and run discovery
    engine = DailyDiscoveryEngine()
    report = engine.generate_daily_report(
        goals=['psychedelics', 'nootropics', 'anxiolytics']
    )
    
    # Extract discoveries from report
    results = report.get('discoveries', [])
    
    # Output results as JSON
    print(json.dumps({
        "success": True,
        "discoveries": results
    }))
except Exception as e:
    print(json.dumps({
        "success": False,
        "error": str(e)
    }))
`;

    // Use venv Python if available to avoid SRE module mismatch
    const venvPython = "/home/ubuntu/pharmasight-admin-dashboard/server/python_modules/venv/bin/python3";
    const pythonCmd = existsSync(venvPython) ? venvPython : "python3";
    
    const python = spawn(pythonCmd, ["-c", pythonScript]);

    let stdout = "";
    let stderr = "";

    python.stdout.on("data", (data) => {
      stdout += data.toString();
    });

    python.stderr.on("data", (data) => {
      stderr += data.toString();
    });

    python.on("close", (code) => {
      if (code !== 0) {
        console.error("[Autonomous Research] Python error:", stderr);
        resolve({
          success: false,
          error: stderr || "Python process exited with error",
        });
        return;
      }

      try {
        const result = JSON.parse(stdout.trim());
        
        // Save results to file for import
        if (result.success && result.discoveries) {
          writeFileSync(OUTPUT_FILE, JSON.stringify(result.discoveries, null, 2));
        }
        
        resolve(result);
      } catch (error) {
        console.error("[Autonomous Research] Failed to parse output:", error);
        resolve({
          success: false,
          error: "Failed to parse research results",
        });
      }
    });
  });
}

/**
 * Get the latest research results from file
 */
export function getLatestResearchResults(): any[] {
  if (!existsSync(OUTPUT_FILE)) {
    return [];
  }

  try {
    const data = readFileSync(OUTPUT_FILE, "utf-8");
    return JSON.parse(data);
  } catch (error) {
    console.error("[Autonomous Research] Failed to read results:", error);
    return [];
  }
}
