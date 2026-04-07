import { spawn } from "child_process";
import { execSync } from "child_process";
import { existsSync } from "fs";
import { join } from "path";
import { fileURLToPath } from "url";
import { dirname } from "path";
import { getPDBFilePath } from "./pdbStorage";
import { getDb } from "./db";
import { dockingParameters } from "../drizzle/schema";
import { eq } from "drizzle-orm";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Detect Python executable at module load time
let PYTHON_EXECUTABLE = "/usr/bin/python3";
try {
  if (existsSync(PYTHON_EXECUTABLE)) {
    console.log("[Docking] Python found at:", PYTHON_EXECUTABLE);
  } else {
    const detected = execSync("which python3", { encoding: "utf-8" }).trim();
    if (detected && existsSync(detected)) {
      PYTHON_EXECUTABLE = detected;
      console.log("[Docking] Python detected at:", PYTHON_EXECUTABLE);
    }
  }
} catch (e) {
  console.log("[Docking] Python detection failed, using default:", PYTHON_EXECUTABLE);
}

export interface DockingResult {
  success: boolean;
  binding_affinity?: number;
  poses?: Array<{
    mode: number;
    affinity: number;
  }>;
  num_poses?: number;
  error?: string;
}

interface DockingParams {
  smiles: string;
  analogId: string;
  targetName: string;
  boxCenter: { x: number; y: number; z: number };
  boxSize: { x: number; y: number; z: number };
  exhaustiveness: number;
  numPoses: number;
}

/**
 * Run molecular docking using AutoDock Vina
 * Calls scripts/dock.py with SMILES, receptor PDB, and docking parameters
 */
export async function runMolecularDocking(
  params: DockingParams
): Promise<DockingResult> {
  return new Promise((resolve, reject) => {
    try {
      console.log("[Docking] Starting docking with Python:", PYTHON_EXECUTABLE);

      // Get the PDB file path for the target
      const receptorPath = getPDBFilePath(params.targetName);

      // Build command arguments for dock.py
      const args = [
        join(__dirname, "../scripts/dock.py"),
        "--smiles",
        params.smiles,
        "--receptor",
        receptorPath,
        "--cx",
        params.boxCenter.x.toString(),
        "--cy",
        params.boxCenter.y.toString(),
        "--cz",
        params.boxCenter.z.toString(),
        "--sx",
        params.boxSize.x.toString(),
        "--sy",
        params.boxSize.y.toString(),
        "--sz",
        params.boxSize.z.toString(),
        "--exhaustiveness",
        params.exhaustiveness.toString(),
        "--num_poses",
        params.numPoses.toString(),
      ];

      console.log("[Docking] Command:", PYTHON_EXECUTABLE, args.join(" "));

      // Spawn Python process with absolute path
      const python = spawn(PYTHON_EXECUTABLE, args, {
        cwd: __dirname,
        timeout: 300000, // 5 minutes
        shell: false,
        env: { ...process.env, PYTHONUNBUFFERED: "1" },
      });

      let stdout = "";
      let stderr = "";

      python.stdout.on("data", (data) => {
        stdout += data.toString();
      });

      python.stderr.on("data", (data) => {
        stderr += data.toString();
        console.error(`[Docking stderr] ${data}`);
      });

      python.on("close", (code) => {
        try {
          if (code !== 0) {
            console.error(
              `[Docking] Python process exited with code ${code}: ${stderr}`
            );
            return resolve({
              success: false,
              error: `Docking failed: ${stderr || "Unknown error"}`,
            });
          }

          // Parse JSON output from dock.py
          const result = JSON.parse(stdout);
          console.log(`[Docking] Result:`, result);

          return resolve(result);
        } catch (parseError: any) {
          console.error(
            `[Docking] Failed to parse output: ${parseError.message}`
          );
          console.error(`[Docking] stdout: ${stdout}`);
          console.error(`[Docking] stderr: ${stderr}`);

          return resolve({
            success: false,
            error: `Failed to parse docking results: ${parseError.message}`,
          });
        }
      });

      python.on("error", (err: Error) => {
        console.error(`[Docking] Spawn error:`, err);
        return resolve({
          success: false,
          error: `Failed to start docking process: ${err.message}`,
        });
      });
    } catch (error: any) {
      console.error("[Docking] Unexpected error:", error);
      return resolve({
        success: false,
        error: `Docking error: ${error.message}`,
      });
    }
  });
}

/**
 * Save docking results to database
 */
export async function saveDockingResults(
  analogId: number,
  results: DockingResult
): Promise<void> {
  // Note: dockingParameters table doesn't have analogId or results columns
  // This function is kept for backward compatibility but doesn't persist results
  // Results should be saved to batchDockingResults table instead
  console.log(`Docking results for analog ${analogId}:`, results);
}
