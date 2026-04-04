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
      // Find the correct Python executable
      let pythonExe = "python3";
      try {
        // Try to get the absolute path to python3
        pythonExe = execSync("which python3", { encoding: "utf-8" }).trim();
        if (!pythonExe || !existsSync(pythonExe)) {
          pythonExe = "/usr/bin/python3";
        }
      } catch (e) {
        // Fallback to common locations
        const commonPaths = ["/usr/bin/python3", "/usr/local/bin/python3", "/opt/python/bin/python3"];
        for (const path of commonPaths) {
          if (existsSync(path)) {
            pythonExe = path;
            break;
          }
        }
      }

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

      // Spawn Python process with absolute path
      const python = spawn(pythonExe, args, {
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

      python.on("error", (err) => {
        console.error(`[Docking] Process error:`, err);
        return resolve({
          success: false,
          error: `Failed to start docking process: ${err.message}`,
        });
      });
    } catch (error: any) {
      console.error(`[Docking] Wrapper error:`, error);
      return resolve({
        success: false,
        error: `Docking wrapper error: ${error.message}`,
      });
    }
  });
}

/**
 * Get default docking parameters for a target
 * Falls back to sensible defaults if not found in DB
 */
export async function getDockingParamsForTarget(
  targetName: string
): Promise<{
  boxCenterX: number;
  boxCenterY: number;
  boxCenterZ: number;
  boxSizeX: number;
  boxSizeY: number;
  boxSizeZ: number;
  exhaustiveness: number;
  numPoses: number;
}> {
  try {
    const db = await getDb();
    if (!db) {
      console.warn("[Docking] Database not available, using defaults");
      return getDefaultDockingParams();
    }

    const params = await db
      .select()
      .from(dockingParameters)
      .where(eq(dockingParameters.targetName, targetName))
      .limit(1);

    if (params && params.length > 0) {
      return {
        boxCenterX: Number(params[0].boxCenterX),
        boxCenterY: Number(params[0].boxCenterY),
        boxCenterZ: Number(params[0].boxCenterZ),
        boxSizeX: Number(params[0].boxSizeX),
        boxSizeY: Number(params[0].boxSizeY),
        boxSizeZ: Number(params[0].boxSizeZ),
        exhaustiveness: Number(params[0].exhaustiveness),
        numPoses: Number(params[0].numPoses),
      };
    }

    console.warn(
      `[Docking] No parameters found for target ${targetName}, using defaults`
    );
    return getDefaultDockingParams();
  } catch (error: any) {
    console.error(`[Docking] Error fetching parameters:`, error);
    return getDefaultDockingParams();
  }
}

/**
 * Get sensible default docking parameters
 */
function getDefaultDockingParams() {
  return {
    boxCenterX: 0,
    boxCenterY: 0,
    boxCenterZ: 0,
    boxSizeX: 25,
    boxSizeY: 25,
    boxSizeZ: 25,
    exhaustiveness: 8,
    numPoses: 9,
  };
}
