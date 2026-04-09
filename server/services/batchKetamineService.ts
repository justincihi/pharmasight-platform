import { getDb } from "../db";
import { analogDiscoveries } from "../../drizzle/schema";
import { like } from "drizzle-orm";
import { runMolecularDocking, DockingResult } from "../molecularDockingWrapper";

export interface KetamineBatchConfig {
  compoundFilter?: string;
  receptors: string[];
  dockingParams?: {
    exhaustiveness?: number;
    numPoses?: number;
    centerX?: number;
    centerY?: number;
    centerZ?: number;
    sizeX?: number;
    sizeY?: number;
    sizeZ?: number;
  };
}

export interface BatchKetamineResult {
  analogId: string;
  compoundName: string;
  smiles: string;
  receptor: string;
  bindingAffinity: number;
  rmsd: number;
  numPoses: number;
  timestamp: Date;
}

/**
 * Fetch all ketamine analogs from database
 */
export async function getKetamineAnalogs(filter?: string) {
  try {
    const db = await getDb();
    if (!db) {
      throw new Error("Database not available");
    }

    let query = db.select().from(analogDiscoveries);
    
    if (filter) {
      query = query.where(like(analogDiscoveries.compoundName, `%${filter}%`)) as any;
    } else {
      // Default: fetch all compounds with "KETAMINE" in name
      query = query.where(like(analogDiscoveries.compoundName, "%KETAMINE%")) as any;
    }

    const results = await (query as any).execute();
    return results;
  } catch (error) {
    console.error("[Batch Ketamine] Error fetching analogs:", error);
    throw error;
  }
}

/**
 * Run batch docking for ketamine analogs
 */
export async function runBatchKetamineDocking(config: KetamineBatchConfig) {
  try {
    // Fetch ketamine analogs
    const ketamineAnalogs = await getKetamineAnalogs(config.compoundFilter);
    
    if (ketamineAnalogs.length === 0) {
      throw new Error("No ketamine analogs found");
    }

    console.log(`[Batch Ketamine] Found ${ketamineAnalogs.length} analogs`);
    console.log(`[Batch Ketamine] Docking against ${config.receptors.length} receptors`);

    const results: BatchKetamineResult[] = [];
    const errors: { analog: string; receptor: string; error: string }[] = [];

    // For each analog
    for (const analog of ketamineAnalogs) {
      // For each receptor
      for (const receptor of config.receptors) {
        try {
          console.log(`[Batch Ketamine] Docking ${analog.compoundName} to ${receptor}...`);

          // Run docking
          const dockingResult = await runMolecularDocking({
            smiles: analog.smiles || "",
            analogId: analog.compoundId,
            targetName: receptor,
            boxCenter: {
              x: config.dockingParams?.centerX ?? 0,
              y: config.dockingParams?.centerY ?? 0,
              z: config.dockingParams?.centerZ ?? 0,
            },
            boxSize: {
              x: config.dockingParams?.sizeX ?? 20,
              y: config.dockingParams?.sizeY ?? 20,
              z: config.dockingParams?.sizeZ ?? 20,
            },
            exhaustiveness: config.dockingParams?.exhaustiveness ?? 8,
            numPoses: config.dockingParams?.numPoses ?? 9,
          });

          if (dockingResult.success && dockingResult.binding_affinity !== undefined && dockingResult.num_poses !== undefined) {
            results.push({
              analogId: analog.compoundId,
              compoundName: analog.compoundName,
              smiles: analog.smiles || "",
              receptor,
              bindingAffinity: dockingResult.binding_affinity,
              rmsd: dockingResult.poses?.[0]?.affinity ?? 0,
              numPoses: dockingResult.num_poses,
              timestamp: new Date(),
            });

            console.log(
              `[Batch Ketamine] ✓ ${analog.compoundName} → ${receptor}: ${dockingResult.binding_affinity.toFixed(2)} kcal/mol`
            );
          } else {
            errors.push({
              analog: analog.compoundName,
              receptor,
              error: dockingResult.error || "Unknown error",
            });

            console.error(
              `[Batch Ketamine] ✗ ${analog.compoundName} → ${receptor}: ${dockingResult.error}`
            );
          }
        } catch (error) {
          const errorMsg = error instanceof Error ? error.message : String(error);
          errors.push({
            analog: analog.compoundName,
            receptor,
            error: errorMsg,
          });

          console.error(
            `[Batch Ketamine] Exception: ${analog.compoundName} → ${receptor}:`,
            error
          );
        }
      }
    }

    return {
      success: results.length > 0,
      resultsCount: results.length,
      errorsCount: errors.length,
      results,
      errors,
    };
  } catch (error) {
    console.error("[Batch Ketamine] Batch docking failed:", error);
    throw error;
  }
}

/**
 * Calculate selectivity scores for a compound across receptors
 */
export function calculateSelectivityScores(
  results: BatchKetamineResult[],
  primaryReceptor: string
) {
  const grouped = results.reduce(
    (acc, result) => {
      if (!acc[result.compoundName]) {
        acc[result.compoundName] = {};
      }
      acc[result.compoundName][result.receptor] = result.bindingAffinity;
      return acc;
    },
    {} as Record<string, Record<string, number>>
  );

  const selectivityScores = Object.entries(grouped).map(([compound, affinities]) => {
    const primaryAffinity = affinities[primaryReceptor] ?? 0;
    const otherAffinities = Object.entries(affinities)
      .filter(([receptor]) => receptor !== primaryReceptor)
      .map(([, affinity]) => affinity);

    // Selectivity = primary affinity - average off-target affinity
    // Higher selectivity score = better (more selective)
    const avgOffTarget = otherAffinities.length > 0
      ? otherAffinities.reduce((a, b) => a + b, 0) / otherAffinities.length
      : 0;

    const selectivityScore = primaryAffinity - avgOffTarget;

    return {
      compound,
      primaryAffinity,
      avgOffTarget,
      selectivityScore,
      affinities,
    };
  });

  return selectivityScores.sort((a, b) => b.selectivityScore - a.selectivityScore);
}

/**
 * Export batch results to CSV format
 */
export function exportBatchResultsToCSV(results: BatchKetamineResult[]): string {
  const headers = [
    "Compound Name",
    "SMILES",
    "Receptor",
    "Binding Affinity (kcal/mol)",
    "RMSD (Å)",
    "Num Poses",
    "Timestamp",
  ];

  const rows = results.map((r) => [
    r.compoundName,
    r.smiles,
    r.receptor,
    r.bindingAffinity.toFixed(2),
    r.rmsd.toFixed(2),
    r.numPoses,
    r.timestamp.toISOString(),
  ]);

  const csvContent = [
    headers.join(","),
    ...rows.map((row) => row.map((cell) => `"${cell}"`).join(",")),
  ].join("\n");

  return csvContent;
}

/**
 * Export batch results to JSON format
 */
export function exportBatchResultsToJSON(results: BatchKetamineResult[]): string {
  return JSON.stringify(results, null, 2);
}
