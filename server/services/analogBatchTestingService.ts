/**
 * Analog Batch Testing Service
 * Orchestrates ADMET, docking, and PK/PD analyses for multiple analogs
 */

export interface AnalogTestingJob {
  id: string;
  analogIds: number[];
  status: "pending" | "running" | "completed" | "failed";
  createdAt: Date;
  startedAt?: Date;
  completedAt?: Date;
  results: AnalogTestResult[];
  error?: string;
}

export interface AnalogTestResult {
  analogId: number;
  smiles: string;
  compoundName: string;
  admetResult?: {
    herg: number;
    hepatotoxicity: number;
    mutagenicity: number;
    carcinogenicity: number;
    confidence: number;
  };
  dockingResult?: {
    bindingAffinity: number;
    rmsd: number;
    poses: number;
  };
  pkpdResult?: {
    tmax: number;
    vd: number;
    halfLife: number;
    clearance: number;
    ec50: number;
    emax: number;
  };
  status: "pending" | "running" | "completed" | "failed";
  error?: string;
}

/**
 * Create a new batch testing job
 */
export async function createBatchTestingJob(
  analogIds: number[]
): Promise<AnalogTestingJob> {
  const jobId = `batch-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;

  const job: AnalogTestingJob = {
    id: jobId,
    analogIds,
    status: "pending",
    createdAt: new Date(),
    results: analogIds.map((id) => ({
      analogId: id,
      smiles: "",
      compoundName: "",
      status: "pending",
    })),
  };

  console.log(`Created batch testing job: ${jobId} with ${analogIds.length} analogs`);
  return job;
}

/**
 * Run ADMET analysis on an analog
 */
export async function runADMETAnalysis(smiles: string): Promise<any> {
  try {
    const { runPythonJson } = await import("../pythonBridge.js");

    const result = await runPythonJson({
      modulePath: "comprehensive_analysis.py",
      args: ["toxicity", smiles],
      timeoutMs: 30000,
    });

    if (!result.ok) {
      throw new Error(result.error?.message || "ADMET analysis failed");
    }

    return result.data;
  } catch (error) {
    console.error("ADMET analysis error:", error);
    throw error;
  }
}

/**
 * Run molecular docking on an analog
 */
export async function runDockingAnalysis(
  smiles: string,
  receptorPath: string = "/home/ubuntu/pharmasight-admin-dashboard/receptors/NMDA_6NQH.pdbqt"
): Promise<any> {
  try {
    const { runPythonJson } = await import("../pythonBridge.js");

    const result = await runPythonJson({
      modulePath: "../scripts/dock.py",
      args: [
        "--smiles",
        smiles,
        "--receptor",
        receptorPath,
        "--cx",
        "0",
        "--cy",
        "0",
        "--cz",
        "0",
        "--sx",
        "20",
        "--sy",
        "20",
        "--sz",
        "20",
      ],
      timeoutMs: 60000,
    });

    if (!result.ok) {
      throw new Error(result.error?.message || "Docking analysis failed");
    }

    return result.data;
  } catch (error) {
    console.error("Docking analysis error:", error);
    throw error;
  }
}

/**
 * Run PK/PD simulation on an analog
 */
export async function runPKPDSimulation(smiles: string): Promise<any> {
  try {
    const { runPythonJson } = await import("../pythonBridge.js");

    const result = await runPythonJson({
      modulePath: "pkpd_pbpk_simulator.py",
      args: ["simulate", smiles],
      timeoutMs: 45000,
    });

    if (!result.ok) {
      throw new Error(result.error?.message || "PK/PD simulation failed");
    }

    return result.data;
  } catch (error) {
    console.error("PK/PD simulation error:", error);
    throw error;
  }
}

/**
 * Run all tests on a single analog
 */
export async function runAllTestsOnAnalog(
  analogId: number,
  smiles: string,
  compoundName: string
): Promise<AnalogTestResult> {
  const result: AnalogTestResult = {
    analogId,
    smiles,
    compoundName,
    status: "running",
  };

  try {
    // Run ADMET
    console.log(`Running ADMET for ${compoundName}...`);
    const admetData = await runADMETAnalysis(smiles);
    result.admetResult = admetData;

    // Run docking
    console.log(`Running docking for ${compoundName}...`);
    const dockingData = await runDockingAnalysis(smiles);
    result.dockingResult = dockingData;

    // Run PK/PD
    console.log(`Running PK/PD for ${compoundName}...`);
    const pkpdData = await runPKPDSimulation(smiles);
    result.pkpdResult = pkpdData;

    result.status = "completed";
    console.log(`Completed all tests for ${compoundName}`);
  } catch (error) {
    result.status = "failed";
    result.error = error instanceof Error ? error.message : "Unknown error";
    console.error(`Failed to test ${compoundName}:`, error);
  }

  return result;
}

/**
 * Process a batch testing job
 */
export async function processBatchTestingJob(
  job: AnalogTestingJob,
  analogData: Array<{ id: number; smiles: string; compoundName: string }>
): Promise<AnalogTestingJob> {
  job.status = "running";
  job.startedAt = new Date();

  try {
    // Process each analog sequentially to avoid overwhelming system
    for (let i = 0; i < job.results.length; i++) {
      const result = job.results[i];
      const analog = analogData.find((a) => a.id === result.analogId);

      if (!analog) {
        result.status = "failed";
        result.error = "Analog data not found";
        continue;
      }

      // Update with analog data
      result.smiles = analog.smiles;
      result.compoundName = analog.compoundName;

      // Run all tests
      const testResult = await runAllTestsOnAnalog(
        analog.id,
        analog.smiles,
        analog.compoundName
      );

      // Merge test results
      job.results[i] = {
        ...result,
        ...testResult,
      };

      // Log progress
      console.log(`Progress: ${i + 1}/${job.results.length} analogs tested`);
    }

    job.status = "completed";
    job.completedAt = new Date();
  } catch (error) {
    job.status = "failed";
    job.error = error instanceof Error ? error.message : "Unknown error";
    job.completedAt = new Date();
  }

  return job;
}

/**
 * Generate summary statistics for batch results
 */
export function generateBatchSummary(job: AnalogTestingJob): {
  totalAnalogs: number;
  completedTests: number;
  failedTests: number;
  avgBindingAffinity: number;
  avgQED: number;
  topCandidates: AnalogTestResult[];
} {
  const completedResults = job.results.filter((r) => r.status === "completed");
  const failedResults = job.results.filter((r) => r.status === "failed");

  const bindingAffinities = completedResults
    .map((r) => r.dockingResult?.bindingAffinity || 0)
    .filter((v) => v !== 0);

  const avgBindingAffinity =
    bindingAffinities.length > 0
      ? bindingAffinities.reduce((a, b) => a + b, 0) / bindingAffinities.length
      : 0;

  // Sort by binding affinity (lower is better)
  const topCandidates = [...completedResults]
    .sort(
      (a, b) =>
        (a.dockingResult?.bindingAffinity || 999) -
        (b.dockingResult?.bindingAffinity || 999)
    )
    .slice(0, 5);

  return {
    totalAnalogs: job.results.length,
    completedTests: completedResults.length,
    failedTests: failedResults.length,
    avgBindingAffinity,
    avgQED: 0.7, // Placeholder
    topCandidates,
  };
}
