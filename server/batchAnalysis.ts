import { z } from "zod";
import { executePythonScriptSafe } from "./_core/pythonBridgeSafe";

export interface BatchTestResult {
  analogId: number;
  compoundName: string;
  status: "pending" | "running" | "completed" | "failed";
  results?: {
    admet?: any;
    docking?: any;
    toxicity?: any;
    pkpd?: any;
  };
  error?: string;
}

export async function runBatchAnalysis(
  analogIds: number[],
  tests: string[],
  getAnalogById: (id: number) => Promise<any>
): Promise<{ results: BatchTestResult[]; completed: number; failed: number; total: number }> {
  const results: BatchTestResult[] = [];
  let completed = 0;
  let failed = 0;

  for (const analogId of analogIds) {
    const analog = await getAnalogById(analogId);
    if (!analog) {
      results.push({
        analogId,
        compoundName: `Unknown-${analogId}`,
        status: "failed",
        error: "Analog not found",
      });
      failed++;
      continue;
    }

    const testResults: any = {};
    let allTestsSucceeded = true;

    for (const test of tests) {
      try {
        let result;
        switch (test) {
          case "admet":
            result = await executePythonScriptSafe("admet_predictor_advanced.py", "predict_admet", [analog.smiles]);
            testResults.admet = result;
            break;
          case "docking":
            result = await executePythonScriptSafe("molecular_docking.py", "perform_docking", [analog.smiles, "default"]);
            testResults.docking = result;
            break;
          case "toxicity":
            result = await executePythonScriptSafe("toxicity_prediction.py", "predict_toxicity", [analog.smiles]);
            testResults.toxicity = result;
            break;
          case "pkpd":
            result = await executePythonScriptSafe("pkpd_pbpk_simulator.py", "simulate_pkpd", [analog.smiles, 100]);
            testResults.pkpd = result;
            break;
        }
      } catch (error) {
        console.error(`Test ${test} failed for analog ${analogId}:`, error);
        allTestsSucceeded = false;
      }
    }

    if (allTestsSucceeded) {
      results.push({
        analogId,
        compoundName: analog.compoundName,
        status: "completed",
        results: testResults,
      });
      completed++;
    } else {
      results.push({
        analogId,
        compoundName: analog.compoundName,
        status: "failed",
        results: testResults,
        error: "Some tests failed",
      });
      failed++;
    }
  }

  return {
    results,
    completed,
    failed,
    total: analogIds.length,
  };
}
