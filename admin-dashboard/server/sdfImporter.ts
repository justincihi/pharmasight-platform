import { spawn } from "child_process";
import { join, dirname } from "path";
import { fileURLToPath } from "url";
import { getDb } from "./db";
import { analogDiscoveries } from "../drizzle/schema";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const PYTHON_MODULES_PATH = join(__dirname, "python_modules");

interface SDFImportResult {
  success: boolean;
  imported: number;
  analogs: any[];
  errors: string[];
}

/**
 * Import analogs from SDF file using Python processor
 */
export async function importFromSDF(sdfPath: string): Promise<SDFImportResult> {
  return new Promise((resolve) => {
    const pythonScript = `
import sys
import json
sys.path.insert(0, '${PYTHON_MODULES_PATH}')

try:
    from nmda_admet_analyzer import analyze_sdf_file
    
    # Process SDF file with ADMET analysis
    result_json = analyze_sdf_file('${sdfPath}')
    analogs = json.loads(result_json)
    
    # Return results
    print(json.dumps({
        'success': True,
        'analogs': analogs
    }))
    
except Exception as e:
    print(json.dumps({
        'success': False,
        'error': str(e)
    }))
    sys.exit(1)
`;

    const python = spawn("python3", ["-c", pythonScript]);
    let stdout = "";
    let stderr = "";

    python.stdout.on("data", (data) => {
      stdout += data.toString();
    });

    python.stderr.on("data", (data) => {
      stderr += data.toString();
    });

    python.on("close", async (code) => {
      if (code !== 0) {
        console.error("[SDF Import] Python error:", stderr);
        resolve({
          success: false,
          imported: 0,
          analogs: [],
          errors: [stderr || "Python script failed"],
        });
        return;
      }

      try {
        // Parse Python output
        const lines = stdout.trim().split("\n");
        const lastLine = lines[lines.length - 1];
        const result = JSON.parse(lastLine);

        if (!result.success) {
          resolve({
            success: false,
            imported: 0,
            analogs: [],
            errors: [result.error],
          });
          return;
        }

        // Import analogs into database
        const db = await getDb();
        if (!db) {
          resolve({
            success: false,
            imported: 0,
            analogs: [],
            errors: ["Database connection failed"],
          });
          return;
        }
        
        const imported: any[] = [];
        const errors: string[] = [];

        for (const analog of result.analogs) {
          try {
            // Generate unique compound ID
            const timestamp = Date.now();
            const randomSuffix = Math.random().toString(36).substring(2, 8).toUpperCase();
            const compoundId = `${analog.compound_name}-${timestamp}-${randomSuffix}`;
            
            // Map analog data to database schema
            const analogData = {
              compoundId,
              compoundName: analog.compound_name,
              parentCompound: analog.mechanism || "NMDA Antagonist",
              smiles: analog.smiles,
              
              // Required scores
              confidenceScore: Math.round(analog.overall_admet_score || 85),
              similarityScore: Math.round(analog.confidence_score || 85),
              safetyScore: Math.round(100 - (analog.hepatotoxicity_risk_score || 20)),
              efficacyScore: Math.round(analog.nmda_binding_potential || 85),
              drugLikenessScore: Math.round(analog.drug_likeness_score || 100),
              
              // Patent and regulatory
              patentStatus: (analog.patent_status || "patent-opportunity") as any,
              marketValue: analog.market_value ? `$${analog.market_value}` : null,
              
              // Therapeutic info
              therapeuticPotential: analog.admet_explanation || "",
              keyDifferences: analog.description || "",
              mechanismOfAction: analog.mechanism || "NMDA receptor antagonist",
              
              // Molecular descriptors
              molecularWeight: analog.molecular_weight ? analog.molecular_weight.toString() : null,
              logP: analog.clogp ? analog.clogp.toString() : null,
              hBondDonors: analog.hbd || 0,
              hBondAcceptors: analog.hba || 0,
              
              // Discovery metadata
              discoveredBy: "sdf-import",
              discoveryMethod: "manual-sdf-upload",
              discoveredAt: new Date(),
              
              // Approval
              approvalStatus: "pending" as any,
              
              createdAt: new Date(),
              updatedAt: new Date(),
            };

            const [inserted] = await db
              .insert(analogDiscoveries)
              .values(analogData)
              .$returningId();

            imported.push({ ...analogData, id: inserted.id });
          } catch (error: any) {
            errors.push(`Failed to import ${analog.compound_name}: ${error.message}`);
          }
        }

        resolve({
          success: true,
          imported: imported.length,
          analogs: imported,
          errors,
        });
      } catch (error: any) {
        console.error("[SDF Import] Parse error:", error);
        resolve({
          success: false,
          imported: 0,
          analogs: [],
          errors: [error.message],
        });
      }
    });
  });
}
