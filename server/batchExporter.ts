import { spawn } from "child_process";
import { fileURLToPath } from "url";
import { join, dirname } from "path";
import { getDb } from "./db";
import { analogDiscoveries } from "../drizzle/schema";
import { sql } from "drizzle-orm";
import { writeFile } from "fs/promises";
import { tmpdir } from "os";
import { randomBytes } from "crypto";

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const PYTHON_MODULES_PATH = join(__dirname, "python_modules");

interface ExportOptions {
  format: "csv" | "sdf";
  filters?: {
    minConfidence?: number;
    patentStatus?: string[];
    therapeuticAreas?: string[];
    dateRange?: { start: Date; end: Date };
  };
}

interface ExportResult {
  success: boolean;
  filePath?: string;
  error?: string;
  count: number;
}

/**
 * Export analogs to CSV format
 */
async function exportToCSV(analogs: any[]): Promise<string> {
  const headers = [
    "Compound ID",
    "Compound Name",
    "SMILES",
    "Parent Compound",
    "Confidence Score",
    "Safety Score",
    "Efficacy Score",
    "Drug Likeness",
    "Patent Status",
    "Therapeutic Potential",
    "Mechanism of Action",
    "Molecular Weight",
    "LogP",
    "H-Bond Donors",
    "H-Bond Acceptors",
    "Discovered At",
    "Approval Status",
  ];

  const rows = analogs.map((analog) => [
    analog.compoundId,
    analog.compoundName,
    analog.smiles,
    analog.parentCompound,
    analog.confidenceScore,
    analog.safetyScore,
    analog.efficacyScore,
    analog.drugLikenessScore,
    analog.patentStatus,
    `"${(analog.therapeuticPotential || "").replace(/"/g, '""')}"`,
    analog.mechanismOfAction || "",
    analog.molecularWeight || "",
    analog.logP || "",
    analog.hBondDonors || 0,
    analog.hBondAcceptors || 0,
    new Date(analog.discoveredAt).toISOString(),
    analog.approvalStatus,
  ]);

  const csvContent = [
    headers.join(","),
    ...rows.map((row) => row.join(",")),
  ].join("\n");

  // Write to temporary file
  const tmpFile = join(tmpdir(), `pharmasight-export-${randomBytes(8).toString("hex")}.csv`);
  await writeFile(tmpFile, csvContent, "utf-8");

  return tmpFile;
}

/**
 * Export analogs to SDF format using Python RDKit
 */
async function exportToSDF(analogs: any[]): Promise<string> {
  return new Promise((resolve, reject) => {
    const tmpFile = join(tmpdir(), `pharmasight-export-${randomBytes(8).toString("hex")}.sdf`);

    // Prepare analog data for Python script
    const analogsJSON = JSON.stringify(analogs);
    const scriptPath = join(PYTHON_MODULES_PATH, 'export_to_sdf.py');

    const python = spawn("python3", [scriptPath, analogsJSON, tmpFile]);
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
        console.error("[SDF Export] Python error:", stderr);
        reject(new Error(stderr || "SDF export failed"));
        return;
      }

      try {
        const lines = stdout.trim().split("\n");
        const lastLine = lines[lines.length - 1];
        const result = JSON.parse(lastLine);

        if (!result.success) {
          reject(new Error(result.error));
          return;
        }

        resolve(result.filePath);
      } catch (error: any) {
        reject(new Error(`Failed to parse Python output: ${error.message}`));
      }
    });
  });
}

/**
 * Export analogs with filtering options
 */
export async function exportAnalogs(options: ExportOptions): Promise<ExportResult> {
  try {
    const db = await getDb();
    if (!db) {
      return {
        success: false,
        error: "Database connection failed",
        count: 0,
      };
    }

    // Build query with filters
    let queryBuilder = db.select().from(analogDiscoveries);

    const conditions: any[] = [];

    if (options.filters?.minConfidence) {
      conditions.push(sql`${analogDiscoveries.confidenceScore} >= ${options.filters.minConfidence}`);
    }

    if (options.filters?.patentStatus && options.filters.patentStatus.length > 0) {
      conditions.push(sql`${analogDiscoveries.patentStatus} IN (${sql.join(options.filters.patentStatus.map(s => sql`${s}`), sql`, `)})`);
    }

    if (options.filters?.dateRange) {
      conditions.push(sql`${analogDiscoveries.discoveredAt} >= ${options.filters.dateRange.start}`);
      conditions.push(sql`${analogDiscoveries.discoveredAt} <= ${options.filters.dateRange.end}`);
    }

    if (conditions.length > 0) {
      queryBuilder = queryBuilder.where(sql`${sql.join(conditions, sql` AND `)}`) as any;
    }

    const analogs = await queryBuilder;

    // Filter by therapeutic area (text search since it's not a separate column)
    let filteredAnalogs = analogs;
    if (options.filters?.therapeuticAreas && options.filters.therapeuticAreas.length > 0) {
      filteredAnalogs = analogs.filter((analog: any) => {
        const potential = (analog.therapeuticPotential || "").toLowerCase();
        return options.filters!.therapeuticAreas!.some((area) =>
          potential.includes(area.toLowerCase())
        );
      });
    }

    if (filteredAnalogs.length === 0) {
      return {
        success: false,
        error: "No analogs match the specified filters",
        count: 0,
      };
    }

    // Export based on format
    let filePath: string;
    if (options.format === "csv") {
      filePath = await exportToCSV(filteredAnalogs);
    } else {
      filePath = await exportToSDF(filteredAnalogs);
    }

    return {
      success: true,
      filePath,
      count: filteredAnalogs.length,
    };
  } catch (error: any) {
    console.error("[Batch Export] Error:", error);
    return {
      success: false,
      error: error.message,
      count: 0,
    };
  }
}
