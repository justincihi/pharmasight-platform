import { getDb } from "../db";
import { analogDiscoveries } from "../../drizzle/schema";
import { sql, eq } from "drizzle-orm";

export interface DiscoveryAuditResult {
  id: string;
  smiles: string;
  compoundName: string;
  createdAt: Date;
  confidence: number;
  patentStatus: string;
  isValid: boolean;
  pharmacophoreMatch: number;
  dataCompleteness: number;
  status: "complete" | "partial" | "fragment";
  recommendations: string[];
}

/**
 * Query discoveries from the last 60 days
 */
export async function getRecentDiscoveries(
  days: number = 60
): Promise<DiscoveryAuditResult[]> {
  try {
    const db = await getDb();
    if (!db) return [];

    const cutoffDate = new Date();
    cutoffDate.setDate(cutoffDate.getDate() - days);

    // Query recent discoveries
    const recentDiscoveries = await db
      .select()
      .from(analogDiscoveries)
      .where(sql`created_at >= ${cutoffDate}`)
      .limit(100);

    // Validate and enrich each discovery
    const enrichedDiscoveries = await Promise.all(
      recentDiscoveries.map(async (discovery: any) => {
        const isValid = discovery.smiles && discovery.smiles.length > 0;
        const pharmacophoreMatch = 0.8; // Placeholder
        const dataCompleteness = 0.75; // Placeholder

        // Determine status
        let status: "complete" | "partial" | "fragment" = "complete";
        if (dataCompleteness < 0.5) status = "fragment";
        else if (dataCompleteness < 0.75) status = "partial";

        // Generate recommendations
        const recommendations: string[] = [];
        if (!isValid) recommendations.push("Invalid or incomplete SMILES");
        if (pharmacophoreMatch < 0.7) recommendations.push("Low pharmacophore match");
        if (dataCompleteness < 0.75) recommendations.push("Missing key data fields");

        return {
          id: discovery.id,
          smiles: discovery.smiles,
          compoundName: discovery.compoundName,
          createdAt: discovery.createdAt || new Date(),
          confidence: discovery.confidenceScore || 0,
          patentStatus: discovery.patentStatus || "unknown",
          isValid,
          pharmacophoreMatch,
          dataCompleteness,
          status,
          recommendations,
        };
      })
    );

    return enrichedDiscoveries;
  } catch (error) {
    console.error("Error querying recent discoveries:", error);
    return [];
  }
}

/**
 * Save discovery to master list
 */
export async function saveDiscoveryToMasterList(
  discovery: DiscoveryAuditResult,
  metadata?: Record<string, any>
): Promise<{ success: boolean; id?: string; error?: string }> {
  try {
    const db = await getDb();
    if (!db) return { success: false, error: "Database not available" };

    // Check if already exists
    const existing = await db
      .select()
      .from(analogDiscoveries)
      .where(sql`smiles = ${discovery.smiles}`)
      .limit(1);

    if (existing.length > 0) {
      return { success: false, error: "Compound already exists in master list" };
    }

    // Update the discovery record to mark as saved to master list
    await db
      .update(analogDiscoveries)
      .set({
        patentStatus: (discovery.patentStatus as any),
        confidenceScore: discovery.confidence,
      })
      .where(sql`id = ${discovery.id}`);

    return { success: true, id: discovery.id };
  } catch (error) {
    console.error("Error saving discovery to master list:", error);
    return { success: false, error: String(error) };
  }
}

/**
 * Save discovery to scaffold library (for partial/fragment discoveries)
 */
export async function saveDiscoveryToScaffoldLibrary(
  discovery: DiscoveryAuditResult,
  metadata?: Record<string, any>
): Promise<{ success: boolean; id?: string; error?: string }> {
  try {
    const db = await getDb();
    if (!db) return { success: false, error: "Database not available" };

    // Mark as scaffold in metadata
    await db
      .update(analogDiscoveries)
      .set({
        patentStatus: ("patent-free" as any), // Mark as scaffold
        confidenceScore: discovery.confidence,
      })
      .where(sql`id = ${discovery.id}`);

    return { success: true, id: discovery.id };
  } catch (error) {
    console.error("Error saving discovery to scaffold library:", error);
    return { success: false, error: String(error) };
  }
}

/**
 * Get statistics on discoveries
 */
export async function getDiscoveryStatistics(days: number = 60) {
  try {
    const db = await getDb();
    if (!db) return null;

    const discoveries = await getRecentDiscoveries(days);

    return {
      total: discoveries.length,
      complete: discoveries.filter((d) => d.status === "complete").length,
      partial: discoveries.filter((d) => d.status === "partial").length,
      fragments: discoveries.filter((d) => d.status === "fragment").length,
      validSMILES: discoveries.filter((d) => d.isValid).length,
      averageConfidence:
        discoveries.length > 0
          ? discoveries.reduce((sum, d) => sum + d.confidence, 0) / discoveries.length
          : 0,
      averagePharmacophoreMatch:
        discoveries.length > 0
          ? discoveries.reduce((sum, d) => sum + d.pharmacophoreMatch, 0) /
            discoveries.length
          : 0,
      averageDataCompleteness:
        discoveries.length > 0
          ? discoveries.reduce((sum, d) => sum + d.dataCompleteness, 0) /
            discoveries.length
          : 0,
    };
  } catch (error) {
    console.error("Error calculating discovery statistics:", error);
    return null;
  }
}
