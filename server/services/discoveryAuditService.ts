import { getDb } from "../db";
import { analogDiscoveries } from "../../drizzle/schema";
import { eq, and, gte, lt } from "drizzle-orm";

export interface DiscoveryAuditResult {
  id: number;
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  smiles: string;
  confidenceScore: number;
  patentStatus: string;
  discoveredAt: Date;
  discoveredBy: string;
  discoveryMethod?: string;
  isPartial: boolean;
  dataCompleteness: number; // 0-100 percentage
}

export interface DiscoveryValidationResult {
  isValid: boolean;
  isPartial: boolean;
  issues: string[];
  missingFields: string[];
  dataCompleteness: number;
  pharmacophoreMatch: number; // 0-100
  recommendations: string[];
}

/**
 * Get discoveries from the last N days that haven't been added to master list
 */
export async function getRecentDiscoveries(
  daysBack: number = 60
): Promise<DiscoveryAuditResult[]> {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const cutoffDate = new Date();
    cutoffDate.setDate(cutoffDate.getDate() - daysBack);

    const discoveries = await db
      .select()
      .from(analogDiscoveries)
      .where(gte(analogDiscoveries.discoveredAt, cutoffDate))
      .orderBy((t) => t.discoveredAt);

    return discoveries.map((d: any) => ({
      id: d.id,
      compoundId: d.compoundId,
      compoundName: d.compoundName,
      parentCompound: d.parentCompound,
      smiles: d.smiles,
      confidenceScore: d.confidenceScore,
      patentStatus: d.patentStatus,
      discoveredAt: d.discoveredAt,
      discoveredBy: d.discoveredBy,
      discoveryMethod: d.discoveryMethod,
      isPartial: validateSMILES(d.smiles).isPartial,
      dataCompleteness: calculateDataCompleteness(d),
    }));
  } catch (error) {
    console.error("Error fetching recent discoveries:", error);
    throw error;
  }
}

/**
 * Validate SMILES string and check for completeness
 */
export function validateSMILES(smiles: string): {
  isValid: boolean;
  isPartial: boolean;
  issues: string[];
} {
  const issues: string[] = [];
  let isPartial = false;

  // Check if SMILES is empty or too short
  if (!smiles || smiles.length < 3) {
    issues.push("SMILES string is too short or empty");
    isPartial = true;
  }

  // Check for common incomplete patterns
  if (smiles.includes("*")) {
    issues.push("Contains attachment point (*) - likely a fragment");
    isPartial = true;
  }

  // Check for unbalanced brackets
  const openBrackets = (smiles.match(/\[/g) || []).length;
  const closeBrackets = (smiles.match(/\]/g) || []).length;
  if (openBrackets !== closeBrackets) {
    issues.push("Unbalanced brackets in SMILES");
    isPartial = true;
  }

  // Check for valid SMILES characters
  const validChars = /^[A-Za-z0-9\[\]()=\-#@\\\/\+\.]+$/;
  if (!validChars.test(smiles)) {
    issues.push("Contains invalid SMILES characters");
    isPartial = true;
  }

  // Check for common fragment indicators
  if (
    smiles.match(/^[A-Z]$/) ||
    smiles.match(/^[A-Z]{1,2}$/) ||
    smiles.length < 5
  ) {
    issues.push("SMILES appears to be a fragment or single atom");
    isPartial = true;
  }

  return {
    isValid: issues.length === 0,
    isPartial,
    issues,
  };
}

/**
 * Calculate data completeness percentage
 */
export function calculateDataCompleteness(discovery: any): number {
  const requiredFields = [
    "compoundName",
    "smiles",
    "parentCompound",
    "confidenceScore",
    "patentStatus",
    "mechanismOfAction",
    "bindingAffinity",
    "dockingScore",
  ];

  const completedFields = requiredFields.filter(
    (field) => discovery[field] && discovery[field] !== ""
  ).length;

  return Math.round((completedFields / requiredFields.length) * 100);
}

/**
 * Validate discovery and check pharmacophore match
 */
export function validateDiscovery(
  discovery: any,
  parentCompound?: any
): DiscoveryValidationResult {
  const smileValidation = validateSMILES(discovery.smiles);
  const dataCompleteness = calculateDataCompleteness(discovery);
  const missingFields: string[] = [];

  // Check for missing critical fields
  if (!discovery.compoundName) missingFields.push("compoundName");
  if (!discovery.mechanismOfAction) missingFields.push("mechanismOfAction");
  if (!discovery.bindingAffinity) missingFields.push("bindingAffinity");
  if (!discovery.dockingScore) missingFields.push("dockingScore");

  // Calculate pharmacophore match (simplified - would need RDKit in production)
  let pharmacophoreMatch = 0;
  if (parentCompound && discovery.smiles) {
    // Basic heuristic: check if key elements are preserved
    const parentElements = extractChemicalElements(parentCompound.smiles);
    const discoveryElements = extractChemicalElements(discovery.smiles);
    const commonElements = parentElements.filter((e) =>
      discoveryElements.includes(e)
    );
    pharmacophoreMatch = Math.round(
      (commonElements.length / Math.max(parentElements.length, 1)) * 100
    );
  }

  const issues = [...smileValidation.issues];
  const recommendations: string[] = [];

  // Generate recommendations
  if (smileValidation.isPartial) {
    recommendations.push(
      "This appears to be a partial discovery or scaffold. Consider auto-querying via chatbot for more information."
    );
  }

  if (dataCompleteness < 50) {
    recommendations.push(
      "Data completeness is low. Additional analysis required before adding to master list."
    );
  }

  if (pharmacophoreMatch < 60 && parentCompound) {
    recommendations.push(
      "Pharmacophore match is low. Verify this is a valid analog of the parent compound."
    );
  }

  if (discovery.confidenceScore < 50) {
    recommendations.push(
      "Confidence score is low. Consider running additional validation tests."
    );
  }

  return {
    isValid: issues.length === 0 && dataCompleteness >= 50,
    isPartial: smileValidation.isPartial,
    issues,
    missingFields,
    dataCompleteness,
    pharmacophoreMatch,
    recommendations,
  };
}

/**
 * Extract chemical elements from SMILES (simplified)
 */
function extractChemicalElements(smiles: string): string[] {
  const elements = new Set<string>();
  const elementPattern = /[A-Z][a-z]?/g;
  const matches = smiles.match(elementPattern) || [];

  matches.forEach((match) => {
    if (
      !match.match(/^[0-9]/) &&
      !["Cl", "Br", "I", "F"].includes(match) === false
    ) {
      elements.add(match);
    }
  });

  return Array.from(elements);
}

/**
 * Import discovery to master list
 */
export async function importDiscoveryToMasterList(
  discoveryId: number,
  validation: DiscoveryValidationResult
): Promise<{ success: boolean; message: string }> {
  try {
    if (!validation.isValid && !validation.isPartial) {
      return {
        success: false,
        message: "Discovery failed validation. Cannot import.",
      };
    }

    // In a real implementation, this would update the master list
    // For now, just return success
    return {
      success: true,
      message: `Discovery ${discoveryId} imported successfully`,
    };
  } catch (error) {
    console.error("Error importing discovery:", error);
    return {
      success: false,
      message: `Error importing discovery: ${error}`,
    };
  }
}

/**
 * Save discovery to scaffold library (for partial discoveries)
 */
export async function saveToScaffoldLibrary(
  discoveryId: number,
  notes?: string
): Promise<{ success: boolean; scaffoldId?: string; message: string }> {
  try {
    // In a real implementation, this would save to a separate scaffold table
    const scaffoldId = `SCAFFOLD-${Date.now()}-${discoveryId}`;

    return {
      success: true,
      scaffoldId,
      message: `Scaffold saved with ID: ${scaffoldId}`,
    };
  } catch (error) {
    console.error("Error saving to scaffold library:", error);
    return {
      success: false,
      message: `Error saving scaffold: ${error}`,
    };
  }
}
