import { readFileSync } from "fs";
import { getDb } from "./db";
import { analogDiscoveries } from "../drizzle/schema";

interface AnalogDiscovery {
  compound_id: string;
  parent_compound: string;
  smiles: string;
  transformation: string;
  confidence_score: number;
  similarity_score: number;
  safety_score: number;
  efficacy_score: number;
  drug_likeness_score: number;
  patent_score: number;
  market_value?: number;
  therapeutic_potential?: string;
  key_differences?: string;
  discovery_date: string;
}

export async function importAnalogDiscoveriesFromFile(filePath: string) {
  console.log(`[Import] Reading analog discoveries from: ${filePath}`);
  
  const fileContent = readFileSync(filePath, "utf-8");
  const data = JSON.parse(fileContent);
  
  // Handle both array format and object format
  const discoveries: AnalogDiscovery[] = Array.isArray(data) 
    ? data 
    : data.discoveries || [];
  
  console.log(`[Import] Found ${discoveries.length} analog discoveries to import`);
  
  const db = await getDb();
  if (!db) {
    throw new Error("Database not available");
  }
  
  let imported = 0;
  let skipped = 0;
  
  for (const discovery of discoveries) {
    try {
      // Determine patent status based on patent_score
      // Higher score = more likely to be patent-free
      let patentStatus: "patent-free" | "patent-opportunity" | "patented" | "unknown" = "patent-free";
      if (discovery.patent_score >= 90) {
        patentStatus = "patent-free";           // 90-100: Definitely patent-free
      } else if (discovery.patent_score >= 70) {
        patentStatus = "patent-opportunity";    // 70-89: Likely patent-free opportunity
      } else if (discovery.patent_score >= 50) {
        patentStatus = "patent-opportunity";    // 50-69: Possible patent-free
      } else {
        patentStatus = "unknown";               // <50: Unknown patent status
      }
      
      await db.insert(analogDiscoveries).values({
        compoundId: discovery.compound_id,
        compoundName: discovery.compound_id, // Use compound_id as name if not provided
        parentCompound: discovery.parent_compound,
        smiles: discovery.smiles,
        confidenceScore: discovery.confidence_score,
        similarityScore: discovery.similarity_score,
        safetyScore: discovery.safety_score,
        efficacyScore: discovery.efficacy_score,
        drugLikenessScore: discovery.drug_likeness_score,
        patentStatus,
        marketValue: discovery.market_value ? `$${discovery.market_value}M` : undefined,
        therapeuticPotential: discovery.therapeutic_potential,
        keyDifferences: discovery.key_differences,
        discoveredBy: "autonomous_system",
        discoveredAt: new Date(discovery.discovery_date),
      }).onDuplicateKeyUpdate({
        set: {
          confidenceScore: discovery.confidence_score,
          similarityScore: discovery.similarity_score,
          safetyScore: discovery.safety_score,
          efficacyScore: discovery.efficacy_score,
        },
      });
      
      imported++;
    } catch (error: any) {
      if (error.code === "ER_DUP_ENTRY") {
        skipped++;
      } else {
        console.error(`[Import] Error importing ${discovery.compound_id}:`, error.message);
      }
    }
  }
  
  console.log(`[Import] Complete: ${imported} imported, ${skipped} skipped (duplicates)`);
  
  return {
    total: discoveries.length,
    imported,
    skipped,
  };
}

// CLI usage
if (import.meta.url === `file://${process.argv[1]}`) {
  const filePath = process.argv[2] || "/home/ubuntu/pharmasight-platform/MASTER_ANALOG_DISCOVERIES.json";
  
  importAnalogDiscoveriesFromFile(filePath)
    .then((result) => {
      console.log(`\n✅ Import successful:`);
      console.log(`   Total: ${result.total}`);
      console.log(`   Imported: ${result.imported}`);
      console.log(`   Skipped: ${result.skipped}`);
      process.exit(0);
    })
    .catch((error) => {
      console.error(`\n❌ Import failed:`, error);
      process.exit(1);
    });
}
