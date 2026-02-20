/**
 * FDA Orange Book API Integration
 * Provides patent and regulatory approval information for pharmaceutical compounds
 */

interface FDAProduct {
  activeIngredient: string;
  tradeName: string;
  applicationType: string;
  approvalDate: string;
  patents: Array<{
    patentNumber: string;
    patentExpireDate: string;
    drugProductFlag: string;
    drugSubstanceFlag: string;
  }>;
  exclusivities: Array<{
    exclusivityCode: string;
    exclusivityDate: string;
  }>;
}

interface FDASearchResult {
  found: boolean;
  products: FDAProduct[];
  patentStatus: "patent-free" | "patented" | "unknown";
  patentNumbers: string[];
  approvalStatus: string | null;
}

/**
 * Search FDA Orange Book for compound information
 */
export async function searchFDAOrangeBook(compoundName: string): Promise<FDASearchResult> {
  try {
    // FDA OpenFDA API endpoint for drug products
    const url = `https://api.fda.gov/drug/drugsfda.json?search=openfda.brand_name:"${encodeURIComponent(compoundName)}"&limit=10`;
    
    const response = await fetch(url);
    
    if (!response.ok) {
      console.warn(`FDA API returned ${response.status} for ${compoundName}`);
      return {
        found: false,
        products: [],
        patentStatus: "unknown",
        patentNumbers: [],
        approvalStatus: null,
      };
    }

    const data = await response.json();
    
    if (!data.results || data.results.length === 0) {
      return {
        found: false,
        products: [],
        patentStatus: "unknown",
        patentNumbers: [],
        approvalStatus: null,
      };
    }

    // Parse FDA results
    const products: FDAProduct[] = data.results.map((result: any) => ({
      activeIngredient: result.products?.[0]?.active_ingredients?.[0]?.name || "Unknown",
      tradeName: result.products?.[0]?.brand_name || "Unknown",
      applicationType: result.application_number || "Unknown",
      approvalDate: result.submissions?.[0]?.submission_status_date || "Unknown",
      patents: [], // FDA API doesn't directly provide patent info in this endpoint
      exclusivities: [],
    }));

    // Determine patent status (simplified - would need Orange Book patents endpoint for full data)
    const hasApproval = products.length > 0;
    const patentStatus: "patent-free" | "patented" | "unknown" = hasApproval ? "patented" : "unknown";

    return {
      found: true,
      products,
      patentStatus,
      patentNumbers: [],
      approvalStatus: hasApproval ? "FDA Approved" : null,
    };
  } catch (error) {
    console.error("FDA Orange Book API error:", error);
    return {
      found: false,
      products: [],
      patentStatus: "unknown",
      patentNumbers: [],
      approvalStatus: null,
    };
  }
}

/**
 * Check if a compound is patent-free based on FDA data
 */
export async function checkPatentStatus(compoundName: string): Promise<{
  status: "patent-free" | "patented" | "unknown";
  patents: string[];
  message: string;
}> {
  const result = await searchFDAOrangeBook(compoundName);
  
  if (!result.found) {
    return {
      status: "patent-free",
      patents: [],
      message: "No FDA approval found - likely patent-free opportunity",
    };
  }

  return {
    status: result.patentStatus,
    patents: result.patentNumbers,
    message: result.approvalStatus || "Patent status unclear",
  };
}
