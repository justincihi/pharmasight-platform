/**
 * Chatbot Analog Generation Handler
 * Processes natural language queries for analog generation and PubChem searches
 */

import { confirmAndFetchSimilars, screenAndFlagForMasterList } from "./pubchemSearchService";
import { generateBRICSAnalogs, rankAnalogs, filterDrugLikeAnalogs } from "./bricsAnalogService";

export interface ChatbotAnalogQuery {
  type: "pubchem_search" | "brics_generation" | "smiles_extraction" | "patent_check";
  compound: string;
  similarity?: number;
  maxResults?: number;
  minConfidence?: number;
  patentFree?: boolean;
}

export interface ChatbotAnalogResponse {
  success: boolean;
  message: string;
  analogs?: any[];
  count?: number;
  saveOptions?: {
    toMasterList: boolean;
    toScaffoldLibrary: boolean;
  };
}

/**
 * Parse natural language query and extract analog generation intent
 */
export function parseAnalogQuery(userMessage: string): ChatbotAnalogQuery | null {
  const message = userMessage.toLowerCase();

  // PubChem search patterns
  if (
    message.includes("find") &&
    (message.includes("analog") || message.includes("similar"))
  ) {
    const compoundMatch = userMessage.match(
      /(?:find|search|look for).*?(?:analog|similar|compound)s?\s+(?:of|for|to)\s+([A-Za-z0-9\-]+)/i
    );
    if (compoundMatch) {
      return {
        type: "pubchem_search",
        compound: compoundMatch[1],
        similarity: 0.7,
        maxResults: 25,
        minConfidence: 85,
        patentFree: message.includes("patent-free") || message.includes("patent free"),
      };
    }
  }

  // BRICS generation patterns
  if (
    message.includes("generate") &&
    (message.includes("analog") || message.includes("variant"))
  ) {
    const compoundMatch = userMessage.match(
      /(?:generate|create|make).*?(?:analog|variant|derivative)s?\s+(?:of|for|from)\s+([A-Za-z0-9\-]+)/i
    );
    if (compoundMatch) {
      return {
        type: "brics_generation",
        compound: compoundMatch[1],
        maxResults: 25,
      };
    }
  }

  // SMILES extraction patterns
  if (message.includes("add") && message.includes("smiles")) {
    const smilesMatch = userMessage.match(/smiles[:\s]+([A-Za-z0-9()[\]\-=#\\/@+.]+)/i);
    if (smilesMatch) {
      return {
        type: "smiles_extraction",
        compound: smilesMatch[1],
      };
    }
  }

  // Patent check patterns
  if (message.includes("patent") && message.includes("check")) {
    const compoundMatch = userMessage.match(
      /(?:check|verify).*?patent.*?(?:for|of|on)\s+([A-Za-z0-9\-]+)/i
    );
    if (compoundMatch) {
      return {
        type: "patent_check",
        compound: compoundMatch[1],
      };
    }
  }

  return null;
}

/**
 * Handle PubChem search query
 */
export async function handlePubChemSearch(
  compound: string,
  similarity: number = 0.7,
  maxResults: number = 25,
  minConfidence: number = 85,
  patentFree: boolean = true
): Promise<ChatbotAnalogResponse> {
  try {
    const { parentSmiles, similar } = await confirmAndFetchSimilars(
      compound,
      similarity,
      maxResults
    );

    // Filter by confidence and patent status
    let filtered = similar.filter((hit) => hit.tanimoto * 100 >= minConfidence);
    if (patentFree) {
      filtered = filtered.filter((hit) => hit.patentFree);
    }

    return {
      success: true,
      message: `Found ${filtered.length} ${patentFree ? "patent-free " : ""}analogs of ${compound} with ${minConfidence}%+ confidence`,
      analogs: filtered,
      count: filtered.length,
      saveOptions: {
        toMasterList: true,
        toScaffoldLibrary: true,
      },
    };
  } catch (error) {
    return {
      success: false,
      message: `Failed to search for analogs: ${error instanceof Error ? error.message : "Unknown error"}`,
    };
  }
}

/**
 * Handle BRICS generation query
 */
export async function handleBRICSGeneration(
  compound: string,
  maxResults: number = 25
): Promise<ChatbotAnalogResponse> {
  try {
    // First, confirm compound via PubChem
    const { parentSmiles } = await confirmAndFetchSimilars(compound, 1.0, 1);

    // Generate BRICS analogs
    const analogs = await generateBRICSAnalogs(parentSmiles, maxResults);

    // Filter and rank
    const drugLike = filterDrugLikeAnalogs(analogs);
    const ranked = rankAnalogs(drugLike);

    return {
      success: true,
      message: `Generated ${ranked.length} novel analogs of ${compound} using BRICS fragmentation`,
      analogs: ranked,
      count: ranked.length,
      saveOptions: {
        toMasterList: true,
        toScaffoldLibrary: true,
      },
    };
  } catch (error) {
    return {
      success: false,
      message: `Failed to generate analogs: ${error instanceof Error ? error.message : "Unknown error"}`,
    };
  }
}

/**
 * Handle SMILES extraction from conversation
 */
export async function handleSMILESExtraction(
  smiles: string
): Promise<ChatbotAnalogResponse> {
  try {
    // Validate SMILES
    const isValid = /^[A-Za-z0-9()[\]\-=#\\/@+.]+$/.test(smiles);
    if (!isValid) {
      return {
        success: false,
        message: `Invalid SMILES string: ${smiles}. Please check the format.`,
      };
    }

    return {
      success: true,
      message: `SMILES extracted: ${smiles}. Would you like to generate analogs or search for similar compounds?`,
      analogs: [{ smiles, name: "Extracted SMILES" }],
      saveOptions: {
        toMasterList: true,
        toScaffoldLibrary: true,
      },
    };
  } catch (error) {
    return {
      success: false,
      message: `Failed to process SMILES: ${error instanceof Error ? error.message : "Unknown error"}`,
    };
  }
}

/**
 * Main chatbot query handler
 */
export async function handleAnalogQuery(
  userMessage: string
): Promise<ChatbotAnalogResponse> {
  const query = parseAnalogQuery(userMessage);

  if (!query) {
    return {
      success: false,
      message: "I didn't understand your request. Try asking me to:\n- Find analogs of [compound]\n- Generate variants of [compound]\n- Check patents for [compound]\n- Add SMILES: [SMILES string]",
    };
  }

  switch (query.type) {
    case "pubchem_search":
      return handlePubChemSearch(
        query.compound,
        query.similarity,
        query.maxResults,
        query.minConfidence,
        query.patentFree
      );

    case "brics_generation":
      return handleBRICSGeneration(query.compound, query.maxResults);

    case "smiles_extraction":
      return handleSMILESExtraction(query.compound);

    case "patent_check":
      return {
        success: true,
        message: `Patent check for ${query.compound} - would query USPTO and Google Patents databases`,
      };

    default:
      return {
        success: false,
        message: "Unknown query type",
      };
  }
}

/**
 * Extract SMILES from conversation history
 */
export function extractSMILESFromHistory(conversationHistory: any[]): string[] {
  const smiles: string[] = [];
  const smilesPattern = /([A-Za-z0-9()[\]\-=#\\/@+.]{10,})/g;

  for (const message of conversationHistory) {
    const content = message.content || message.text || "";
    const matches = content.match(smilesPattern);
    if (matches) {
      smiles.push(...matches);
    }
  }

  return Array.from(new Set(smiles)); // Remove duplicates
}
