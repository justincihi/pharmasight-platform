import { invokeLLM } from "./_core/llm";

/**
 * Retrosynthesis AI Module
 * Generates synthetic routes for target molecules using LLM-powered retrosynthesis
 */

export interface SynthesisStep {
  stepNumber: number;
  reaction: string;
  reagents: string[];
  conditions: string;
  yield: string;
  difficulty: "easy" | "moderate" | "difficult";
  estimatedCost: number;
  notes: string;
}

export interface SynthesisRoute {
  routeId: string;
  targetSmiles: string;
  targetName: string;
  totalSteps: number;
  overallYield: string;
  totalCost: number;
  feasibilityScore: number;
  difficulty: "easy" | "moderate" | "difficult";
  estimatedTime: string;
  steps: SynthesisStep[];
  startingMaterials: Array<{
    name: string;
    smiles: string;
    availability: "commercial" | "synthesize";
    cost: number;
  }>;
  summary: string;
}

/**
 * Generate synthesis routes for a target molecule
 */
export async function generateSynthesisRoutes(
  smiles: string,
  compoundName: string,
  numRoutes: number = 2
): Promise<SynthesisRoute[]> {
  try {
    const prompt = `You are an expert synthetic organic chemist. Generate ${numRoutes} practical synthesis routes for the following compound:

Compound Name: ${compoundName}
SMILES: ${smiles}

For each route, provide:
1. Starting materials (commercially available or easy to synthesize)
2. Step-by-step reactions with reagents and conditions
3. Expected yields for each step
4. Difficulty assessment (easy/moderate/difficult)
5. Estimated cost per step (in USD)
6. Overall feasibility score (0-100)
7. Estimated total synthesis time

Focus on:
- Using commercially available starting materials when possible
- Practical reaction conditions (avoid extreme temperatures/pressures)
- High-yielding, well-established reactions
- Cost-effective reagents
- Scalability to at least 1g scale

Return the response as a JSON array of synthesis routes.`;

    const response = await invokeLLM({
      messages: [
        {
          role: "system",
          content:
            "You are an expert synthetic organic chemist specializing in retrosynthesis and practical organic synthesis. Provide detailed, realistic synthesis routes.",
        },
        { role: "user", content: prompt },
      ],
      response_format: {
        type: "json_schema",
        json_schema: {
          name: "synthesis_routes",
          strict: true,
          schema: {
            type: "object",
            properties: {
              routes: {
                type: "array",
                items: {
                  type: "object",
                  properties: {
                    routeId: { type: "string" },
                    totalSteps: { type: "integer" },
                    overallYield: { type: "string" },
                    totalCost: { type: "number" },
                    feasibilityScore: { type: "integer" },
                    difficulty: {
                      type: "string",
                      enum: ["easy", "moderate", "difficult"],
                    },
                    estimatedTime: { type: "string" },
                    steps: {
                      type: "array",
                      items: {
                        type: "object",
                        properties: {
                          stepNumber: { type: "integer" },
                          reaction: { type: "string" },
                          reagents: {
                            type: "array",
                            items: { type: "string" },
                          },
                          conditions: { type: "string" },
                          yield: { type: "string" },
                          difficulty: {
                            type: "string",
                            enum: ["easy", "moderate", "difficult"],
                          },
                          estimatedCost: { type: "number" },
                          notes: { type: "string" },
                        },
                        required: [
                          "stepNumber",
                          "reaction",
                          "reagents",
                          "conditions",
                          "yield",
                          "difficulty",
                          "estimatedCost",
                          "notes",
                        ],
                        additionalProperties: false,
                      },
                    },
                    startingMaterials: {
                      type: "array",
                      items: {
                        type: "object",
                        properties: {
                          name: { type: "string" },
                          smiles: { type: "string" },
                          availability: {
                            type: "string",
                            enum: ["commercial", "synthesize"],
                          },
                          cost: { type: "number" },
                        },
                        required: ["name", "smiles", "availability", "cost"],
                        additionalProperties: false,
                      },
                    },
                    summary: { type: "string" },
                  },
                  required: [
                    "routeId",
                    "totalSteps",
                    "overallYield",
                    "totalCost",
                    "feasibilityScore",
                    "difficulty",
                    "estimatedTime",
                    "steps",
                    "startingMaterials",
                    "summary",
                  ],
                  additionalProperties: false,
                },
              },
            },
            required: ["routes"],
            additionalProperties: false,
          },
        },
      },
    });

    const content = response.choices[0]?.message?.content;
    if (!content || typeof content !== "string") {
      throw new Error("No response from LLM");
    }

    const parsed = JSON.parse(content);
    const routes: SynthesisRoute[] = parsed.routes.map((route: any) => ({
      ...route,
      targetSmiles: smiles,
      targetName: compoundName,
    }));

    return routes;
  } catch (error) {
    console.error("[Retrosynthesis] Error generating routes:", error);
    throw new Error(`Failed to generate synthesis routes: ${error}`);
  }
}

/**
 * Estimate reagent costs based on common lab pricing
 */
export function estimateReagentCost(reagent: string): number {
  // Simple cost estimation based on common reagents
  const costMap: Record<string, number> = {
    // Solvents (per liter)
    "dichloromethane": 50,
    "tetrahydrofuran": 60,
    "dimethylformamide": 70,
    "acetonitrile": 55,
    "methanol": 30,
    "ethanol": 25,
    "water": 5,
    
    // Common reagents (per 100g)
    "sodium borohydride": 150,
    "lithium aluminum hydride": 200,
    "palladium on carbon": 500,
    "triethylamine": 40,
    "pyridine": 35,
    "acetic anhydride": 45,
    
    // Expensive reagents
    "grubbs catalyst": 2000,
    "rhodium catalyst": 3000,
  };

  const reagentLower = reagent.toLowerCase();
  for (const [key, cost] of Object.entries(costMap)) {
    if (reagentLower.includes(key)) {
      return cost;
    }
  }

  // Default estimate
  return 100;
}
