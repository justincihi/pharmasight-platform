import { invokeLLM } from './_core/llm';

export interface SynthesisStep {
  stepNumber: number;
  reaction: string;
  reagents: string[];
  conditions: string;
  estimatedYield: number; // percentage
  estimatedCost: number; // USD
  difficulty: 'easy' | 'moderate' | 'hard';
  notes: string;
}

export interface SynthesisRoute {
  targetCompound: string;
  smiles: string;
  totalSteps: number;
  totalEstimatedCost: number;
  totalEstimatedYield: number;
  estimatedTime: string;
  steps: SynthesisStep[];
  startingMaterials: string[];
  overallDifficulty: 'easy' | 'moderate' | 'hard';
  commercialAvailability: {
    material: string;
    available: boolean;
    supplier?: string;
    catalogPrice?: number;
  }[];
}

/**
 * Generate synthesis route using AI-powered retrosynthesis
 */
export async function generateSynthesisRoute(
  compoundName: string,
  smiles: string
): Promise<SynthesisRoute> {
  const prompt = `You are an expert synthetic organic chemist. Generate a detailed retrosynthetic analysis and forward synthesis route for the following compound:

**Target Compound:** ${compoundName}
**SMILES:** ${smiles}

Please provide:
1. A step-by-step synthesis route (3-7 steps recommended)
2. For each step, include:
   - Reaction type and mechanism
   - Required reagents and catalysts
   - Reaction conditions (temperature, solvent, time)
   - Estimated yield (%)
   - Estimated cost per step (USD)
   - Difficulty level (easy/moderate/hard)
   - Important notes or precautions

3. Starting materials (commercially available)
4. Overall synthesis difficulty
5. Estimated total time
6. Commercial availability of starting materials

Format your response as a structured JSON object matching this schema:
{
  "targetCompound": string,
  "smiles": string,
  "totalSteps": number,
  "totalEstimatedCost": number,
  "totalEstimatedYield": number,
  "estimatedTime": string,
  "steps": [
    {
      "stepNumber": number,
      "reaction": string,
      "reagents": string[],
      "conditions": string,
      "estimatedYield": number,
      "estimatedCost": number,
      "difficulty": "easy" | "moderate" | "hard",
      "notes": string
    }
  ],
  "startingMaterials": string[],
  "overallDifficulty": "easy" | "moderate" | "hard",
  "commercialAvailability": [
    {
      "material": string,
      "available": boolean,
      "supplier": string,
      "catalogPrice": number
    }
  ]
}`;

  try {
    const response = await invokeLLM({
      messages: [
        {
          role: 'system',
          content: 'You are an expert synthetic organic chemist specializing in retrosynthetic analysis and route optimization. Provide detailed, practical synthesis routes.',
        },
        {
          role: 'user',
          content: prompt,
        },
      ],
      response_format: {
        type: 'json_schema',
        json_schema: {
          name: 'synthesis_route',
          strict: true,
          schema: {
            type: 'object',
            properties: {
              targetCompound: { type: 'string' },
              smiles: { type: 'string' },
              totalSteps: { type: 'integer' },
              totalEstimatedCost: { type: 'number' },
              totalEstimatedYield: { type: 'number' },
              estimatedTime: { type: 'string' },
              steps: {
                type: 'array',
                items: {
                  type: 'object',
                  properties: {
                    stepNumber: { type: 'integer' },
                    reaction: { type: 'string' },
                    reagents: { type: 'array', items: { type: 'string' } },
                    conditions: { type: 'string' },
                    estimatedYield: { type: 'number' },
                    estimatedCost: { type: 'number' },
                    difficulty: { type: 'string', enum: ['easy', 'moderate', 'hard'] },
                    notes: { type: 'string' },
                  },
                  required: ['stepNumber', 'reaction', 'reagents', 'conditions', 'estimatedYield', 'estimatedCost', 'difficulty', 'notes'],
                  additionalProperties: false,
                },
              },
              startingMaterials: { type: 'array', items: { type: 'string' } },
              overallDifficulty: { type: 'string', enum: ['easy', 'moderate', 'hard'] },
              commercialAvailability: {
                type: 'array',
                items: {
                  type: 'object',
                  properties: {
                    material: { type: 'string' },
                    available: { type: 'boolean' },
                    supplier: { type: 'string' },
                    catalogPrice: { type: 'number' },
                  },
                  required: ['material', 'available'],
                  additionalProperties: false,
                },
              },
            },
            required: ['targetCompound', 'smiles', 'totalSteps', 'totalEstimatedCost', 'totalEstimatedYield', 'estimatedTime', 'steps', 'startingMaterials', 'overallDifficulty', 'commercialAvailability'],
            additionalProperties: false,
          },
        },
      },
    });

    const content = response.choices[0].message.content;
    if (!content) {
      throw new Error('No response from LLM');
    }

    const contentStr = typeof content === 'string' ? content : JSON.stringify(content);
    const route = JSON.parse(contentStr) as SynthesisRoute;
    return route;
  } catch (error: any) {
    console.error('Synthesis route generation failed:', error);
    throw new Error(`Failed to generate synthesis route: ${error.message}`);
  }
}

/**
 * Optimize existing synthesis route for cost or yield
 */
export async function optimizeSynthesisRoute(
  route: SynthesisRoute,
  optimizationGoal: 'cost' | 'yield' | 'time'
): Promise<SynthesisRoute> {
  const prompt = `Optimize the following synthesis route for ${optimizationGoal}:

**Current Route:**
${JSON.stringify(route, null, 2)}

Please provide an optimized version that:
${optimizationGoal === 'cost' ? '- Minimizes total cost while maintaining reasonable yields\n- Suggests cheaper reagents or alternative reactions' : ''}
${optimizationGoal === 'yield' ? '- Maximizes overall yield\n- Suggests higher-yielding conditions or protecting group strategies' : ''}
${optimizationGoal === 'time' ? '- Minimizes total synthesis time\n- Suggests faster reactions or one-pot procedures' : ''}

Return the optimized route in the same JSON format.`;

  try {
    const response = await invokeLLM({
      messages: [
        {
          role: 'system',
          content: 'You are an expert synthetic organic chemist specializing in route optimization.',
        },
        {
          role: 'user',
          content: prompt,
        },
      ],
    });

    const content = response.choices[0].message.content;
    if (!content) {
      throw new Error('No response from LLM');
    }

    const contentStr = typeof content === 'string' ? content : JSON.stringify(content);

    // Extract JSON from response
    const jsonMatch = contentStr.match(/\{[\s\S]*\}/);
    if (!jsonMatch) {
      throw new Error('No JSON found in response');
    }

    const optimizedRoute = JSON.parse(jsonMatch[0]) as SynthesisRoute;
    return optimizedRoute;
  } catch (error: any) {
    console.error('Route optimization failed:', error);
    throw new Error(`Failed to optimize route: ${error.message}`);
  }
}
