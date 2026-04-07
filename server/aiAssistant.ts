/**
 * PharmaSight AI Assistant
 * 
 * LLM-powered chatbot for natural language queries about analog discoveries
 * Provides structure-activity relationship analysis and research insights
 */

import { invokeLLM } from './_core/llm';
import { getDb } from './db';
import { analogDiscoveries, metabolites } from '../drizzle/schema';
import { desc, sql } from 'drizzle-orm';

interface ChatMessage {
  role: 'system' | 'user' | 'assistant';
  content: string;
}

interface AIAssistantResponse {
  message: string;
  data?: any;
  suggestions?: string[];
}

/**
 * Process user query and generate AI response
 */
export async function processAIQuery(
  userQuery: string,
  conversationHistory: ChatMessage[] = []
): Promise<AIAssistantResponse> {
  try {
    // Step 1: Analyze query intent and extract parameters
    const queryIntent = await analyzeQueryIntent(userQuery);
    
    // Step 2: Fetch relevant data from database
    const relevantData = await fetchRelevantData(queryIntent);
    
    // Step 3: Generate AI response with context
    const response = await generateAIResponse(userQuery, relevantData, conversationHistory);
    
    // Step 4: Generate follow-up suggestions
    const suggestions = generateSuggestions(queryIntent);

    return {
      message: response,
      data: relevantData,
      suggestions,
    };
  } catch (error) {
    console.error('AI Assistant error:', error);
    return {
      message: "I apologize, but I encountered an error processing your query. Please try rephrasing your question or contact support if the issue persists.",
      suggestions: [
        "What are the top 5 analogs discovered this week?",
        "Show me patent-free compounds with >85% confidence",
        "Compare ADMET profiles of recent discoveries",
      ],
    };
  }
}

/**
 * Analyze user query to determine intent and extract parameters
 */
async function analyzeQueryIntent(query: string): Promise<any> {
  const lowerQuery = query.toLowerCase();
  
  return {
    type: detectQueryType(lowerQuery),
    filters: extractFilters(lowerQuery),
    timeframe: extractTimeframe(lowerQuery),
    sortBy: extractSortCriteria(lowerQuery),
    limit: extractLimit(lowerQuery),
  };
}

function detectQueryType(query: string): string {
  if (query.includes('top') || query.includes('best') || query.includes('highest')) {
    return 'ranking';
  }
  if (query.includes('compare') || query.includes('vs') || query.includes('versus')) {
    return 'comparison';
  }
  if (query.includes('patent') || query.includes('ip') || query.includes('freedom')) {
    return 'patent';
  }
  if (query.includes('admet') || query.includes('safety') || query.includes('toxicity')) {
    return 'admet';
  }
  if (query.includes('synthesis') || query.includes('cost') || query.includes('route')) {
    return 'synthesis';
  }
  if (query.includes('metabolite') || query.includes('metabolism')) {
    return 'metabolite';
  }
  return 'general';
}

function extractFilters(query: string): any {
  const filters: any = {};
  
  // Extract confidence threshold
  const confidenceMatch = query.match(/(\d+)%?\s*confidence/);
  if (confidenceMatch) {
    filters.minConfidence = parseInt(confidenceMatch[1]);
  }
  
  // Extract patent-free requirement
  if (query.includes('patent-free') || query.includes('patent free')) {
    filters.patentFree = true;
  }
  
  // Extract therapeutic area
  if (query.includes('cns') || query.includes('neurological')) {
    filters.therapeuticArea = 'CNS Disorders';
  }
  
  return filters;
}

function extractTimeframe(query: string): string | null {
  if (query.includes('today')) return 'today';
  if (query.includes('week') || query.includes('7 days')) return 'week';
  if (query.includes('month') || query.includes('30 days')) return 'month';
  if (query.includes('year')) return 'year';
  return null;
}

function extractSortCriteria(query: string): string {
  if (query.includes('confidence')) return 'confidenceScore';
  if (query.includes('safety')) return 'safetyScore';
  if (query.includes('efficacy')) return 'efficacyScore';
  if (query.includes('recent') || query.includes('latest') || query.includes('new')) return 'createdAt';
  return 'confidenceScore';
}

function extractLimit(query: string): number {
  const numberMatch = query.match(/top\s+(\d+)|(\d+)\s+analogs?/);
  if (numberMatch) {
    return parseInt(numberMatch[1] || numberMatch[2]);
  }
  return 5;
}

/**
 * Fetch relevant data from database based on query intent
 */
async function fetchRelevantData(queryIntent: any): Promise<any> {
  const db = await getDb();
  if (!db) return null;

  try {
    let query = db.select().from(analogDiscoveries);

    // Apply filters
    const conditions = [];
    if (queryIntent.filters.minConfidence) {
      conditions.push(sql`${analogDiscoveries.confidenceScore} >= ${queryIntent.filters.minConfidence}`);
    }
    // Note: Additional filters like patentFree and therapeuticArea can be added when schema is updated

    // Apply timeframe filter
    if (queryIntent.timeframe) {
      const timeCondition = getTimeframeCondition(queryIntent.timeframe);
      if (timeCondition) {
        conditions.push(timeCondition);
      }
    }

    // Build query with conditions
    if (conditions.length > 0) {
      query = query.where(sql.join(conditions, sql` AND `)) as any;
    }

    // Apply sorting
    query = query.orderBy(desc(analogDiscoveries.confidenceScore)) as any;

    // Apply limit
    query = query.limit(queryIntent.limit) as any;

    const results = await query;
    return results;
  } catch (error) {
    console.error('Data fetch error:', error);
    return null;
  }
}

function getTimeframeCondition(timeframe: string): any {
  const now = new Date();
  let startDate: Date;

  switch (timeframe) {
    case 'today':
      startDate = new Date(now.setHours(0, 0, 0, 0));
      break;
    case 'week':
      startDate = new Date(now.setDate(now.getDate() - 7));
      break;
    case 'month':
      startDate = new Date(now.setMonth(now.getMonth() - 1));
      break;
    case 'year':
      startDate = new Date(now.setFullYear(now.getFullYear() - 1));
      break;
    default:
      return null;
  }

  return sql`${analogDiscoveries.createdAt} >= ${startDate}`;
}

/**
 * Generate AI response using LLM with context
 */
async function generateAIResponse(
  userQuery: string,
  data: any,
  conversationHistory: ChatMessage[]
): Promise<string> {
  const systemPrompt = `You are PharmaSight AI Assistant, an expert in pharmaceutical research and drug discovery. You help researchers analyze analog discoveries, interpret ADMET data, and provide insights on structure-activity relationships.

Your responses should be:
- Professional and scientifically accurate
- Concise but informative
- Data-driven when possible
- Helpful for decision-making

When presenting analog data, highlight key metrics like confidence scores, safety profiles, and patent status.`;

  const dataContext = data && data.length > 0
    ? `\n\nRelevant analog data:\n${JSON.stringify(data.slice(0, 3), null, 2)}`
    : '';

  const messages: ChatMessage[] = [
    { role: 'system', content: systemPrompt },
    ...conversationHistory,
    { role: 'user', content: userQuery + dataContext },
  ];

  const response = await invokeLLM({
    messages: messages as any,
  });

  const content = response.choices[0].message.content;
  return typeof content === 'string' ? content : "I couldn't generate a response. Please try again.";
}

/**
 * Generate follow-up suggestions based on query type
 */
function generateSuggestions(queryIntent: any): string[] {
  const suggestions: Record<string, string[]> = {
    ranking: [
      "Compare ADMET profiles of these top analogs",
      "Show synthesis routes for the highest-confidence compounds",
      "Check patent status for these discoveries",
    ],
    comparison: [
      "Generate a detailed PDF report for these analogs",
      "Show metabolite predictions for comparison",
      "Run docking analysis on all targets",
    ],
    patent: [
      "Find similar patent-free compounds",
      "Show freedom-to-operate analysis",
      "Estimate IP value for these discoveries",
    ],
    admet: [
      "Compare with known drugs in the same class",
      "Show metabolite safety profiles",
      "Identify potential toxicity risks",
    ],
    synthesis: [
      "Compare synthesis costs across analogs",
      "Show alternative synthesis routes",
      "Estimate production feasibility",
    ],
    metabolite: [
      "Show parent-metabolite ADMET comparison",
      "Identify active metabolites",
      "Predict metabolic stability",
    ],
    general: [
      "What are the top 5 analogs discovered this week?",
      "Show me patent-free compounds with >85% confidence",
      "Compare ADMET profiles of recent discoveries",
    ],
  };

  return suggestions[queryIntent.type] || suggestions.general;
}

/**
 * Get sample questions for user guidance
 */
export function getSampleQuestions(): string[] {
  return [
    "What are the top 5 analogs discovered this week?",
    "Show me patent-free compounds with >85% confidence",
    "Compare ADMET profiles of recent discoveries",
    "Which analogs have the best safety scores?",
    "Show synthesis routes for high-confidence CNS compounds",
    "What metabolites were predicted for the latest discoveries?",
  ];
}
