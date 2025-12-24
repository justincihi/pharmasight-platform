import { invokeLLM } from "./_core/llm";
import { GoogleGenerativeAI } from "@google/generative-ai";
import Anthropic from "@anthropic-ai/sdk";

/**
 * Multi-LLM Router
 * Unified interface for OpenAI, Gemini, Claude, and Perplexity
 */

export type LLMProvider = "openai" | "gemini" | "claude" | "perplexity";

export interface LLMMessage {
  role: "system" | "user" | "assistant";
  content: string;
}

export interface LLMResponse {
  content: string;
  provider: LLMProvider;
  model: string;
}

/**
 * OpenAI GPT-4 (via built-in invokeLLM)
 */
async function callOpenAI(messages: LLMMessage[]): Promise<LLMResponse> {
  const response = await invokeLLM({
    messages: messages.map((msg) => ({
      role: msg.role,
      content: msg.content,
    })),
  });

  const messageContent = response.choices[0]?.message?.content;
  const content = typeof messageContent === 'string' ? messageContent : '';

  return {
    content,
    provider: "openai",
    model: "gpt-4o",
  };
}

/**
 * Google Gemini
 */
async function callGemini(messages: LLMMessage[]): Promise<LLMResponse> {
  const apiKey = process.env.GEMINI_API_KEY;
  if (!apiKey) {
    throw new Error("GEMINI_API_KEY not found in environment");
  }

  const genAI = new GoogleGenerativeAI(apiKey);
  const model = genAI.getGenerativeModel({ model: "gemini-2.0-flash-exp" });

  // Convert messages to Gemini format
  const systemMessage = messages.find((m) => m.role === "system");
  const chatMessages = messages.filter((m) => m.role !== "system");

  const chat = model.startChat({
    history: chatMessages.slice(0, -1).map((msg) => ({
      role: msg.role === "assistant" ? "model" : "user",
      parts: [{ text: msg.content }],
    })),
    systemInstruction: systemMessage?.content,
  });

  const lastMessage = chatMessages[chatMessages.length - 1];
  const result = await chat.sendMessage(lastMessage?.content || "");
  const response = await result.response;

  return {
    content: response.text(),
    provider: "gemini",
    model: "gemini-2.0-flash-exp",
  };
}

/**
 * Anthropic Claude
 */
async function callClaude(messages: LLMMessage[]): Promise<LLMResponse> {
  const apiKey = process.env.ANTHROPIC_API_KEY;
  if (!apiKey) {
    throw new Error("ANTHROPIC_API_KEY not found in environment");
  }

  const anthropic = new Anthropic({
    apiKey,
  });

  // Extract system message
  const systemMessage = messages.find((m) => m.role === "system")?.content || "";
  const chatMessages = messages
    .filter((m) => m.role !== "system")
    .map((msg) => ({
      role: msg.role as "user" | "assistant",
      content: msg.content,
    }));

  const response = await anthropic.messages.create({
    model: "claude-3-5-sonnet-20241022",
    max_tokens: 4096,
    system: systemMessage,
    messages: chatMessages,
  });

  const content =
    response.content[0]?.type === "text" ? response.content[0].text : "";

  return {
    content,
    provider: "claude",
    model: "claude-3-5-sonnet-20241022",
  };
}

/**
 * Perplexity (for research queries)
 */
async function callPerplexity(messages: LLMMessage[]): Promise<LLMResponse> {
  const apiKey = process.env.PERPLEXITY_API_KEY;
  if (!apiKey) {
    throw new Error("PERPLEXITY_API_KEY not found in environment");
  }

  const response = await fetch("https://api.perplexity.ai/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      Authorization: `Bearer ${apiKey}`,
    },
    body: JSON.stringify({
      model: "sonar",
      messages: messages.map((msg) => ({
        role: msg.role,
        content: msg.content,
      })),
    }),
  });

  if (!response.ok) {
    throw new Error(`Perplexity API error: ${response.statusText}`);
  }

  const data = await response.json();

  return {
    content: data.choices[0]?.message?.content || "",
    provider: "perplexity",
    model: "sonar",
  };
}

/**
 * Unified LLM call with automatic provider selection
 */
export async function callLLM(
  messages: LLMMessage[],
  provider?: LLMProvider
): Promise<LLMResponse> {
  // Auto-select provider based on query type if not specified
  if (!provider) {
    const lastMessage = messages[messages.length - 1]?.content.toLowerCase() || "";

    // Use Perplexity for research queries
    if (
      lastMessage.includes("research") ||
      lastMessage.includes("latest") ||
      lastMessage.includes("recent studies") ||
      lastMessage.includes("clinical trials")
    ) {
      provider = process.env.PERPLEXITY_API_KEY ? "perplexity" : "gemini";
    }
    // Use Gemini for multimodal or long context
    else if (
      lastMessage.includes("analyze") ||
      lastMessage.includes("compare") ||
      lastMessage.length > 1000
    ) {
      provider = "gemini";
    }
    // Use Claude for deep reasoning
    else if (
      lastMessage.includes("explain") ||
      lastMessage.includes("why") ||
      lastMessage.includes("mechanism")
    ) {
      provider = "claude";
    }
    // Default to OpenAI
    else {
      provider = "openai";
    }
  }

  try {
    switch (provider) {
      case "openai":
        return await callOpenAI(messages);
      case "gemini":
        return await callGemini(messages);
      case "claude":
        return await callClaude(messages);
      case "perplexity":
        return await callPerplexity(messages);
      default:
        return await callOpenAI(messages);
    }
  } catch (error: any) {
    console.error(`[MultiLLM] Error calling ${provider}:`, error.message);

    // Fallback to OpenAI if provider fails
    if (provider !== "openai") {
      console.log("[MultiLLM] Falling back to OpenAI");
      return await callOpenAI(messages);
    }

    throw error;
  }
}

/**
 * Get available LLM providers
 */
export function getAvailableProviders(): LLMProvider[] {
  const providers: LLMProvider[] = ["openai"]; // Always available via built-in

  if (process.env.GEMINI_API_KEY) providers.push("gemini");
  if (process.env.ANTHROPIC_API_KEY) providers.push("claude");
  if (process.env.PERPLEXITY_API_KEY) providers.push("perplexity");

  return providers;
}
