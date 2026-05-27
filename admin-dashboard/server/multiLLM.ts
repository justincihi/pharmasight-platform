import { invokeLLM } from "./_core/llm";
import { GoogleGenerativeAI } from "@google/generative-ai";
import Anthropic from "@anthropic-ai/sdk";

/**
 * Multi-LLM Router
 * Unified interface for OpenAI, Gemini, Claude, Perplexity, and xAI (Grok)
 *
 * Required environment variables (set in .env, never committed to source control):
 *   OPENAI_API_KEY       – OpenAI (GPT-4o); falls back to built-in Forge proxy
 *   GEMINI_API_KEY       – Google Gemini
 *   ANTHROPIC_API_KEY    – Anthropic Claude
 *   PERPLEXITY_API_KEY   – Perplexity Sonar (SONAR_API_KEY is also accepted for backwards compatibility)
 *   XAI_API_KEY          – xAI Grok
 */

export type LLMProvider = "openai" | "gemini" | "claude" | "perplexity" | "xai";

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
 * OpenAI GPT-4o
 * Uses OPENAI_API_KEY if set, otherwise falls back to the built-in Forge proxy.
 */
async function callOpenAI(messages: LLMMessage[]): Promise<LLMResponse> {
  const directApiKey = process.env.OPENAI_API_KEY;

  if (directApiKey) {
    // Direct OpenAI API call
    const response = await fetch("https://api.openai.com/v1/chat/completions", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${directApiKey}`,
      },
      body: JSON.stringify({
        model: "gpt-4o",
        messages: messages.map((msg) => ({
          role: msg.role,
          content: msg.content,
        })),
        max_tokens: 4096,
      }),
    });

    if (!response.ok) {
      throw new Error(`OpenAI API error: ${response.statusText}`);
    }

    const data = await response.json();
    const content = data.choices?.[0]?.message?.content;
    if (typeof content !== "string") {
      throw new Error("OpenAI API returned an unexpected response structure");
    }
    return {
      content,
      provider: "openai",
      model: "gpt-4o",
    };
  }

  // Fallback: built-in Forge proxy
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
 * Perplexity Sonar (for research queries)
 * Reads PERPLEXITY_API_KEY, with SONAR_API_KEY as a fallback for backwards compatibility.
 */
async function callPerplexity(messages: LLMMessage[]): Promise<LLMResponse> {
  const apiKey = process.env.PERPLEXITY_API_KEY || process.env.SONAR_API_KEY;
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
 * xAI Grok
 * xAI exposes an OpenAI-compatible endpoint at https://api.x.ai/v1
 */
async function callXAI(messages: LLMMessage[]): Promise<LLMResponse> {
  const apiKey = process.env.XAI_API_KEY;
  if (!apiKey) {
    throw new Error("XAI_API_KEY not found in environment");
  }

  const response = await fetch("https://api.x.ai/v1/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      Authorization: `Bearer ${apiKey}`,
    },
    body: JSON.stringify({
      model: "grok-3-latest",
      messages: messages.map((msg) => ({
        role: msg.role,
        content: msg.content,
      })),
      max_tokens: 4096,
    }),
  });

  if (!response.ok) {
    throw new Error(`xAI API error: ${response.statusText}`);
  }

  const data = await response.json();
  const content = data.choices?.[0]?.message?.content;
  if (typeof content !== "string") {
    throw new Error("xAI API returned an unexpected response structure");
  }

  return {
    content,
    provider: "xai",
    model: "grok-3-latest",
  };
}

/**
 * Unified LLM call with automatic provider selection.
 * Accepts an explicit provider or auto-selects one based on the query content.
 * Falls back to OpenAI if the selected provider fails.
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
      provider = (process.env.PERPLEXITY_API_KEY || process.env.SONAR_API_KEY) ? "perplexity" : "gemini";
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
      case "xai":
        return await callXAI(messages);
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
 * Get available LLM providers based on configured environment variables.
 */
export function getAvailableProviders(): LLMProvider[] {
  const providers: LLMProvider[] = ["openai"]; // Always available via built-in or OPENAI_API_KEY

  if (process.env.GEMINI_API_KEY) providers.push("gemini");
  if (process.env.ANTHROPIC_API_KEY) providers.push("claude");
  if (process.env.PERPLEXITY_API_KEY || process.env.SONAR_API_KEY) providers.push("perplexity");
  if (process.env.XAI_API_KEY) providers.push("xai");

  return providers;
}

/**
 * Get the configuration status for every supported LLM provider.
 * Indicates whether the required API key environment variable is set.
 * Does NOT expose the actual key values.
 */
export function getProvidersStatus(): Record<LLMProvider, { configured: boolean; keyVar: string }> {
  return {
    openai: {
      configured: !!(process.env.OPENAI_API_KEY || process.env.BUILT_IN_FORGE_API_KEY),
      keyVar: "OPENAI_API_KEY",
    },
    gemini: {
      configured: !!process.env.GEMINI_API_KEY,
      keyVar: "GEMINI_API_KEY",
    },
    claude: {
      configured: !!process.env.ANTHROPIC_API_KEY,
      keyVar: "ANTHROPIC_API_KEY",
    },
    perplexity: {
      configured: !!(process.env.PERPLEXITY_API_KEY || process.env.SONAR_API_KEY),
      keyVar: "PERPLEXITY_API_KEY",
    },
    xai: {
      configured: !!process.env.XAI_API_KEY,
      keyVar: "XAI_API_KEY",
    },
  };
}
