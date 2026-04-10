import * as fs from "fs/promises";
import * as path from "path";

export interface ChatMessageLog {
  timestamp: Date;
  role: "user" | "assistant" | "system";
  content: string;
  metadata?: {
    executionTime?: number;
    tokensUsed?: number;
    model?: string;
    error?: string;
  };
}

const LOGS_DIR = path.join(process.cwd(), "logs", "conversations");

/**
 * Initialize logs directory
 */
async function ensureLogsDirectory(): Promise<void> {
  try {
    await fs.mkdir(LOGS_DIR, { recursive: true });
  } catch (error) {
    console.error("Error creating logs directory:", error);
  }
}

/**
 * Save conversation to markdown file
 */
export async function saveConversationToMarkdown(
  userId: string,
  userName: string,
  messages: ChatMessageLog[],
  topic?: string
): Promise<string> {
  try {
    await ensureLogsDirectory();

    const timestamp = new Date().toISOString().replace(/[:.]/g, "-");
    const filename = `conversation-${userId}-${timestamp}.md`;
    const filepath = path.join(LOGS_DIR, filename);

    // Calculate duration
    const startTime = messages[0]?.timestamp || new Date();
    const endTime = messages[messages.length - 1]?.timestamp || new Date();
    const duration = Math.round(
      (endTime.getTime() - startTime.getTime()) / 1000
    );

    // Build markdown content
    let markdown = `# Conversation Log\n\n`;
    markdown += `**User:** ${userName} (ID: ${userId})\n`;
    markdown += `**Date:** ${new Date().toLocaleString()}\n`;
    markdown += `**Duration:** ${duration} seconds\n`;
    markdown += `**Total Messages:** ${messages.length}\n`;
    if (topic) markdown += `**Topic:** ${topic}\n`;
    markdown += `\n---\n\n`;

    // Add messages
    messages.forEach((msg, idx) => {
      markdown += `## Message ${idx + 1}\n\n`;
      markdown += `**Role:** ${msg.role}\n`;
      markdown += `**Time:** ${msg.timestamp.toLocaleString()}\n`;

      if (msg.metadata) {
        markdown += `**Metadata:**\n`;
        if (msg.metadata.executionTime)
          markdown += `- Execution Time: ${msg.metadata.executionTime}ms\n`;
        if (msg.metadata.tokensUsed)
          markdown += `- Tokens Used: ${msg.metadata.tokensUsed}\n`;
        if (msg.metadata.model)
          markdown += `- Model: ${msg.metadata.model}\n`;
        if (msg.metadata.error)
          markdown += `- Error: ${msg.metadata.error}\n`;
        markdown += `\n`;
      }

      markdown += `**Content:**\n\n${msg.content}\n\n`;
      markdown += `---\n\n`;
    });

    // Write to file
    await fs.writeFile(filepath, markdown, "utf-8");

    console.log(`Conversation saved to: ${filepath}`);
    return filepath;
  } catch (error) {
    console.error("Error saving conversation to markdown:", error);
    throw error;
  }
}

/**
 * Get all conversation logs from file system
 */
export async function getAllConversationLogs(): Promise<string[]> {
  try {
    await ensureLogsDirectory();
    const files = await fs.readdir(LOGS_DIR);
    return files
      .filter((f) => f.endsWith(".md"))
      .map((f) => path.join(LOGS_DIR, f));
  } catch (error) {
    console.error("Error reading conversation logs:", error);
    return [];
  }
}

/**
 * Read a specific conversation log
 */
export async function readConversationLog(filepath: string): Promise<string> {
  try {
    return await fs.readFile(filepath, "utf-8");
  } catch (error) {
    console.error("Error reading conversation log:", error);
    throw error;
  }
}

/**
 * Search conversation logs by keyword
 */
export async function searchConversationLogs(
  keyword: string
): Promise<{ filepath: string; matches: number }[]> {
  try {
    const logs = await getAllConversationLogs();
    const results: { filepath: string; matches: number }[] = [];

    for (const filepath of logs) {
      const content = await readConversationLog(filepath);
      const matches = (content.match(new RegExp(keyword, "gi")) || []).length;
      if (matches > 0) {
        results.push({ filepath, matches });
      }
    }

    return results.sort((a, b) => b.matches - a.matches);
  } catch (error) {
    console.error("Error searching conversation logs:", error);
    return [];
  }
}

/**
 * Generate conversation summary
 */
export function generateConversationSummary(
  messages: ChatMessageLog[]
): string {
  const userMessages = messages.filter((m) => m.role === "user");
  const assistantMessages = messages.filter((m) => m.role === "assistant");

  const summary = `
**Conversation Summary:**
- Total Messages: ${messages.length}
- User Messages: ${userMessages.length}
- Assistant Responses: ${assistantMessages.length}
- Average Response Length: ${Math.round(
    assistantMessages.reduce((sum, m) => sum + m.content.length, 0) /
      Math.max(assistantMessages.length, 1)
  )} characters
- Topics Discussed: ${extractTopics(messages).join(", ")}
  `.trim();

  return summary;
}

/**
 * Extract topics from conversation (simplified)
 */
function extractTopics(messages: ChatMessageLog[]): string[] {
  const topics = new Set<string>();
  const keywords = [
    "docking",
    "admet",
    "toxicity",
    "pkpd",
    "analog",
    "discovery",
    "scaffold",
    "receptor",
    "binding",
    "affinity",
    "pharmacophore",
    "smiles",
    "compound",
  ];

  messages.forEach((msg) => {
    const content = msg.content.toLowerCase();
    keywords.forEach((keyword) => {
      if (content.includes(keyword)) {
        topics.add(keyword);
      }
    });
  });

  return Array.from(topics);
}
