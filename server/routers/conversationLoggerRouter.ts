/**
 * Conversation Logger Router - tRPC procedures for chatbot logging
 */

import { router, protectedProcedure } from "../_core/trpc";
import { z } from "zod";
import { promises as fs } from "fs";
import { join } from "path";

const LOGS_DIR = "/home/ubuntu/pharmasight-admin-dashboard/conversation-logs";

export const conversationLoggerRouter = router({
  /**
   * Save a conversation to markdown file
   */
  saveConversation: protectedProcedure
    .input(
      z.object({
        conversationId: z.string(),
        messages: z.array(
          z.object({
            role: z.enum(["user", "assistant"]),
            content: z.string(),
            timestamp: z.number().optional(), // Unix timestamp in ms
          })
        ),
        topic: z.string().optional(),
        metadata: z.record(z.string(), z.any()).optional(),
      })
    )
    .mutation(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Ensure logs directory exists
        try {
          await fs.mkdir(LOGS_DIR, { recursive: true });
        } catch (e) {
          // Directory might already exist
        }

        // Generate markdown content
        const timestamp = new Date().toISOString();
        const filename = `${input.conversationId}-${Date.now()}.md`;
        const filepath = join(LOGS_DIR, filename);

        let markdown = `# Conversation Log\n\n`;
        markdown += `**ID:** ${input.conversationId}\n`;
        markdown += `**Date:** ${timestamp}\n`;
        markdown += `**User:** ${ctx.user?.name || ctx.user?.openId}\n`;
        markdown += `**Topic:** ${input.topic || "General"}\n`;

        if (input.metadata) {
          markdown += `\n## Metadata\n\n`;
          markdown += `\`\`\`json\n`;
          markdown += JSON.stringify(input.metadata, null, 2);
          markdown += `\n\`\`\`\n`;
        }

        markdown += `\n## Conversation\n\n`;

        // Add messages
        for (const msg of input.messages) {
          const role = msg.role === "user" ? "👤 User" : "🤖 Assistant";
          const time = msg.timestamp
            ? new Date(msg.timestamp).toLocaleTimeString()
            : "";
          markdown += `### ${role}${time ? ` (${time})` : ""}\n\n`;
          markdown += `${msg.content}\n\n`;
        }

        // Save to file
        await fs.writeFile(filepath, markdown, "utf-8");

        return {
          success: true,
          message: "Conversation saved successfully",
          filename,
          path: filepath,
        };
      } catch (error) {
        console.error("Error saving conversation:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * List all saved conversations
   */
  listConversations: protectedProcedure
    .input(
      z.object({
        limit: z.number().default(50),
        offset: z.number().default(0),
      })
    )
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Ensure logs directory exists
        try {
          await fs.mkdir(LOGS_DIR, { recursive: true });
        } catch (e) {
          // Directory might already exist
        }

        // List files
        const files = await fs.readdir(LOGS_DIR);
        const mdFiles = files.filter((f) => f.endsWith(".md"));

        // Sort by date (newest first)
        const sortedFiles = mdFiles.sort().reverse();

        // Apply pagination
        const paginatedFiles = sortedFiles.slice(
          input.offset,
          input.offset + input.limit
        );

        // Get file stats
        const conversations = await Promise.all(
          paginatedFiles.map(async (filename) => {
            const filepath = join(LOGS_DIR, filename);
            const stats = await fs.stat(filepath);
            return {
              filename,
              createdAt: stats.birthtime,
              modifiedAt: stats.mtime,
              size: stats.size,
            };
          })
        );

        return {
          success: true,
          conversations,
          total: mdFiles.length,
          limit: input.limit,
          offset: input.offset,
        };
      } catch (error) {
        console.error("Error listing conversations:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
          conversations: [],
          total: 0,
        };
      }
    }),

  /**
   * Get a specific conversation
   */
  getConversation: protectedProcedure
    .input(z.object({ filename: z.string() }))
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Sanitize filename to prevent path traversal
        if (
          input.filename.includes("..") ||
          input.filename.includes("/") ||
          !input.filename.endsWith(".md")
        ) {
          throw new Error("Invalid filename");
        }

        const filepath = join(LOGS_DIR, input.filename);
        const content = await fs.readFile(filepath, "utf-8");

        return {
          success: true,
          filename: input.filename,
          content,
        };
      } catch (error) {
        console.error("Error reading conversation:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * Search conversations by topic or content
   */
  searchConversations: protectedProcedure
    .input(
      z.object({
        query: z.string(),
        limit: z.number().default(20),
      })
    )
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Ensure logs directory exists
        try {
          await fs.mkdir(LOGS_DIR, { recursive: true });
        } catch (e) {
          // Directory might already exist
        }

        const files = await fs.readdir(LOGS_DIR);
        const mdFiles = files.filter((f) => f.endsWith(".md"));

        const results = [];
        const searchLower = input.query.toLowerCase();

        for (const filename of mdFiles) {
          if (results.length >= input.limit) break;

          try {
            const filepath = join(LOGS_DIR, filename);
            const content = await fs.readFile(filepath, "utf-8");

            if (content.toLowerCase().includes(searchLower)) {
              results.push({
                filename,
                preview: content.substring(0, 200),
              });
            }
          } catch (e) {
            // Skip files that can't be read
          }
        }

        return {
          success: true,
          results,
          query: input.query,
          count: results.length,
        };
      } catch (error) {
        console.error("Error searching conversations:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
          results: [],
          count: 0,
        };
      }
    }),

  /**
   * Delete a conversation
   */
  deleteConversation: protectedProcedure
    .input(z.object({ filename: z.string() }))
    .mutation(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Sanitize filename
        if (
          input.filename.includes("..") ||
          input.filename.includes("/") ||
          !input.filename.endsWith(".md")
        ) {
          throw new Error("Invalid filename");
        }

        const filepath = join(LOGS_DIR, input.filename);
        await fs.unlink(filepath);

        return {
          success: true,
          message: "Conversation deleted successfully",
        };
      } catch (error) {
        console.error("Error deleting conversation:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),

  /**
   * Export conversation as markdown
   */
  exportConversation: protectedProcedure
    .input(z.object({ filename: z.string() }))
    .query(async ({ input, ctx }) => {
      try {
        if (ctx.user?.role !== "admin" && ctx.user?.openId !== process.env.OWNER_OPEN_ID) {
          throw new Error("Unauthorized: Admin access required");
        }

        // Sanitize filename
        if (
          input.filename.includes("..") ||
          input.filename.includes("/") ||
          !input.filename.endsWith(".md")
        ) {
          throw new Error("Invalid filename");
        }

        const filepath = join(LOGS_DIR, input.filename);
        const content = await fs.readFile(filepath, "utf-8");

        return {
          success: true,
          filename: input.filename,
          content,
          downloadName: `conversation-${Date.now()}.md`,
        };
      } catch (error) {
        console.error("Error exporting conversation:", error);
        return {
          success: false,
          error: error instanceof Error ? error.message : "Unknown error",
        };
      }
    }),
});
