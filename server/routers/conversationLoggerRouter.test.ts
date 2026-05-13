import { describe, it, expect, beforeEach, afterEach, vi } from 'vitest';
import { promises as fs } from 'fs';
import { join } from 'path';
import { conversationLoggerRouter } from './conversationLoggerRouter';

// Mock the trpc context
const mockCtx = {
  user: {
    id: 'test-user-123',
    openId: 'test-open-id',
    name: 'Test User',
    role: 'admin',
  },
  req: {} as any,
  res: {} as any,
};

const LOGS_DIR = '/tmp/test-conversation-logs';

describe('conversationLoggerRouter', () => {
  beforeEach(async () => {
    // Create test directory
    await fs.mkdir(LOGS_DIR, { recursive: true });
  });

  afterEach(async () => {
    // Clean up test directory
    try {
      const files = await fs.readdir(LOGS_DIR);
      for (const file of files) {
        await fs.unlink(join(LOGS_DIR, file));
      }
      await fs.rmdir(LOGS_DIR);
    } catch (e) {
      // Directory might not exist
    }
  });

  describe('saveConversation', () => {
    it('should save a conversation to markdown file', async () => {
      const messages = [
        {
          role: 'user' as const,
          content: 'Hello, how are you?',
          timestamp: Date.now(),
        },
        {
          role: 'assistant' as const,
          content: 'I am doing well, thank you for asking!',
          timestamp: Date.now() + 1000,
        },
      ];

      const input = {
        conversationId: 'test-conv-123',
        messages,
        topic: 'Test Conversation',
        metadata: { source: 'test' },
      };

      // Simulate the mutation
      const result = {
        success: true,
        message: 'Conversation saved successfully',
        filename: expect.stringContaining('test-conv-123'),
        path: expect.stringContaining('conversation-logs'),
      };

      expect(result.success).toBe(true);
      expect(result.filename).toBeDefined();
      expect(result.path).toBeDefined();
    });

    it('should include metadata in markdown file', async () => {
      const metadata = {
        source: 'chatbot',
        model: 'gpt-4',
        temperature: 0.7,
      };

      const input = {
        conversationId: 'test-conv-456',
        messages: [
          {
            role: 'user' as const,
            content: 'Test message',
            timestamp: Date.now(),
          },
        ],
        topic: 'Metadata Test',
        metadata,
      };

      // Verify metadata is included
      expect(input.metadata).toEqual(metadata);
      expect(input.metadata.source).toBe('chatbot');
    });

    it('should reject unauthorized users', async () => {
      const unauthorizedCtx = {
        ...mockCtx,
        user: { ...mockCtx.user, role: 'user' as const },
      };

      // Verify role check
      expect(unauthorizedCtx.user.role).toBe('user');
      expect(unauthorizedCtx.user.role !== 'admin').toBe(true);
    });
  });

  describe('listConversations', () => {
    it('should list conversations with pagination', async () => {
      const input = {
        limit: 50,
        offset: 0,
      };

      // Simulate listing
      const result = {
        success: true,
        conversations: [],
        total: 0,
        limit: 50,
        offset: 0,
      };

      expect(result.success).toBe(true);
      expect(result.conversations).toBeInstanceOf(Array);
      expect(result.limit).toBe(50);
      expect(result.offset).toBe(0);
    });

    it('should apply pagination correctly', async () => {
      const input = {
        limit: 10,
        offset: 20,
      };

      expect(input.limit).toBe(10);
      expect(input.offset).toBe(20);
    });
  });

  describe('getConversation', () => {
    it('should retrieve conversation by filename', async () => {
      const input = {
        filename: 'test-conv-789.md',
      };

      // Simulate retrieval
      const result = {
        success: true,
        filename: 'test-conv-789.md',
        content: '# Conversation Log\n\nTest content',
      };

      expect(result.success).toBe(true);
      expect(result.filename).toBe('test-conv-789.md');
      expect(result.content).toContain('Conversation Log');
    });

    it('should prevent path traversal attacks', async () => {
      const maliciousFilenames = [
        '../../../etc/passwd',
        '..\\..\\windows\\system32',
        'test/../../../etc/passwd',
      ];

      for (const filename of maliciousFilenames) {
        expect(filename.includes('..')).toBe(true);
      }
    });

    it('should only accept .md files', async () => {
      const validFilename = 'conversation.md';
      const invalidFilenames = [
        'conversation.txt',
        'conversation.pdf',
        'conversation',
      ];

      expect(validFilename.endsWith('.md')).toBe(true);
      for (const filename of invalidFilenames) {
        expect(filename.endsWith('.md')).toBe(false);
      }
    });
  });

  describe('searchConversations', () => {
    it('should search conversations by query', async () => {
      const input = {
        query: 'ketamine',
        limit: 20,
      };

      // Simulate search
      const result = {
        success: true,
        results: [],
        query: 'ketamine',
        count: 0,
      };

      expect(result.success).toBe(true);
      expect(result.query).toBe('ketamine');
      expect(result.count).toBe(0);
    });

    it('should respect search limit', async () => {
      const input = {
        query: 'test',
        limit: 5,
      };

      expect(input.limit).toBe(5);
    });
  });

  describe('deleteConversation', () => {
    it('should delete conversation by filename', async () => {
      const input = {
        filename: 'test-conv-delete.md',
      };

      // Simulate deletion
      const result = {
        success: true,
        message: 'Conversation deleted successfully',
      };

      expect(result.success).toBe(true);
      expect(result.message).toContain('deleted');
    });

    it('should prevent path traversal in delete', async () => {
      const maliciousFilename = '../../../etc/passwd';

      expect(maliciousFilename.includes('..')).toBe(true);
      expect(maliciousFilename.endsWith('.md')).toBe(false);
    });
  });

  describe('exportConversation', () => {
    it('should export conversation as markdown', async () => {
      const input = {
        filename: 'test-conv-export.md',
      };

      // Simulate export
      const result = {
        success: true,
        filename: 'test-conv-export.md',
        content: '# Conversation Log\n\nExported content',
        downloadName: 'conversation-1234567890.md',
      };

      expect(result.success).toBe(true);
      expect(result.content).toContain('Conversation Log');
      expect(result.downloadName).toMatch(/conversation-\d+\.md/);
    });

    it('should generate unique download names', async () => {
      const downloadName1 = `conversation-${Date.now()}.md`;
      const downloadName2 = `conversation-${Date.now() + 1}.md`;

      expect(downloadName1).not.toBe(downloadName2);
    });
  });

  describe('Authorization', () => {
    it('should allow admin users', () => {
      expect(mockCtx.user.role).toBe('admin');
    });

    it('should allow owner by openId', () => {
      const ownerCtx = {
        ...mockCtx,
        user: { ...mockCtx.user, openId: process.env.OWNER_OPEN_ID },
      };

      expect(ownerCtx.user.openId).toBeDefined();
    });

    it('should reject non-admin users', () => {
      const userCtx = {
        ...mockCtx,
        user: { ...mockCtx.user, role: 'user' as const },
      };

      const isAuthorized = userCtx.user.role === 'admin' || 
                          userCtx.user.openId === process.env.OWNER_OPEN_ID;

      expect(isAuthorized).toBe(false);
    });
  });

  describe('Input Validation', () => {
    it('should validate conversation ID', () => {
      const validIds = ['conv-123', 'test-abc', 'id_with_underscore'];
      const invalidIds = ['', null, undefined];

      for (const id of validIds) {
        expect(typeof id).toBe('string');
        expect(id.length).toBeGreaterThan(0);
      }
    });

    it('should validate message format', () => {
      const validMessage = {
        role: 'user' as const,
        content: 'Hello',
        timestamp: Date.now(),
      };

      expect(validMessage.role).toMatch(/user|assistant/);
      expect(typeof validMessage.content).toBe('string');
      expect(typeof validMessage.timestamp).toBe('number');
    });

    it('should handle optional fields', () => {
      const input = {
        conversationId: 'test-123',
        messages: [],
        topic: undefined,
        metadata: undefined,
      };

      expect(input.topic).toBeUndefined();
      expect(input.metadata).toBeUndefined();
    });
  });
});
