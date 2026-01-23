import { describe, it, expect, beforeAll } from "vitest";
import { appRouter } from "./routers";
import type { TrpcContext } from "./_core/context";

// Mock admin user context
const createAdminContext = (): TrpcContext => ({
  user: {
    id: 1,
    openId: "test-admin",
    email: "admin@test.com",
    name: "Test Admin",
    loginMethod: "test",
    role: "admin",
    createdAt: new Date(),
    updatedAt: new Date(),
    lastSignedIn: new Date(),
  },
  req: {
    protocol: "https",
    headers: {},
  } as any,
  res: {
    clearCookie: () => {},
  } as any,
});

// Mock regular user context
const createUserContext = (): TrpcContext => ({
  user: {
    id: 2,
    openId: "test-user",
    email: "user@test.com",
    name: "Test User",
    loginMethod: "test",
    role: "user",
    createdAt: new Date(),
    updatedAt: new Date(),
    lastSignedIn: new Date(),
  },
  req: {
    protocol: "https",
    headers: {},
  } as any,
  res: {
    clearCookie: () => {},
  } as any,
});

describe("Bookmark System", () => {
  let adminCaller: any;
  let userCaller: any;

  beforeAll(() => {
    adminCaller = appRouter.createCaller(createAdminContext());
    userCaller = appRouter.createCaller(createUserContext());
  });

  describe("Bookmark CRUD Operations", () => {
    it("should allow users to get their bookmarks", async () => {
      const result = await adminCaller.bookmarks.getAll({ limit: 50 });
      expect(Array.isArray(result)).toBe(true);
    });

    it("should allow users to create a bookmark", async () => {
      const result = await adminCaller.bookmarks.create({
        analogId: 1,
        title: "Test Bookmark",
        notes: "Test notes",
        category: "review-later",
      });
      expect(result.success).toBe(true);
      expect(result.bookmarkId).toBeDefined();
    });

    it("should allow users to check if an analog is bookmarked", async () => {
      const result = await adminCaller.bookmarks.isBookmarked({ analogId: 1 });
      expect(typeof result).toBe("boolean");
    });

    it("should allow users to toggle bookmark status", async () => {
      const result = await adminCaller.bookmarks.toggle({
        analogId: 999,
        title: "Toggle Test",
      });
      expect(result).toHaveProperty("bookmarked");
      expect(typeof result.bookmarked).toBe("boolean");
    });
  });

  describe("Regular User Access", () => {
    it("should allow regular users to get their bookmarks", async () => {
      const result = await userCaller.bookmarks.getAll({ limit: 50 });
      expect(Array.isArray(result)).toBe(true);
    });

    it("should allow regular users to create bookmarks", async () => {
      const result = await userCaller.bookmarks.create({
        analogId: 2,
        title: "User Bookmark",
        category: "promising",
      });
      expect(result.success).toBe(true);
    });
  });

  describe("Database Functions", () => {
    it("should export createBookmark function", async () => {
      const { createBookmark } = await import("./db");
      expect(createBookmark).toBeDefined();
      expect(typeof createBookmark).toBe("function");
    });

    it("should export getUserBookmarks function", async () => {
      const { getUserBookmarks } = await import("./db");
      expect(getUserBookmarks).toBeDefined();
      expect(typeof getUserBookmarks).toBe("function");
    });

    it("should export getBookmarkById function", async () => {
      const { getBookmarkById } = await import("./db");
      expect(getBookmarkById).toBeDefined();
      expect(typeof getBookmarkById).toBe("function");
    });

    it("should export updateBookmark function", async () => {
      const { updateBookmark } = await import("./db");
      expect(updateBookmark).toBeDefined();
      expect(typeof updateBookmark).toBe("function");
    });

    it("should export deleteBookmark function", async () => {
      const { deleteBookmark } = await import("./db");
      expect(deleteBookmark).toBeDefined();
      expect(typeof deleteBookmark).toBe("function");
    });

    it("should export isAnalogBookmarked function", async () => {
      const { isAnalogBookmarked } = await import("./db");
      expect(isAnalogBookmarked).toBeDefined();
      expect(typeof isAnalogBookmarked).toBe("function");
    });

    it("should export getBookmarkByAnalogId function", async () => {
      const { getBookmarkByAnalogId } = await import("./db");
      expect(getBookmarkByAnalogId).toBeDefined();
      expect(typeof getBookmarkByAnalogId).toBe("function");
    });
  });

  describe("Bookmark Categories", () => {
    it("should accept high-priority category", async () => {
      const result = await adminCaller.bookmarks.create({
        analogId: 100,
        title: "High Priority Test",
        category: "high-priority",
      });
      expect(result.success).toBe(true);
    });

    it("should accept promising category", async () => {
      const result = await adminCaller.bookmarks.create({
        analogId: 101,
        title: "Promising Test",
        category: "promising",
      });
      expect(result.success).toBe(true);
    });

    it("should accept archived category", async () => {
      const result = await adminCaller.bookmarks.create({
        analogId: 102,
        title: "Archived Test",
        category: "archived",
      });
      expect(result.success).toBe(true);
    });
  });
});
