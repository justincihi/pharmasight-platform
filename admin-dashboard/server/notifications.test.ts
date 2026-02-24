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

describe("Real-Time Notification System", () => {
  let adminCaller: any;
  let userCaller: any;

  beforeAll(() => {
    adminCaller = appRouter.createCaller(createAdminContext());
    userCaller = appRouter.createCaller(createUserContext());
  });

  describe("Admin Access", () => {
    it("should allow admin to get recent notifications", async () => {
      const result = await adminCaller.notifications.getRecent({ limit: 10 });
      expect(Array.isArray(result)).toBe(true);
    });

    it("should allow admin to get unread count", async () => {
      const result = await adminCaller.notifications.getUnreadCount();
      expect(typeof result).toBe("number");
      expect(result).toBeGreaterThanOrEqual(0);
    });

    it("should allow admin to poll for new notifications", async () => {
      const result = await adminCaller.notifications.pollNew({ since: null });
      expect(result).toHaveProperty("notifications");
      expect(result).toHaveProperty("lastChecked");
      expect(Array.isArray(result.notifications)).toBe(true);
    });

    it("should allow admin to mark all as read", async () => {
      const result = await adminCaller.notifications.markAllAsRead();
      expect(result.success).toBe(true);
    });
  });

  describe("Access Control", () => {
    it("should deny non-admin access to notifications", async () => {
      await expect(
        userCaller.notifications.getRecent({ limit: 10 })
      ).rejects.toThrow("Unauthorized");
    });

    it("should deny non-admin access to unread count", async () => {
      await expect(
        userCaller.notifications.getUnreadCount()
      ).rejects.toThrow("Unauthorized");
    });

    it("should deny non-admin access to poll new", async () => {
      await expect(
        userCaller.notifications.pollNew({ since: null })
      ).rejects.toThrow("Unauthorized");
    });

    it("should deny non-admin access to mark all as read", async () => {
      await expect(
        userCaller.notifications.markAllAsRead()
      ).rejects.toThrow("Unauthorized");
    });
  });

  describe("Database Functions", () => {
    it("should export getUnreadNotificationCount function", async () => {
      const { getUnreadNotificationCount } = await import("./db");
      expect(getUnreadNotificationCount).toBeDefined();
      expect(typeof getUnreadNotificationCount).toBe("function");
    });

    it("should export markNotificationAsRead function", async () => {
      const { markNotificationAsRead } = await import("./db");
      expect(markNotificationAsRead).toBeDefined();
      expect(typeof markNotificationAsRead).toBe("function");
    });

    it("should export markAllNotificationsAsRead function", async () => {
      const { markAllNotificationsAsRead } = await import("./db");
      expect(markAllNotificationsAsRead).toBeDefined();
      expect(typeof markAllNotificationsAsRead).toBe("function");
    });

    it("should export getNewNotifications function", async () => {
      const { getNewNotifications } = await import("./db");
      expect(getNewNotifications).toBeDefined();
      expect(typeof getNewNotifications).toBe("function");
    });
  });
});
