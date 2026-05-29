import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { Bell, Check, CheckCheck, FlaskConical, FileText, AlertTriangle, Info, Bookmark, BookmarkCheck } from "lucide-react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from "@/components/ui/popover";
import { toast } from "sonner";
import { cn } from "@/lib/utils";
import { Link } from "wouter";

interface Notification {
  id: number;
  title: string;
  message: string;
  notificationType: "new-discovery" | "high-confidence" | "patent-alert" | "system";
  isRead: number;
  createdAt: Date | string;
  analogId?: number | null;
}

// Stable mount timestamp — created once per module load, never changes
const MOUNT_SINCE = new Date().toISOString();

export function NotificationBell() {
  const [isOpen, setIsOpen] = useState(false);
  const [bookmarkedIds, setBookmarkedIds] = useState<Set<number>>(new Set());

  // Track previous unread count to detect new arrivals (ref = no re-render)
  const previousCountRef = useRef<number>(0);
  // Track whether we've already toasted for the current batch (ref = no re-render)
  const toastedRef = useRef<boolean>(false);

  const utils = trpc.useUtils();

  // Unread badge — polls every 30s, no input object (stable)
  const { data: unreadCount = 0 } = trpc.notifications.getUnreadCount.useQuery(
    undefined,
    { refetchInterval: 30000, staleTime: 15000 }
  );

  // Recent list — only fetched when popover is open
  const { data: notifications = [] } = trpc.notifications.getRecent.useQuery(
    { limit: 20 },
    { enabled: isOpen, staleTime: 5000 }
  );

  // Poll for new notifications since mount — STABLE input, never changes
  const { data: newNotificationsData } = trpc.notifications.pollNew.useQuery(
    { since: MOUNT_SINCE },
    { refetchInterval: 30000, staleTime: 15000 }
  );

  // Bookmarks — only fetched when popover is open
  const { data: bookmarks = [] } = trpc.bookmarks.getAll.useQuery(
    { limit: 100 },
    { enabled: isOpen, staleTime: 10000 }
  );

  // Sync bookmarked IDs
  useEffect(() => {
    const ids = new Set<number>();
    (bookmarks as any[]).forEach((b) => {
      if (b.analogId) ids.add(b.analogId);
    });
    setBookmarkedIds(ids);
  }, [bookmarks]);

  // Show toast only when unread count genuinely increases
  useEffect(() => {
    if (unreadCount > previousCountRef.current && !toastedRef.current) {
      toastedRef.current = true;
      const newOnes = (newNotificationsData?.notifications ?? []) as Notification[];
      const unread = newOnes.filter((n) => n.isRead === 0);
      if (unread.length > 0) {
        const latest = unread[0];
        toast(getNotificationIcon(latest.notificationType) + " " + latest.title, {
          description: latest.message.substring(0, 100) + (latest.message.length > 100 ? "..." : ""),
          duration: 5000,
        });
      }
    }
    if (unreadCount !== previousCountRef.current) {
      previousCountRef.current = unreadCount;
      toastedRef.current = false;
    }
    // Only re-run when unreadCount changes — newNotificationsData is intentionally excluded
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [unreadCount]);

  const markAsReadMutation = trpc.notifications.markAsRead.useMutation({
    onSuccess: () => {
      utils.notifications.getUnreadCount.invalidate();
      utils.notifications.getRecent.invalidate();
    },
  });

  const markAllAsReadMutation = trpc.notifications.markAllAsRead.useMutation({
    onSuccess: () => {
      utils.notifications.getUnreadCount.invalidate();
      utils.notifications.getRecent.invalidate();
      toast.success("All notifications marked as read");
    },
  });

  const toggleBookmarkMutation = trpc.bookmarks.toggle.useMutation({
    onSuccess: (result, variables) => {
      if (result.bookmarked) {
        setBookmarkedIds((prev) => { const s = new Set(prev); s.add(variables.analogId); return s; });
        toast.success("Discovery saved to bookmarks");
      } else {
        setBookmarkedIds((prev) => { const s = new Set(prev); s.delete(variables.analogId); return s; });
        toast.success("Removed from bookmarks");
      }
      utils.bookmarks.getAll.invalidate();
    },
  });

  const getNotificationIcon = useCallback((type: string) => {
    switch (type) {
      case "new-discovery": return "🔬";
      case "high-confidence": return "⭐";
      case "patent-alert": return "📋";
      case "system": return "ℹ️";
      default: return "🔔";
    }
  }, []);

  const getNotificationIconComponent = useCallback((type: string) => {
    switch (type) {
      case "new-discovery": return <FlaskConical className="h-4 w-4 text-blue-500" />;
      case "high-confidence": return <AlertTriangle className="h-4 w-4 text-yellow-500" />;
      case "patent-alert": return <FileText className="h-4 w-4 text-purple-500" />;
      case "system": return <Info className="h-4 w-4 text-gray-500" />;
      default: return <Bell className="h-4 w-4" />;
    }
  }, []);

  const formatTimeAgo = useCallback((dateInput: Date | string) => {
    const date = typeof dateInput === "string" ? new Date(dateInput) : dateInput;
    const diffMs = Date.now() - date.getTime();
    const diffMins = Math.floor(diffMs / 60000);
    const diffHours = Math.floor(diffMs / 3600000);
    const diffDays = Math.floor(diffMs / 86400000);
    if (diffMins < 1) return "Just now";
    if (diffMins < 60) return `${diffMins}m ago`;
    if (diffHours < 24) return `${diffHours}h ago`;
    if (diffDays < 7) return `${diffDays}d ago`;
    return date.toLocaleDateString();
  }, []);

  const handleMarkAsRead = useCallback((notificationId: number, e: React.MouseEvent) => {
    e.stopPropagation();
    markAsReadMutation.mutate({ notificationId });
  }, [markAsReadMutation]);

  const handleMarkAllAsRead = useCallback(() => {
    markAllAsReadMutation.mutate();
  }, [markAllAsReadMutation]);

  const handleToggleBookmark = useCallback((notification: Notification, e: React.MouseEvent) => {
    e.stopPropagation();
    if (notification.analogId) {
      toggleBookmarkMutation.mutate({ analogId: notification.analogId, title: notification.title });
    }
  }, [toggleBookmarkMutation]);

  const notificationItems = useMemo(() => {
    return (notifications as Notification[]).map((notification) => {
      const isBookmarked = notification.analogId ? bookmarkedIds.has(notification.analogId) : false;
      return (
        <div
          key={notification.id}
          className={cn(
            "flex gap-3 px-4 py-3 hover:bg-muted/50 cursor-pointer transition-colors",
            notification.isRead === 0 && "bg-blue-50/50 dark:bg-blue-950/20"
          )}
          onClick={() => {
            if (notification.isRead === 0) markAsReadMutation.mutate({ notificationId: notification.id });
            if (notification.analogId) window.location.href = `/analog/${notification.analogId}`;
          }}
        >
          <div className="flex-shrink-0 mt-1">
            {getNotificationIconComponent(notification.notificationType)}
          </div>
          <div className="flex-1 min-w-0">
            <div className="flex items-start justify-between gap-2">
              <p className={cn("text-sm truncate", notification.isRead === 0 && "font-medium")}>
                {notification.title}
              </p>
              <div className="flex items-center gap-1 flex-shrink-0">
                {notification.analogId && (
                  <Button
                    variant="ghost"
                    size="icon"
                    className={cn("h-6 w-6", isBookmarked && "text-yellow-500")}
                    onClick={(e) => handleToggleBookmark(notification, e)}
                    title={isBookmarked ? "Remove from bookmarks" : "Save to bookmarks"}
                  >
                    {isBookmarked ? <BookmarkCheck className="h-3.5 w-3.5" /> : <Bookmark className="h-3.5 w-3.5" />}
                  </Button>
                )}
                {notification.isRead === 0 && (
                  <Button
                    variant="ghost"
                    size="icon"
                    className="h-6 w-6"
                    onClick={(e) => handleMarkAsRead(notification.id, e)}
                    title="Mark as read"
                  >
                    <Check className="h-3.5 w-3.5" />
                  </Button>
                )}
              </div>
            </div>
            <p className="text-xs text-muted-foreground line-clamp-2 mt-0.5">{notification.message}</p>
            <p className="text-xs text-muted-foreground mt-1">{formatTimeAgo(notification.createdAt)}</p>
          </div>
        </div>
      );
    });
  }, [notifications, bookmarkedIds, handleMarkAsRead, handleToggleBookmark, getNotificationIconComponent, formatTimeAgo, markAsReadMutation]);

  return (
    <Popover open={isOpen} onOpenChange={setIsOpen}>
      <PopoverTrigger asChild>
        <Button variant="ghost" size="icon" className="relative" aria-label="Notifications">
          <Bell className="h-5 w-5" />
          {unreadCount > 0 && (
            <span className="absolute -top-1 -right-1 flex h-5 w-5 items-center justify-center rounded-full bg-red-500 text-[10px] font-bold text-white">
              {unreadCount > 99 ? "99+" : unreadCount}
            </span>
          )}
        </Button>
      </PopoverTrigger>
      <PopoverContent className="w-96 p-0" align="end">
        <div className="flex items-center justify-between border-b px-4 py-3">
          <h3 className="font-semibold">Notifications</h3>
          <div className="flex items-center gap-2">
            <Link href="/bookmarks">
              <Button
                variant="ghost"
                size="sm"
                className="h-auto py-1 px-2 text-xs text-muted-foreground hover:text-foreground"
                onClick={() => setIsOpen(false)}
              >
                <Bookmark className="h-3 w-3 mr-1" />
                Saved
              </Button>
            </Link>
            {unreadCount > 0 && (
              <Button
                variant="ghost"
                size="sm"
                className="h-auto py-1 px-2 text-xs text-muted-foreground hover:text-foreground"
                onClick={handleMarkAllAsRead}
              >
                <CheckCheck className="h-3 w-3 mr-1" />
                Mark all read
              </Button>
            )}
          </div>
        </div>
        <div className="max-h-[400px] overflow-y-auto">
          {notifications.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-8 text-muted-foreground">
              <Bell className="h-8 w-8 mb-2 opacity-50" />
              <p className="text-sm">No notifications yet</p>
              <p className="text-xs">New discoveries will appear here</p>
            </div>
          ) : (
            <div className="divide-y">{notificationItems}</div>
          )}
        </div>
        {notifications.length > 0 && (
          <div className="border-t px-4 py-2 flex gap-2">
            <Link href="/bookmarks" className="flex-1">
              <Button variant="outline" size="sm" className="w-full text-xs" onClick={() => setIsOpen(false)}>
                <Bookmark className="h-3 w-3 mr-1" />
                View Saved ({bookmarks.length})
              </Button>
            </Link>
            <Button variant="ghost" size="sm" className="flex-1 text-xs" onClick={() => setIsOpen(false)}>
              View all notifications
            </Button>
          </div>
        )}
      </PopoverContent>
    </Popover>
  );
}

export default NotificationBell;
