import { useState, useEffect, useRef } from "react";
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

export function NotificationBell() {
  const [isOpen, setIsOpen] = useState(false);
  const [lastChecked, setLastChecked] = useState<string | null>(null);
  const [bookmarkedIds, setBookmarkedIds] = useState<Set<number>>(new Set());
  
  const previousCountRef = useRef<number>(0);
  const utils = trpc.useUtils();

  // Get unread count for badge
  const { data: unreadCount = 0, refetch: refetchCount } = trpc.notifications.getUnreadCount.useQuery(
    undefined,
    {
      refetchInterval: 10000, // Poll every 10 seconds
      enabled: true,
    }
  );

  // Get recent notifications
  const { data: notifications = [], refetch: refetchNotifications } = trpc.notifications.getRecent.useQuery(
    { limit: 20 },
    {
      enabled: isOpen,
    }
  );

  // Poll for new notifications
  const { data: newNotificationsData } = trpc.notifications.pollNew.useQuery(
    { since: lastChecked },
    {
      refetchInterval: 15000, // Poll every 15 seconds
      enabled: !!lastChecked,
    }
  );

  // Get user's bookmarks to check which notifications are bookmarked
  const { data: bookmarks = [] } = trpc.bookmarks.getAll.useQuery(
    { limit: 100 },
    {
      enabled: isOpen,
    }
  );

  // Update bookmarked IDs when bookmarks change
  useEffect(() => {
    const ids = new Set<number>();
    bookmarks.forEach((b: any) => {
      if (b.analogId) ids.add(b.analogId);
    });
    setBookmarkedIds(ids);
  }, [bookmarks]);

  // Mark as read mutation
  const markAsReadMutation = trpc.notifications.markAsRead.useMutation({
    onSuccess: () => {
      refetchCount();
      refetchNotifications();
    },
  });

  // Mark all as read mutation
  const markAllAsReadMutation = trpc.notifications.markAllAsRead.useMutation({
    onSuccess: () => {
      refetchCount();
      refetchNotifications();
      toast.success("All notifications marked as read");
    },
  });

  // Bookmark toggle mutation
  const toggleBookmarkMutation = trpc.bookmarks.toggle.useMutation({
    onSuccess: (result, variables) => {
      if (result.bookmarked) {
        setBookmarkedIds(prev => {
          const newSet = new Set(prev);
          newSet.add(variables.analogId);
          return newSet;
        });
        toast.success("Discovery saved to bookmarks");
      } else {
        setBookmarkedIds(prev => {
          const newSet = new Set(prev);
          newSet.delete(variables.analogId);
          return newSet;
        });
        toast.success("Removed from bookmarks");
      }
      utils.bookmarks.getAll.invalidate();
    },
  });

  // Initialize lastChecked on mount
  useEffect(() => {
    setLastChecked(new Date().toISOString());
  }, []);

  // Show toast for new notifications
  useEffect(() => {
    if (newNotificationsData?.notifications && newNotificationsData.notifications.length > 0) {
      const newOnes = newNotificationsData.notifications.filter(
        (n: Notification) => n.isRead === 0
      );
      
      if (newOnes.length > 0 && unreadCount > previousCountRef.current) {
        const latest = newOnes[0] as Notification;
        toast(getNotificationIcon(latest.notificationType) + " " + latest.title, {
          description: latest.message.substring(0, 100) + (latest.message.length > 100 ? "..." : ""),
          duration: 5000,
        });
      }
      
      setLastChecked(newNotificationsData.lastChecked);
    }
    previousCountRef.current = unreadCount;
  }, [newNotificationsData, unreadCount]);

  const getNotificationIcon = (type: string) => {
    switch (type) {
      case "new-discovery":
        return "🔬";
      case "high-confidence":
        return "⭐";
      case "patent-alert":
        return "📋";
      case "system":
        return "ℹ️";
      default:
        return "🔔";
    }
  };

  const getNotificationIconComponent = (type: string) => {
    switch (type) {
      case "new-discovery":
        return <FlaskConical className="h-4 w-4 text-blue-500" />;
      case "high-confidence":
        return <AlertTriangle className="h-4 w-4 text-yellow-500" />;
      case "patent-alert":
        return <FileText className="h-4 w-4 text-purple-500" />;
      case "system":
        return <Info className="h-4 w-4 text-gray-500" />;
      default:
        return <Bell className="h-4 w-4" />;
    }
  };

  const formatTimeAgo = (dateInput: Date | string) => {
    const date = typeof dateInput === 'string' ? new Date(dateInput) : dateInput;
    const now = new Date();
    const diffMs = now.getTime() - date.getTime();
    const diffMins = Math.floor(diffMs / 60000);
    const diffHours = Math.floor(diffMs / 3600000);
    const diffDays = Math.floor(diffMs / 86400000);

    if (diffMins < 1) return "Just now";
    if (diffMins < 60) return `${diffMins}m ago`;
    if (diffHours < 24) return `${diffHours}h ago`;
    if (diffDays < 7) return `${diffDays}d ago`;
    return date.toLocaleDateString();
  };

  const handleMarkAsRead = (notificationId: number, e: React.MouseEvent) => {
    e.stopPropagation();
    markAsReadMutation.mutate({ notificationId });
  };

  const handleMarkAllAsRead = () => {
    markAllAsReadMutation.mutate();
  };

  const handleToggleBookmark = (notification: Notification, e: React.MouseEvent) => {
    e.stopPropagation();
    if (notification.analogId) {
      toggleBookmarkMutation.mutate({
        analogId: notification.analogId,
        title: notification.title,
      });
    }
  };

  return (
    <Popover open={isOpen} onOpenChange={setIsOpen}>
      <PopoverTrigger asChild>
        <Button
          variant="ghost"
          size="icon"
          className="relative"
          aria-label="Notifications"
        >
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
            <div className="divide-y">
              {(notifications as Notification[]).map((notification) => {
                const isBookmarked = notification.analogId ? bookmarkedIds.has(notification.analogId) : false;
                
                return (
                  <div
                    key={notification.id}
                    className={cn(
                      "flex gap-3 px-4 py-3 hover:bg-muted/50 cursor-pointer transition-colors",
                      notification.isRead === 0 && "bg-blue-50/50 dark:bg-blue-950/20"
                    )}
                    onClick={() => {
                      if (notification.isRead === 0) {
                        markAsReadMutation.mutate({ notificationId: notification.id });
                      }
                      // Navigate to analog if applicable
                      if (notification.analogId) {
                        window.location.href = `/analog/${notification.analogId}`;
                      }
                    }}
                  >
                    <div className="flex-shrink-0 mt-1">
                      {getNotificationIconComponent(notification.notificationType)}
                    </div>
                    <div className="flex-1 min-w-0">
                      <div className="flex items-start justify-between gap-2">
                        <p className={cn(
                          "text-sm truncate",
                          notification.isRead === 0 && "font-medium"
                        )}>
                          {notification.title}
                        </p>
                        <div className="flex items-center gap-1 flex-shrink-0">
                          {/* Bookmark button */}
                          {notification.analogId && (
                            <Button
                              variant="ghost"
                              size="icon"
                              className={cn(
                                "h-6 w-6",
                                isBookmarked && "text-yellow-500"
                              )}
                              onClick={(e) => handleToggleBookmark(notification, e)}
                              title={isBookmarked ? "Remove from bookmarks" : "Save to bookmarks"}
                            >
                              {isBookmarked ? (
                                <BookmarkCheck className="h-3.5 w-3.5" />
                              ) : (
                                <Bookmark className="h-3.5 w-3.5" />
                              )}
                            </Button>
                          )}
                          {/* Mark as read button */}
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
                      <p className="text-xs text-muted-foreground line-clamp-2 mt-0.5">
                        {notification.message}
                      </p>
                      <p className="text-xs text-muted-foreground mt-1">
                        {formatTimeAgo(notification.createdAt)}
                      </p>
                    </div>
                  </div>
                );
              })}
            </div>
          )}
        </div>
        {notifications.length > 0 && (
          <div className="border-t px-4 py-2 flex gap-2">
            <Link href="/bookmarks" className="flex-1">
              <Button
                variant="outline"
                size="sm"
                className="w-full text-xs"
                onClick={() => setIsOpen(false)}
              >
                <Bookmark className="h-3 w-3 mr-1" />
                View Saved ({bookmarks.length})
              </Button>
            </Link>
            <Button
              variant="ghost"
              size="sm"
              className="flex-1 text-xs"
              onClick={() => setIsOpen(false)}
            >
              View all notifications
            </Button>
          </div>
        )}
      </PopoverContent>
    </Popover>
  );
}

export default NotificationBell;
