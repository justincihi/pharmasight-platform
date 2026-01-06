import { useState, useEffect, useRef } from "react";
import { Bell, Check, CheckCheck, X, FlaskConical, FileText, AlertTriangle, Info } from "lucide-react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from "@/components/ui/popover";
import { toast } from "sonner";
import { cn } from "@/lib/utils";

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
  
  const previousCountRef = useRef<number>(0);

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
  }, [newNotificationsData, unreadCount, toast]);

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
      <PopoverContent className="w-80 p-0" align="end">
        <div className="flex items-center justify-between border-b px-4 py-3">
          <h3 className="font-semibold">Notifications</h3>
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
        <div className="max-h-[400px] overflow-y-auto">
          {notifications.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-8 text-muted-foreground">
              <Bell className="h-8 w-8 mb-2 opacity-50" />
              <p className="text-sm">No notifications yet</p>
              <p className="text-xs">New discoveries will appear here</p>
            </div>
          ) : (
            <div className="divide-y">
              {(notifications as Notification[]).map((notification) => (
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
                      {notification.isRead === 0 && (
                        <Button
                          variant="ghost"
                          size="icon"
                          className="h-5 w-5 flex-shrink-0"
                          onClick={(e) => handleMarkAsRead(notification.id, e)}
                        >
                          <Check className="h-3 w-3" />
                        </Button>
                      )}
                    </div>
                    <p className="text-xs text-muted-foreground line-clamp-2 mt-0.5">
                      {notification.message}
                    </p>
                    <p className="text-xs text-muted-foreground mt-1">
                      {formatTimeAgo(notification.createdAt)}
                    </p>
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
        {notifications.length > 0 && (
          <div className="border-t px-4 py-2">
            <Button
              variant="ghost"
              size="sm"
              className="w-full text-xs"
              onClick={() => {
                setIsOpen(false);
                // Could navigate to a full notifications page
              }}
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
