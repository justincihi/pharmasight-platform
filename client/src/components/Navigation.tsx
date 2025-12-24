import { Link, useLocation } from "wouter";
import { Button } from "@/components/ui/button";
import { Home, LayoutDashboard, BarChart3, TestTube2, Layers, LogOut } from "lucide-react";
import { useAuth } from "@/_core/hooks/useAuth";
import { trpc } from "@/lib/trpc";
import { toast } from "sonner";

export function Navigation() {
  const [location] = useLocation();
  const { user, isAuthenticated } = useAuth();
  const logoutMutation = trpc.auth.logout.useMutation({
    onSuccess: () => {
      window.location.href = "/";
      toast.success("Logged out successfully");
    },
  });

  const navItems = [
    { path: "/", label: "Home", icon: Home },
    { path: "/admin/dashboard", label: "Dashboard", icon: LayoutDashboard, adminOnly: true },
    { path: "/analytics", label: "Analytics", icon: BarChart3, adminOnly: true },
    { path: "/testing", label: "Testing", icon: TestTube2, adminOnly: true },
    { path: "/batch", label: "Batch", icon: Layers, adminOnly: true },
  ];

  const isActive = (path: string) => {
    if (path === "/") return location === "/";
    return location.startsWith(path);
  };

  return (
    <nav className="bg-white border-b border-gray-200 sticky top-0 z-50">
      <div className="container mx-auto px-4">
        <div className="flex items-center justify-between h-16">
          {/* Logo */}
          <Link href="/">
            <a className="flex items-center gap-2 hover:opacity-80 transition-opacity">
              <div className="w-8 h-8 bg-blue-600 rounded-lg flex items-center justify-center">
                <span className="text-white font-bold text-sm">PS</span>
              </div>
              <span className="font-bold text-lg">PharmaSight™</span>
            </a>
          </Link>

          {/* Navigation Links */}
          <div className="flex items-center gap-1">
            {navItems.map((item) => {
              // Hide admin-only items if not authenticated or not admin
              if (item.adminOnly && (!isAuthenticated || user?.role !== "admin")) {
                return null;
              }

              const Icon = item.icon;
              const active = isActive(item.path);

              return (
                <Link key={item.path} href={item.path}>
                  <a>
                    <Button
                      variant={active ? "default" : "ghost"}
                      size="sm"
                      className="gap-2"
                    >
                      <Icon className="w-4 h-4" />
                      {item.label}
                    </Button>
                  </a>
                </Link>
              );
            })}

            {/* User Menu */}
            {isAuthenticated && (
              <div className="ml-4 flex items-center gap-2 pl-4 border-l border-gray-200">
                <span className="text-sm text-gray-600">
                  {user?.name || user?.email}
                </span>
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={() => logoutMutation.mutate()}
                  disabled={logoutMutation.isPending}
                >
                  <LogOut className="w-4 h-4" />
                </Button>
              </div>
            )}
          </div>
        </div>
      </div>
    </nav>
  );
}
