import { motion } from "framer-motion";
import { ReactNode } from "react";
import { cn } from "@/lib/utils";

interface GlassCardProps {
  children: ReactNode;
  className?: string;
  hover?: boolean;
  gradient?: boolean;
  animate?: boolean;
}

export function GlassCard({
  children,
  className,
  hover = false,
  gradient = false,
  animate = true,
}: GlassCardProps) {
  const baseClasses = gradient
    ? "glass-strong rounded-xl border border-glass-border bg-gradient-to-br from-blue-500/5 to-indigo-500/5"
    : "glass-strong rounded-xl border border-glass-border";

  const hoverClasses = hover
    ? "hover:border-blue-500/50 hover:shadow-lg transition-all cursor-pointer"
    : "";

  if (!animate) {
    return (
      <div className={cn(baseClasses, hoverClasses, className)}>
        {children}
      </div>
    );
  }

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      whileHover={hover ? { y: -4, scale: 1.02 } : undefined}
      className={cn(baseClasses, hoverClasses, className)}
    >
      {children}
    </motion.div>
  );
}

interface StatCardProps {
  title: string;
  value: string | number;
  description?: string;
  icon?: ReactNode;
  trend?: {
    value: number;
    isPositive: boolean;
  };
  gradient?: string;
}

export function StatCard({
  title,
  value,
  description,
  icon,
  trend,
  gradient = "from-blue-500 to-indigo-500",
}: StatCardProps) {
  return (
    <GlassCard hover gradient>
      <div className="p-6">
        <div className="flex items-start justify-between mb-4">
          <div>
            <p className="text-sm font-medium text-muted-foreground mb-1">
              {title}
            </p>
            <p className="text-3xl font-bold text-foreground">{value}</p>
          </div>
          {icon && (
            <div
              className={`w-12 h-12 rounded-xl bg-gradient-to-br ${gradient} flex items-center justify-center shadow-lg`}
            >
              {icon}
            </div>
          )}
        </div>
        {description && (
          <p className="text-sm text-muted-foreground">{description}</p>
        )}
        {trend && (
          <div className="flex items-center gap-1 mt-2">
            <span
              className={cn(
                "text-sm font-medium",
                trend.isPositive ? "text-green-500" : "text-red-500"
              )}
            >
              {trend.isPositive ? "+" : ""}
              {trend.value}%
            </span>
            <span className="text-xs text-muted-foreground">vs last month</span>
          </div>
        )}
      </div>
    </GlassCard>
  );
}
