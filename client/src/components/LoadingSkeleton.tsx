import { cn } from "@/lib/utils";

interface LoadingSkeletonProps {
  className?: string;
  count?: number;
  height?: string;
}

export function LoadingSkeleton({
  className,
  count = 1,
  height = "h-4",
}: LoadingSkeletonProps) {
  return (
    <>
      {Array.from({ length: count }).map((_, index) => (
        <div
          key={index}
          className={cn(
            "glass rounded shimmer",
            height,
            className
          )}
        />
      ))}
    </>
  );
}

export function CardSkeleton() {
  return (
    <div className="glass-strong rounded-xl border border-glass-border p-6">
      <LoadingSkeleton height="h-6" className="w-1/3 mb-4" />
      <LoadingSkeleton height="h-10" className="w-1/2 mb-2" />
      <LoadingSkeleton height="h-4" className="w-full" />
    </div>
  );
}

export function TableSkeleton({ rows = 5 }: { rows?: number }) {
  return (
    <div className="space-y-3">
      {Array.from({ length: rows }).map((_, index) => (
        <div key={index} className="flex gap-4">
          <LoadingSkeleton height="h-12" className="w-1/4" />
          <LoadingSkeleton height="h-12" className="w-1/2" />
          <LoadingSkeleton height="h-12" className="w-1/4" />
        </div>
      ))}
    </div>
  );
}
