import { Badge } from "@/components/ui/badge";
import { Beaker, Clock } from "lucide-react";
import { cn } from "@/lib/utils";

interface SAScoreData {
  sa_score: number;
  difficulty: string;
  estimated_steps: string;
  recommendation: string;
  num_atoms?: number;
  num_rings?: number;
  num_stereocenters?: number;
  num_rotatable_bonds?: number;
}

interface SyntheticAccessibilityBadgeProps {
  saData: SAScoreData;
  showDetails?: boolean;
  className?: string;
}

export function SyntheticAccessibilityBadge({ 
  saData, 
  showDetails = false,
  className 
}: SyntheticAccessibilityBadgeProps) {
  const getDifficultyColor = (difficulty: string) => {
    switch (difficulty.toLowerCase()) {
      case 'easy':
        return 'bg-green-100 text-green-800 border-green-200';
      case 'moderate':
        return 'bg-blue-100 text-blue-800 border-blue-200';
      case 'challenging':
        return 'bg-yellow-100 text-yellow-800 border-yellow-200';
      case 'very difficult':
        return 'bg-red-100 text-red-800 border-red-200';
      default:
        return 'bg-gray-100 text-gray-800 border-gray-200';
    }
  };

  const getScoreColor = (score: number) => {
    if (score <= 3) return 'text-green-600';
    if (score <= 5) return 'text-blue-600';
    if (score <= 7) return 'text-yellow-600';
    return 'text-red-600';
  };

  if (!showDetails) {
    return (
      <div className={cn("inline-flex items-center gap-2", className)}>
        <Beaker className="h-4 w-4 text-muted-foreground" />
        <span className="text-sm font-medium">SA Score:</span>
        <span className={cn("text-sm font-bold", getScoreColor(saData.sa_score))}>
          {saData.sa_score.toFixed(1)}
        </span>
        <span className={cn(
          "text-xs px-2 py-0.5 rounded-full border",
          getDifficultyColor(saData.difficulty)
        )}>
          {saData.difficulty}
        </span>
      </div>
    );
  }

  return (
    <div className={cn("border rounded-lg p-4 space-y-3", className)}>
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Beaker className="h-5 w-5" />
          <h4 className="font-semibold">Synthetic Accessibility</h4>
        </div>
        <div className="flex items-center gap-2">
          <span className={cn("text-2xl font-bold", getScoreColor(saData.sa_score))}>
            {saData.sa_score.toFixed(1)}
          </span>
          <span className="text-sm text-muted-foreground">/ 10</span>
        </div>
      </div>

      <div className="flex items-center gap-2">
        <div className="flex-1 bg-gray-200 rounded-full h-2">
          <div
            className={cn(
              "h-2 rounded-full transition-all",
              saData.sa_score <= 3 && "bg-green-500",
              saData.sa_score > 3 && saData.sa_score <= 5 && "bg-blue-500",
              saData.sa_score > 5 && saData.sa_score <= 7 && "bg-yellow-500",
              saData.sa_score > 7 && "bg-red-500"
            )}
            style={{ width: `${(saData.sa_score / 10) * 100}%` }}
          />
        </div>
      </div>

      <div className="grid grid-cols-2 gap-3">
        <div className="space-y-1">
          <div className="flex items-center gap-1 text-xs text-muted-foreground">
            <span className="font-medium">Difficulty:</span>
          </div>
          <div className={cn(
            "text-sm px-2 py-1 rounded border inline-block",
            getDifficultyColor(saData.difficulty)
          )}>
            {saData.difficulty}
          </div>
        </div>

        <div className="space-y-1">
          <div className="flex items-center gap-1 text-xs text-muted-foreground">
            <Clock className="h-3 w-3" />
            <span className="font-medium">Est. Steps:</span>
          </div>
          <div className="text-sm font-medium">
            {saData.estimated_steps}
          </div>
        </div>
      </div>

      <p className="text-sm text-muted-foreground">
        {saData.recommendation}
      </p>

      {(saData.num_stereocenters !== undefined || saData.num_rings !== undefined) && (
        <div className="grid grid-cols-2 gap-2 pt-2 border-t text-xs">
          {saData.num_rings !== undefined && (
            <div>
              <span className="text-muted-foreground">Rings:</span>{' '}
              <span className="font-medium">{saData.num_rings}</span>
            </div>
          )}
          {saData.num_stereocenters !== undefined && (
            <div>
              <span className="text-muted-foreground">Stereocenters:</span>{' '}
              <span className="font-medium">{saData.num_stereocenters}</span>
            </div>
          )}
          {saData.num_rotatable_bonds !== undefined && (
            <div>
              <span className="text-muted-foreground">Rotatable Bonds:</span>{' '}
              <span className="font-medium">{saData.num_rotatable_bonds}</span>
            </div>
          )}
          {saData.num_atoms !== undefined && (
            <div>
              <span className="text-muted-foreground">Heavy Atoms:</span>{' '}
              <span className="font-medium">{saData.num_atoms}</span>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default SyntheticAccessibilityBadge;
