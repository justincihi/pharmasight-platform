import { AlertCircle } from 'lucide-react';
import { Badge } from '@/components/ui/badge';
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip';

interface DemoModeBadgeProps {
  source?: 'python' | 'fallback' | 'mock' | 'api';
  message?: string;
  variant?: 'default' | 'secondary' | 'outline' | 'destructive';
  showIcon?: boolean;
}

/**
 * DemoModeBadge Component
 * Displays when analysis results are using mock/fallback data instead of real Python services
 */
export default function DemoModeBadge({
  source = 'fallback',
  message,
  variant = 'secondary',
  showIcon = true,
}: DemoModeBadgeProps) {
  const getSourceInfo = (src: string) => {
    switch (src) {
      case 'python':
        return {
          label: 'Live Analysis',
          tooltip: 'Results from Python microservice',
        };
      case 'api':
        return {
          label: 'API Data',
          tooltip: 'Results from external API',
        };
      case 'fallback':
      case 'mock':
        return {
          label: 'Demo Mode',
          tooltip:
            message ||
            'Using simulated data. Python service unavailable. Results for visualization only.',
        };
      default:
        return {
          label: 'Demo Mode',
          tooltip: 'Using simulated data',
        };
    }
  };

  const info = getSourceInfo(source);

  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger asChild>
          <Badge
            variant={variant}
            className={`flex items-center gap-1 cursor-help ${
              source === 'python' || source === 'api' ? '' : 'animate-pulse'
            }`}
          >
            {showIcon && source !== 'python' && source !== 'api' && (
              <AlertCircle className="w-3 h-3" />
            )}
            {info.label}
          </Badge>
        </TooltipTrigger>
        <TooltipContent className="max-w-xs">
          <p className="text-sm">{info.tooltip}</p>
          {source !== 'python' && source !== 'api' && (
            <p className="text-xs mt-2 text-amber-200">
              ⚠️ For demonstration purposes only
            </p>
          )}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

/**
 * DemoModeWarning Component
 * Shows a prominent warning when demo mode is active
 */
export function DemoModeWarning({ source }: { source: string }) {
  if (source === 'python' || source === 'api') {
    return null;
  }

  return (
    <div className="p-3 mb-4 rounded-lg bg-amber-50 border border-amber-200 flex gap-3">
      <AlertCircle className="w-5 h-5 text-amber-600 flex-shrink-0 mt-0.5" />
      <div>
        <h4 className="font-semibold text-amber-900 text-sm">Demo Mode Active</h4>
        <p className="text-sm text-amber-800 mt-1">
          Results are simulated for visualization purposes. Connect a Python microservice for real analysis.
        </p>
      </div>
    </div>
  );
}
