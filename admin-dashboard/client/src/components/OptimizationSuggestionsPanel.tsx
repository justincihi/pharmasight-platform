import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Lightbulb, TrendingUp, TrendingDown, ArrowRight, Sparkles } from "lucide-react";
import { cn } from "@/lib/utils";

interface OptimizationSuggestion {
  modification: string;
  category: string;
  rationale: string;
  original_smiles: string;
  optimized_smiles: string;
  property_changes: {
    mw_change: number;
    logP_change: number;
    tpsa_change: number;
  };
  priority: number;
}

interface OptimizationSuggestionsPanelProps {
  suggestions: OptimizationSuggestion[];
  onApplySuggestion?: (suggestion: OptimizationSuggestion) => void;
  onViewDetails?: (suggestion: OptimizationSuggestion) => void;
  compact?: boolean;
}

export function OptimizationSuggestionsPanel({
  suggestions,
  onApplySuggestion,
  onViewDetails,
  compact = false
}: OptimizationSuggestionsPanelProps) {
  const getCategoryColor = (category: string) => {
    const categoryMap: Record<string, string> = {
      'reduce_lipophilicity': 'bg-blue-100 text-blue-800 border-blue-200',
      'reduce_molecular_weight': 'bg-purple-100 text-purple-800 border-purple-200',
      'improve_solubility': 'bg-cyan-100 text-cyan-800 border-cyan-200',
      'reduce_hERG_risk': 'bg-red-100 text-red-800 border-red-200',
      'improve_metabolic_stability': 'bg-green-100 text-green-800 border-green-200',
      'reduce_toxicity': 'bg-orange-100 text-orange-800 border-orange-200',
    };
    return categoryMap[category] || 'bg-gray-100 text-gray-800 border-gray-200';
  };

  const formatCategory = (category: string) => {
    return category
      .split('_')
      .map(word => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');
  };

  const renderPropertyChange = (label: string, value: number, unit: string = '') => {
    if (Math.abs(value) < 0.01) return null;
    
    const isPositive = value > 0;
    const Icon = isPositive ? TrendingUp : TrendingDown;
    const colorClass = isPositive ? 'text-red-600' : 'text-green-600';
    
    return (
      <div className="flex items-center gap-1 text-sm">
        <Icon className={cn("h-3 w-3", colorClass)} />
        <span className="text-muted-foreground">{label}:</span>
        <span className={cn("font-medium", colorClass)}>
          {isPositive ? '+' : ''}{value.toFixed(2)}{unit}
        </span>
      </div>
    );
  };

  if (suggestions.length === 0) {
    return (
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Lightbulb className="h-5 w-5" />
            Optimization Suggestions
          </CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-sm text-muted-foreground text-center py-8">
            No optimization suggestions available. The molecule may already have optimal properties.
          </p>
        </CardContent>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Sparkles className="h-5 w-5" />
          AI-Driven Optimization Suggestions
        </CardTitle>
        <CardDescription>
          Structural modifications to improve drug-like properties
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div className="space-y-3">
          {suggestions.slice(0, compact ? 3 : 10).map((suggestion, index) => (
            <div
              key={index}
              className="border rounded-lg p-4 hover:border-primary/50 transition-colors"
            >
              <div className="flex items-start justify-between mb-2">
                <div className="flex-1">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="text-sm font-semibold">#{index + 1}</span>
                    <h4 className="font-semibold">{suggestion.modification}</h4>
                    <Badge className={cn("text-xs", getCategoryColor(suggestion.category))}>
                      {formatCategory(suggestion.category)}
                    </Badge>
                  </div>
                  <p className="text-sm text-muted-foreground">
                    {suggestion.rationale}
                  </p>
                </div>
                <div className="flex items-center gap-1 ml-2">
                  <div className="text-right">
                    <div className="text-xs text-muted-foreground">Priority</div>
                    <div className={cn(
                      "text-lg font-bold",
                      suggestion.priority >= 2 && "text-green-600",
                      suggestion.priority >= 1 && suggestion.priority < 2 && "text-blue-600",
                      suggestion.priority < 1 && "text-gray-600"
                    )}>
                      {suggestion.priority.toFixed(1)}
                    </div>
                  </div>
                </div>
              </div>

              <div className="grid grid-cols-3 gap-2 mt-3 p-2 bg-muted/30 rounded">
                {renderPropertyChange('MW', suggestion.property_changes.mw_change, ' Da')}
                {renderPropertyChange('LogP', suggestion.property_changes.logP_change)}
                {renderPropertyChange('TPSA', suggestion.property_changes.tpsa_change, ' Ų')}
              </div>

              {!compact && (
                <div className="flex items-center gap-2 mt-3">
                  {onApplySuggestion && (
                    <Button
                      size="sm"
                      onClick={() => onApplySuggestion(suggestion)}
                      className="flex items-center gap-1"
                    >
                      <ArrowRight className="h-3 w-3" />
                      Create Optimized Analog
                    </Button>
                  )}
                  {onViewDetails && (
                    <Button
                      size="sm"
                      variant="outline"
                      onClick={() => onViewDetails(suggestion)}
                    >
                      View Details
                    </Button>
                  )}
                </div>
              )}
            </div>
          ))}
        </div>

        {compact && suggestions.length > 3 && (
          <div className="mt-4 text-center">
            <Button variant="outline" size="sm">
              View All {suggestions.length} Suggestions
            </Button>
          </div>
        )}
      </CardContent>
    </Card>
  );
}

export default OptimizationSuggestionsPanel;
