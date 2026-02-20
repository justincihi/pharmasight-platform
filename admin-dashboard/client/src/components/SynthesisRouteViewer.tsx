import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { ArrowDown, DollarSign, Clock, Beaker, ShoppingCart, AlertCircle } from 'lucide-react';

interface SynthesisStep {
  stepNumber: number;
  reaction: string;
  reagents: string[];
  conditions: string;
  estimatedYield: number;
  estimatedCost: number;
  difficulty: 'easy' | 'moderate' | 'hard';
  notes: string;
}

interface SynthesisRoute {
  targetCompound: string;
  smiles: string;
  totalSteps: number;
  totalEstimatedCost: number;
  totalEstimatedYield: number;
  estimatedTime: string;
  steps: SynthesisStep[];
  startingMaterials: string[];
  overallDifficulty: 'easy' | 'moderate' | 'hard';
  commercialAvailability: {
    material: string;
    available: boolean;
    supplier?: string;
    catalogPrice?: number;
  }[];
}

interface SynthesisRouteViewerProps {
  route: SynthesisRoute;
  onOptimize?: (goal: 'cost' | 'yield' | 'time') => void;
}

export default function SynthesisRouteViewer({ route, onOptimize }: SynthesisRouteViewerProps) {
  const getDifficultyColor = (difficulty: string) => {
    switch (difficulty) {
      case 'easy':
        return 'bg-green-100 text-green-800 border-green-300';
      case 'moderate':
        return 'bg-yellow-100 text-yellow-800 border-yellow-300';
      case 'hard':
        return 'bg-red-100 text-red-800 border-red-300';
      default:
        return 'bg-gray-100 text-gray-800 border-gray-300';
    }
  };

  return (
    <div className="space-y-6">
      {/* Overview Card */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center justify-between">
            <span>Synthesis Route Overview</span>
            <Badge className={getDifficultyColor(route.overallDifficulty)}>
              {route.overallDifficulty.toUpperCase()}
            </Badge>
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
            <div className="flex items-center gap-3 p-3 bg-blue-50 rounded-lg">
              <Beaker className="w-8 h-8 text-blue-600" />
              <div>
                <p className="text-sm text-gray-600">Total Steps</p>
                <p className="text-2xl font-bold text-blue-900">{route.totalSteps}</p>
              </div>
            </div>
            <div className="flex items-center gap-3 p-3 bg-green-50 rounded-lg">
              <DollarSign className="w-8 h-8 text-green-600" />
              <div>
                <p className="text-sm text-gray-600">Est. Cost</p>
                <p className="text-2xl font-bold text-green-900">${route.totalEstimatedCost}</p>
              </div>
            </div>
            <div className="flex items-center gap-3 p-3 bg-purple-50 rounded-lg">
              <AlertCircle className="w-8 h-8 text-purple-600" />
              <div>
                <p className="text-sm text-gray-600">Overall Yield</p>
                <p className="text-2xl font-bold text-purple-900">{route.totalEstimatedYield}%</p>
              </div>
            </div>
            <div className="flex items-center gap-3 p-3 bg-orange-50 rounded-lg">
              <Clock className="w-8 h-8 text-orange-600" />
              <div>
                <p className="text-sm text-gray-600">Est. Time</p>
                <p className="text-2xl font-bold text-orange-900">{route.estimatedTime}</p>
              </div>
            </div>
          </div>

          {/* Optimization Buttons */}
          {onOptimize && (
            <div className="mt-4 flex gap-2">
              <Button variant="outline" size="sm" onClick={() => onOptimize('cost')}>
                Optimize for Cost
              </Button>
              <Button variant="outline" size="sm" onClick={() => onOptimize('yield')}>
                Optimize for Yield
              </Button>
              <Button variant="outline" size="sm" onClick={() => onOptimize('time')}>
                Optimize for Time
              </Button>
            </div>
          )}
        </CardContent>
      </Card>

      {/* Starting Materials */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <ShoppingCart className="w-5 h-5" />
            Starting Materials
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-2">
            {route.commercialAvailability.map((item, idx) => (
              <div
                key={idx}
                className={`p-3 rounded-lg border ${
                  item.available
                    ? 'bg-green-50 border-green-200'
                    : 'bg-red-50 border-red-200'
                }`}
              >
                <div className="flex justify-between items-center">
                  <div>
                    <p className="font-medium">{item.material}</p>
                    {item.supplier && (
                      <p className="text-sm text-gray-600">Supplier: {item.supplier}</p>
                    )}
                  </div>
                  <div className="text-right">
                    <Badge variant={item.available ? 'default' : 'destructive'}>
                      {item.available ? 'Available' : 'Not Available'}
                    </Badge>
                    {item.catalogPrice && (
                      <p className="text-sm font-medium mt-1">${item.catalogPrice}</p>
                    )}
                  </div>
                </div>
              </div>
            ))}
          </div>
        </CardContent>
      </Card>

      {/* Synthesis Steps */}
      <div className="space-y-4">
        <h3 className="text-lg font-semibold">Synthesis Steps</h3>
        {route.steps.map((step, idx) => (
          <div key={idx}>
            <Card>
              <CardHeader>
                <div className="flex items-center justify-between">
                  <CardTitle className="text-base">
                    Step {step.stepNumber}: {step.reaction}
                  </CardTitle>
                  <Badge className={getDifficultyColor(step.difficulty)}>
                    {step.difficulty}
                  </Badge>
                </div>
              </CardHeader>
              <CardContent className="space-y-3">
                {/* Reagents */}
                <div>
                  <p className="text-sm font-medium text-gray-700">Reagents:</p>
                  <div className="flex flex-wrap gap-2 mt-1">
                    {step.reagents.map((reagent, ridx) => (
                      <Badge key={ridx} variant="outline">
                        {reagent}
                      </Badge>
                    ))}
                  </div>
                </div>

                {/* Conditions */}
                <div>
                  <p className="text-sm font-medium text-gray-700">Conditions:</p>
                  <p className="text-sm text-gray-600 mt-1">{step.conditions}</p>
                </div>

                {/* Metrics */}
                <div className="grid grid-cols-3 gap-4 pt-2 border-t">
                  <div>
                    <p className="text-xs text-gray-500">Yield</p>
                    <p className="text-lg font-semibold text-purple-600">{step.estimatedYield}%</p>
                  </div>
                  <div>
                    <p className="text-xs text-gray-500">Cost</p>
                    <p className="text-lg font-semibold text-green-600">${step.estimatedCost}</p>
                  </div>
                  <div>
                    <p className="text-xs text-gray-500">Difficulty</p>
                    <p className="text-lg font-semibold capitalize">{step.difficulty}</p>
                  </div>
                </div>

                {/* Notes */}
                {step.notes && (
                  <div className="p-2 bg-yellow-50 border border-yellow-200 rounded">
                    <p className="text-xs font-medium text-yellow-800">⚠️ Notes:</p>
                    <p className="text-sm text-yellow-700 mt-1">{step.notes}</p>
                  </div>
                )}
              </CardContent>
            </Card>

            {/* Arrow between steps */}
            {idx < route.steps.length - 1 && (
              <div className="flex justify-center py-2">
                <ArrowDown className="w-6 h-6 text-gray-400" />
              </div>
            )}
          </div>
        ))}
      </div>

      {/* Final Product */}
      <Card className="border-2 border-green-500">
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-green-700">
            <Beaker className="w-5 h-5" />
            Final Product: {route.targetCompound}
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="p-4 bg-green-50 rounded-lg">
            <p className="text-sm font-medium text-gray-700">SMILES:</p>
            <p className="text-sm font-mono text-gray-900 mt-1">{route.smiles}</p>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
