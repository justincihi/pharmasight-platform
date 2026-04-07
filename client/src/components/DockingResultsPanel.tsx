import React, { useState } from 'react';
import { Card } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { ChevronDown, ChevronUp, Zap } from 'lucide-react';

interface DockingPose {
  mode: number;
  affinity: number;
  rmsd: number;
}

interface DockingResults {
  success: boolean;
  binding_affinity: number;
  poses: DockingPose[];
  num_poses: number;
  error?: string;
}

interface DockingResultsPanelProps {
  results: DockingResults;
  compoundName?: string;
  receptorName?: string;
  onVisualize?: (pose: DockingPose) => void;
}

export const DockingResultsPanel: React.FC<DockingResultsPanelProps> = ({
  results,
  compoundName = 'Compound',
  receptorName = 'Receptor',
  onVisualize,
}) => {
  const [expandedPose, setExpandedPose] = useState<number | null>(0);
  const [sortBy, setSortBy] = useState<'affinity' | 'rmsd'>('affinity');

  if (!results.success) {
    return (
      <Card className="p-4 bg-red-50 border-red-200">
        <div className="flex items-center gap-2 text-red-800">
          <Zap className="w-5 h-5" />
          <div>
            <p className="font-semibold">Docking Failed</p>
            <p className="text-sm">{results.error || 'Unknown error'}</p>
          </div>
        </div>
      </Card>
    );
  }

  const sortedPoses = [...results.poses].sort((a, b) => {
    if (sortBy === 'affinity') {
      return a.affinity - b.affinity; // Lower (more negative) is better
    }
    return a.rmsd - b.rmsd; // Lower RMSD is better
  });

  const bestPose = sortedPoses[0];

  return (
    <Card className="p-6">
      <div className="mb-6">
        <h2 className="text-2xl font-bold mb-2">Docking Results</h2>
        <div className="grid grid-cols-3 gap-4">
          <div className="bg-blue-50 p-3 rounded-lg">
            <p className="text-sm text-gray-600">Compound</p>
            <p className="font-semibold text-blue-900">{compoundName}</p>
          </div>
          <div className="bg-green-50 p-3 rounded-lg">
            <p className="text-sm text-gray-600">Receptor</p>
            <p className="font-semibold text-green-900">{receptorName}</p>
          </div>
          <div className="bg-purple-50 p-3 rounded-lg">
            <p className="text-sm text-gray-600">Poses Generated</p>
            <p className="font-semibold text-purple-900">{results.num_poses}</p>
          </div>
        </div>
      </div>

      <div className="mb-6 p-4 bg-gradient-to-r from-amber-50 to-orange-50 rounded-lg border border-amber-200">
        <p className="text-sm text-gray-600 mb-1">Best Binding Affinity</p>
        <div className="flex items-baseline gap-2">
          <span className="text-3xl font-bold text-amber-900">{bestPose.affinity.toFixed(2)}</span>
          <span className="text-lg text-amber-700">kcal/mol</span>
        </div>
        <p className="text-xs text-amber-600 mt-2">
          {bestPose.affinity < -7 && '✓ Strong binding'}
          {bestPose.affinity >= -7 && bestPose.affinity < -5 && '~ Moderate binding'}
          {bestPose.affinity >= -5 && '⚠ Weak binding'}
        </p>
      </div>

      <Tabs defaultValue="poses" className="w-full">
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="poses">Poses ({results.num_poses})</TabsTrigger>
          <TabsTrigger value="analysis">Analysis</TabsTrigger>
        </TabsList>

        <TabsContent value="poses" className="space-y-3 mt-4">
          <div className="flex gap-2 mb-4">
            <Button
              variant={sortBy === 'affinity' ? 'default' : 'outline'}
              size="sm"
              onClick={() => setSortBy('affinity')}
            >
              Sort by Affinity
            </Button>
            <Button
              variant={sortBy === 'rmsd' ? 'default' : 'outline'}
              size="sm"
              onClick={() => setSortBy('rmsd')}
            >
              Sort by RMSD
            </Button>
          </div>

          {sortedPoses.map((pose, idx) => (
            <div
              key={pose.mode}
              className="border rounded-lg overflow-hidden hover:shadow-md transition-shadow"
            >
              <button
                onClick={() => setExpandedPose(expandedPose === pose.mode ? null : pose.mode)}
                className="w-full p-4 bg-gray-50 hover:bg-gray-100 flex items-center justify-between"
              >
                <div className="flex items-center gap-4 text-left">
                  <div className="font-mono font-bold text-lg w-12">#{pose.mode}</div>
                  <div>
                    <p className="font-semibold">
                      {pose.affinity.toFixed(2)} kcal/mol
                    </p>
                    <p className="text-sm text-gray-600">
                      RMSD: {pose.rmsd.toFixed(2)} Å
                    </p>
                  </div>
                </div>
                {expandedPose === pose.mode ? (
                  <ChevronUp className="w-5 h-5" />
                ) : (
                  <ChevronDown className="w-5 h-5" />
                )}
              </button>

              {expandedPose === pose.mode && (
                <div className="p-4 bg-white border-t">
                  <div className="grid grid-cols-2 gap-4 mb-4">
                    <div>
                      <p className="text-sm text-gray-600">Binding Affinity</p>
                      <p className="text-xl font-bold text-blue-600">{pose.affinity.toFixed(2)}</p>
                      <p className="text-xs text-gray-500 mt-1">kcal/mol</p>
                    </div>
                    <div>
                      <p className="text-sm text-gray-600">RMSD (Lower Bound)</p>
                      <p className="text-xl font-bold text-green-600">{pose.rmsd.toFixed(2)}</p>
                      <p className="text-xs text-gray-500 mt-1">Ångströms</p>
                    </div>
                  </div>

                  {onVisualize && (
                    <Button
                      onClick={() => onVisualize(pose)}
                      className="w-full mt-4"
                      variant="default"
                    >
                      View 3D Structure
                    </Button>
                  )}
                </div>
              )}
            </div>
          ))}
        </TabsContent>

        <TabsContent value="analysis" className="mt-4">
          <div className="space-y-4">
            <div className="p-4 bg-blue-50 rounded-lg border border-blue-200">
              <h4 className="font-semibold text-blue-900 mb-2">Binding Affinity Interpretation</h4>
              <ul className="text-sm text-blue-800 space-y-1">
                <li>• &lt; -7 kcal/mol: Strong binding (likely active)</li>
                <li>• -7 to -5 kcal/mol: Moderate binding (possible activity)</li>
                <li>• &gt; -5 kcal/mol: Weak binding (unlikely to be active)</li>
              </ul>
            </div>

            <div className="p-4 bg-green-50 rounded-lg border border-green-200">
              <h4 className="font-semibold text-green-900 mb-2">RMSD Interpretation</h4>
              <ul className="text-sm text-green-800 space-y-1">
                <li>• &lt; 2 Å: Excellent pose reproducibility</li>
                <li>• 2-3 Å: Good pose reproducibility</li>
                <li>• &gt; 3 Å: Poses may differ significantly</li>
              </ul>
            </div>

            <div className="p-4 bg-purple-50 rounded-lg border border-purple-200">
              <h4 className="font-semibold text-purple-900 mb-2">Best Pose Summary</h4>
              <p className="text-sm text-purple-800">
                Pose #{bestPose.mode} shows the best binding affinity at{' '}
                <span className="font-bold">{bestPose.affinity.toFixed(2)} kcal/mol</span> with an
                RMSD of <span className="font-bold">{bestPose.rmsd.toFixed(2)} Å</span>.
              </p>
            </div>
          </div>
        </TabsContent>
      </Tabs>
    </Card>
  );
};

export default DockingResultsPanel;
