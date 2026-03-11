import { useState } from 'react';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Checkbox } from '@/components/ui/checkbox';
import { trpc } from '@/lib/trpc';
import { Loader2, Plus, Trash2, Save } from 'lucide-react';

const COMMON_TARGETS = [
  'NMDA Receptor',
  '5HT2A Receptor',
  'Dopamine D2 Receptor',
  'Serotonin Transporter',
  'Norepinephrine Transporter',
  'GABA-A Receptor',
  'Opioid Receptor',
  'Acetylcholine Receptor',
];

const PRESET_CONFIGS: Record<string, any> = {
  'NMDA Receptor': {
    boxCenterX: '0.0',
    boxCenterY: '0.0',
    boxCenterZ: '0.0',
    boxSizeX: '25',
    boxSizeY: '25',
    boxSizeZ: '25',
    exhaustiveness: 16,
    numPoses: 9,
  },
  '5HT2A Receptor': {
    boxCenterX: '5.0',
    boxCenterY: '5.0',
    boxCenterZ: '5.0',
    boxSizeX: '20',
    boxSizeY: '20',
    boxSizeZ: '20',
    exhaustiveness: 12,
    numPoses: 9,
  },
  'Dopamine D2 Receptor': {
    boxCenterX: '0.0',
    boxCenterY: '0.0',
    boxCenterZ: '0.0',
    boxSizeX: '22',
    boxSizeY: '22',
    boxSizeZ: '22',
    exhaustiveness: 14,
    numPoses: 9,
  },
};

interface DockingParam {
  id?: number;
  name: string;
  targetName: string;
  boxCenterX: string;
  boxCenterY: string;
  boxCenterZ: string;
  boxSizeX: string;
  boxSizeY: string;
  boxSizeZ: string;
  exhaustiveness: number;
  numPoses: number;
  isDefault: boolean;
}

export function DockingParametersPanel() {
  const [selectedTarget, setSelectedTarget] = useState('');
  const [params, setParams] = useState<DockingParam>({
    name: '',
    targetName: '',
    boxCenterX: '0',
    boxCenterY: '0',
    boxCenterZ: '0',
    boxSizeX: '20',
    boxSizeY: '20',
    boxSizeZ: '20',
    exhaustiveness: 8,
    numPoses: 9,
    isDefault: false,
  });
  const [isLoading, setIsLoading] = useState(false);

  const listQuery = trpc.dockingParams.list.useQuery({ targetName: selectedTarget });
  const createMutation = trpc.dockingParams.create.useMutation();
  const deleteMutation = trpc.dockingParams.delete.useMutation();

  const handlePresetSelect = (target: string) => {
    setSelectedTarget(target);
    const preset = PRESET_CONFIGS[target];
    if (preset) {
      setParams((prev) => ({
        ...prev,
        targetName: target,
        ...preset,
      }));
    }
  };

  const handleSaveParams = async () => {
    if (!params.name || !params.targetName) {
      alert('❌ Please provide a name and select a target');
      return;
    }

    setIsLoading(true);
    try {
      await createMutation.mutateAsync({
        name: params.name,
        targetName: params.targetName,
        boxCenterX: params.boxCenterX,
        boxCenterY: params.boxCenterY,
        boxCenterZ: params.boxCenterZ,
        boxSizeX: params.boxSizeX,
        boxSizeY: params.boxSizeY,
        boxSizeZ: params.boxSizeZ,
        exhaustiveness: params.exhaustiveness,
        numPoses: params.numPoses,
        isDefault: params.isDefault,
      });

      alert(`✅ Parameters saved: ${params.name}`);
      setParams({
        name: '',
        targetName: '',
        boxCenterX: '0',
        boxCenterY: '0',
        boxCenterZ: '0',
        boxSizeX: '20',
        boxSizeY: '20',
        boxSizeZ: '20',
        exhaustiveness: 8,
        numPoses: 9,
        isDefault: false,
      });
      listQuery.refetch();
    } catch (error: any) {
      alert(`❌ Failed to save: ${error.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  const handleDelete = async (id: number) => {
    if (!confirm('Are you sure you want to delete this configuration?')) return;

    try {
      await deleteMutation.mutateAsync({ id });
      alert('✅ Configuration deleted');
      listQuery.refetch();
    } catch (error: any) {
      alert(`❌ Failed to delete: ${error.message}`);
    }
  };

  return (
    <div className="space-y-6">
      {/* Preset Configurations */}
      <Card>
        <CardHeader>
          <CardTitle>Preset Configurations</CardTitle>
          <CardDescription>Load pre-configured parameters for common targets</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 gap-2">
            {COMMON_TARGETS.map((target) => (
              <Button
                key={target}
                variant={selectedTarget === target ? 'default' : 'outline'}
                onClick={() => handlePresetSelect(target)}
                className="text-sm"
              >
                {target}
              </Button>
            ))}
          </div>
        </CardContent>
      </Card>

      {/* Parameter Editor */}
      <Card>
        <CardHeader>
          <CardTitle>Edit Docking Parameters</CardTitle>
          <CardDescription>Customize box center, size, and search parameters</CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Name */}
          <div className="space-y-2">
            <Label htmlFor="param-name">Configuration Name</Label>
            <Input
              id="param-name"
              placeholder="e.g., NMDA-High-Precision"
              value={params.name}
              onChange={(e) => setParams({ ...params, name: e.target.value })}
              disabled={isLoading}
            />
          </div>

          {/* Target Selection */}
          <div className="space-y-2">
            <Label htmlFor="param-target">Target Protein</Label>
            <Select value={params.targetName} onValueChange={(value) => setParams({ ...params, targetName: value })}>
              <SelectTrigger id="param-target">
                <SelectValue placeholder="Select target" />
              </SelectTrigger>
              <SelectContent>
                {COMMON_TARGETS.map((target) => (
                  <SelectItem key={target} value={target}>
                    {target}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Box Center */}
          <div className="grid grid-cols-3 gap-2">
            <div className="space-y-2">
              <Label htmlFor="center-x">Box Center X</Label>
              <Input
                id="center-x"
                type="number"
                step="0.1"
                value={params.boxCenterX}
                onChange={(e) => setParams({ ...params, boxCenterX: e.target.value })}
                disabled={isLoading}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="center-y">Box Center Y</Label>
              <Input
                id="center-y"
                type="number"
                step="0.1"
                value={params.boxCenterY}
                onChange={(e) => setParams({ ...params, boxCenterY: e.target.value })}
                disabled={isLoading}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="center-z">Box Center Z</Label>
              <Input
                id="center-z"
                type="number"
                step="0.1"
                value={params.boxCenterZ}
                onChange={(e) => setParams({ ...params, boxCenterZ: e.target.value })}
                disabled={isLoading}
              />
            </div>
          </div>

          {/* Box Size */}
          <div className="grid grid-cols-3 gap-2">
            <div className="space-y-2">
              <Label htmlFor="size-x">Box Size X (Å)</Label>
              <Input
                id="size-x"
                type="number"
                step="1"
                min="10"
                max="50"
                value={params.boxSizeX}
                onChange={(e) => setParams({ ...params, boxSizeX: e.target.value })}
                disabled={isLoading}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="size-y">Box Size Y (Å)</Label>
              <Input
                id="size-y"
                type="number"
                step="1"
                min="10"
                max="50"
                value={params.boxSizeY}
                onChange={(e) => setParams({ ...params, boxSizeY: e.target.value })}
                disabled={isLoading}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="size-z">Box Size Z (Å)</Label>
              <Input
                id="size-z"
                type="number"
                step="1"
                min="10"
                max="50"
                value={params.boxSizeZ}
                onChange={(e) => setParams({ ...params, boxSizeZ: e.target.value })}
                disabled={isLoading}
              />
            </div>
          </div>

          {/* Search Parameters */}
          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="exhaustiveness">Exhaustiveness (1-32)</Label>
              <Input
                id="exhaustiveness"
                type="number"
                min="1"
                max="32"
                value={params.exhaustiveness}
                onChange={(e) => setParams({ ...params, exhaustiveness: parseInt(e.target.value) || 8 })}
                disabled={isLoading}
              />
              <p className="text-xs text-muted-foreground">Higher = more thorough but slower</p>
            </div>
            <div className="space-y-2">
              <Label htmlFor="num-poses">Number of Poses (1-20)</Label>
              <Input
                id="num-poses"
                type="number"
                min="1"
                max="20"
                value={params.numPoses}
                onChange={(e) => setParams({ ...params, numPoses: parseInt(e.target.value) || 9 })}
                disabled={isLoading}
              />
              <p className="text-xs text-muted-foreground">Predicted binding poses</p>
            </div>
          </div>

          {/* Default Checkbox */}
          <div className="flex items-center space-x-2">
            <Checkbox
              id="is-default"
              checked={params.isDefault}
              onCheckedChange={(checked) => setParams({ ...params, isDefault: checked as boolean })}
              disabled={isLoading}
            />
            <Label htmlFor="is-default" className="cursor-pointer">
              Set as default for this target
            </Label>
          </div>

          {/* Save Button */}
          <Button onClick={handleSaveParams} disabled={isLoading} className="w-full">
            {isLoading && <Loader2 className="w-4 h-4 mr-2 animate-spin" />}
            {isLoading ? 'Saving...' : 'Save Configuration'}
          </Button>
        </CardContent>
      </Card>

      {/* Saved Configurations */}
      {selectedTarget && (
        <Card>
          <CardHeader>
            <CardTitle>Saved Configurations for {selectedTarget}</CardTitle>
            <CardDescription>Your custom parameter sets for this target</CardDescription>
          </CardHeader>
          <CardContent>
            {listQuery.isLoading ? (
              <div className="flex justify-center py-8">
                <Loader2 className="w-6 h-6 animate-spin text-muted-foreground" />
              </div>
            ) : listQuery.data && listQuery.data.length > 0 ? (
              <div className="space-y-2">
                {listQuery.data.map((config: any) => (
                  <div key={config.id} className="flex items-center justify-between p-3 border rounded-lg">
                    <div>
                      <p className="font-medium">{config.name}</p>
                      <p className="text-sm text-muted-foreground">
                        Exhaustiveness: {config.exhaustiveness}, Poses: {config.numPoses}
                        {config.isDefault === 1 && ' • Default'}
                      </p>
                    </div>
                    <Button
                      variant="ghost"
                      size="sm"
                      onClick={() => handleDelete(config.id)}
                      disabled={isLoading}
                    >
                      <Trash2 className="w-4 h-4" />
                    </Button>
                  </div>
                ))}
              </div>
            ) : (
              <p className="text-center text-muted-foreground py-8">No configurations saved yet</p>
            )}
          </CardContent>
        </Card>
      )}
    </div>
  );
}
