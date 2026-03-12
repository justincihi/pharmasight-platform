import { useState } from 'react';
import { Button } from '@/components/ui/button';
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from '@/components/ui/dialog';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { trpc } from '@/lib/trpc';
import { Loader2, Upload } from 'lucide-react';
import { COMMON_TARGETS, TOAST_MESSAGES } from '@/../../shared/dockingConstants';

interface PDBUploadDialogProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  onSuccess?: () => void;
}

export function PDBUploadDialog({ open, onOpenChange, onSuccess }: PDBUploadDialogProps) {

  const [fileName, setFileName] = useState('');
  const [targetName, setTargetName] = useState('');
  const [customTarget, setCustomTarget] = useState('');
  const [description, setDescription] = useState('');
  const [fileContent, setFileContent] = useState('');
  const [isLoading, setIsLoading] = useState(false);

  const uploadMutation = trpc.pdb.upload.useMutation();

  const handleFileSelect = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setFileName(file.name);

    // Read file content
    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      setFileContent(content);
    };
    reader.readAsText(file);
  };

  const handleSubmit = async () => {
    if (!fileContent || !fileName || !targetName) {
      alert('❌ Missing Information: Please provide file, filename, and target name');
      return;
    }

    const finalTargetName = targetName === 'Custom Target' ? customTarget : targetName;
    if (!finalTargetName) {
      alert('❌ Missing Target Name: Please enter a custom target name');
      return;
    }

    setIsLoading(true);
    try {
      await uploadMutation.mutateAsync({
        fileContent,
        fileName,
        targetName: finalTargetName,
        description: description,
      });

      alert(`✅ Upload Successful: PDB file "${fileName}" uploaded for ${finalTargetName}`);

      // Reset form
      setFileName('');
      setTargetName('');
      setCustomTarget('');
      setDescription('');
      setFileContent('');
      onOpenChange(false);
      onSuccess?.();
      } catch (error: any) {
      alert(`❌ Upload Failed: ${error.message || 'Failed to upload PDB file'}`);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-w-md">
        <DialogHeader>
          <DialogTitle>Upload PDB Receptor File</DialogTitle>
          <DialogDescription>
            Upload a protein receptor structure (PDB format) for molecular docking simulations
          </DialogDescription>
        </DialogHeader>

        <div className="space-y-4">
          {/* File Upload */}
          <div className="space-y-2">
            <Label htmlFor="pdb-file">PDB File</Label>
            <div className="flex items-center gap-2">
              <Input
                id="pdb-file"
                type="file"
                accept=".pdb,.PDB"
                onChange={handleFileSelect}
                disabled={isLoading}
                className="flex-1"
              />
              <Upload className="w-4 h-4 text-muted-foreground" />
            </div>
            {fileName && (
              <p className="text-sm text-muted-foreground">Selected: {fileName}</p>
            )}
          </div>

          {/* Target Selection */}
          <div className="space-y-2">
            <Label htmlFor="target">Target Protein</Label>
            <Select value={targetName} onValueChange={setTargetName} disabled={isLoading}>
              <SelectTrigger id="target">
                <SelectValue placeholder="Select target protein" />
              </SelectTrigger>
              <SelectContent>
                {COMMON_TARGETS.map((target) => (
                  <SelectItem key={target.name} value={target.name}>
                    {target.description}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Custom Target Input */}
          {targetName === 'Custom Target' && (
            <div className="space-y-2">
              <Label htmlFor="custom-target">Custom Target Name</Label>
              <Input
                id="custom-target"
                placeholder="e.g., Sigma-1 Receptor"
                value={customTarget}
                onChange={(e) => setCustomTarget(e.target.value)}
                disabled={isLoading}
              />
            </div>
          )}

          {/* Description */}
          <div className="space-y-2">
            <Label htmlFor="description">Description (Optional)</Label>
            <Textarea
              id="description"
              placeholder="Add notes about this receptor structure..."
              value={description}
              onChange={(e) => setDescription(e.target.value)}
              disabled={isLoading}
              rows={3}
              className="resize-none"
            />
          </div>

          {/* Action Buttons */}
          <div className="flex gap-2 justify-end pt-4">
            <Button
              variant="outline"
              onClick={() => onOpenChange(false)}
              disabled={isLoading}
            >
              Cancel
            </Button>
            <Button
              onClick={handleSubmit}
              disabled={isLoading || !fileContent}
            >
              {isLoading && <Loader2 className="w-4 h-4 mr-2 animate-spin" />}
              {isLoading ? 'Uploading...' : 'Upload'}
            </Button>
          </div>
        </div>
      </DialogContent>
    </Dialog>
  );
}
