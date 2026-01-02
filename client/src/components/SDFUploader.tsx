import { useState, useCallback } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Upload, FileText, CheckCircle, XCircle, Loader2 } from 'lucide-react';
import { trpc } from '@/lib/trpc';

interface UploadResult {
  success: boolean;
  imported: number;
  analogs: any[];
  error?: string;
}

export default function SDFUploader() {
  const [isDragging, setIsDragging] = useState(false);
  const [isUploading, setIsUploading] = useState(false);
  const [uploadResult, setUploadResult] = useState<UploadResult | null>(null);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);

  const importMutation = trpc.analog.importFromSDF.useMutation();

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
  }, []);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);

    const files = Array.from(e.dataTransfer.files);
    const sdfFile = files.find(f => f.name.endsWith('.sdf'));

    if (sdfFile) {
      setSelectedFile(sdfFile);
      setUploadResult(null);
    } else {
      alert('❌ Please drop a valid SDF file (.sdf extension)');
    }
  }, []);

  const handleFileSelect = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file && file.name.endsWith('.sdf')) {
      setSelectedFile(file);
      setUploadResult(null);
    } else {
      alert('❌ Please select a valid SDF file (.sdf extension)');
    }
  }, []);

  const handleUpload = async () => {
    if (!selectedFile) return;

    setIsUploading(true);
    setUploadResult(null);

    try {
      // Read file content
      const text = await selectedFile.text();

      // Save to temporary location (in production, upload to server)
      // For now, we'll pass the content directly
      const blob = new Blob([text], { type: 'chemical/x-mdl-sdfile' });
      const formData = new FormData();
      formData.append('file', blob, selectedFile.name);

      // Upload file to server
      const response = await fetch('/api/upload-sdf', {
        method: 'POST',
        body: formData,
      });

      if (!response.ok) {
        throw new Error('Upload failed');
      }

      const { filePath } = await response.json();

      // Import SDF via tRPC
      const result = await importMutation.mutateAsync({ sdfPath: filePath });

      setUploadResult(result as UploadResult);
    } catch (error: any) {
      setUploadResult({
        success: false,
        imported: 0,
        analogs: [],
        error: error.message,
      });
    } finally {
      setIsUploading(false);
    }
  };

  const handleClear = () => {
    setSelectedFile(null);
    setUploadResult(null);
  };

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Upload className="w-5 h-5" />
          Import SDF Files
        </CardTitle>
        <p className="text-sm text-gray-600 mt-1">
          Upload molecular structure files to import ketamine analogs and other compounds
        </p>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Drag and Drop Zone */}
        <div
          onDragOver={handleDragOver}
          onDragLeave={handleDragLeave}
          onDrop={handleDrop}
          className={`
            border-2 border-dashed rounded-lg p-8 text-center transition-colors
            ${isDragging ? 'border-blue-500 bg-blue-50' : 'border-gray-300 bg-gray-50'}
            ${selectedFile ? 'border-green-500 bg-green-50' : ''}
          `}
        >
          {selectedFile ? (
            <div className="space-y-2">
              <FileText className="w-12 h-12 mx-auto text-green-600" />
              <p className="font-semibold text-gray-900">{selectedFile.name}</p>
              <p className="text-sm text-gray-600">
                {(selectedFile.size / 1024).toFixed(2)} KB
              </p>
            </div>
          ) : (
            <div className="space-y-2">
              <Upload className="w-12 h-12 mx-auto text-gray-400" />
              <p className="text-gray-600">
                Drag and drop an SDF file here, or click to browse
              </p>
              <input
                type="file"
                accept=".sdf"
                onChange={handleFileSelect}
                className="hidden"
                id="sdf-file-input"
              />
              <label htmlFor="sdf-file-input" className="cursor-pointer">
                <Button variant="outline" className="mt-2" type="button">
                  Browse Files
                </Button>
              </label>
            </div>
          )}
        </div>

        {/* Action Buttons */}
        {selectedFile && !uploadResult && (
          <div className="flex gap-2">
            <Button
              onClick={handleUpload}
              disabled={isUploading}
              className="flex-1"
            >
              {isUploading ? (
                <>
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                  Processing...
                </>
              ) : (
                <>
                  <Upload className="w-4 h-4 mr-2" />
                  Import Analogs
                </>
              )}
            </Button>
            <Button
              onClick={handleClear}
              variant="outline"
              disabled={isUploading}
            >
              Clear
            </Button>
          </div>
        )}

        {/* Upload Result */}
        {uploadResult && (
          <div
            className={`p-4 rounded-lg ${
              uploadResult.success
                ? 'bg-green-50 border border-green-200'
                : 'bg-red-50 border border-red-200'
            }`}
          >
            <div className="flex items-start gap-3">
              {uploadResult.success ? (
                <CheckCircle className="w-5 h-5 text-green-600 flex-shrink-0 mt-0.5" />
              ) : (
                <XCircle className="w-5 h-5 text-red-600 flex-shrink-0 mt-0.5" />
              )}
              <div className="flex-1">
                <p
                  className={`font-semibold ${
                    uploadResult.success ? 'text-green-900' : 'text-red-900'
                  }`}
                >
                  {uploadResult.success
                    ? `✅ Successfully imported ${uploadResult.imported} analog(s)`
                    : '❌ Import failed'}
                </p>
                {uploadResult.error && (
                  <p className="text-sm text-red-700 mt-1">{uploadResult.error}</p>
                )}
                {uploadResult.success && uploadResult.analogs.length > 0 && (
                  <div className="mt-2 space-y-1">
                    <p className="text-sm font-medium text-green-800">Imported compounds:</p>
                    <ul className="text-sm text-green-700 space-y-0.5">
                      {uploadResult.analogs.slice(0, 5).map((analog, idx) => (
                        <li key={idx}>
                          • {analog.compoundName} (Confidence: {analog.confidenceScore}%)
                        </li>
                      ))}
                      {uploadResult.analogs.length > 5 && (
                        <li className="text-green-600">
                          ... and {uploadResult.analogs.length - 5} more
                        </li>
                      )}
                    </ul>
                  </div>
                )}
              </div>
            </div>
            {uploadResult.success && (
              <Button
                onClick={handleClear}
                variant="outline"
                size="sm"
                className="mt-3"
              >
                Upload Another File
              </Button>
            )}
          </div>
        )}

        {/* Info */}
        <div className="text-xs text-gray-500 space-y-1">
          <p>• Supported format: SDF (Structure Data File)</p>
          <p>• Maximum file size: 10 MB</p>
          <p>• The system will automatically compute ADMET properties and molecular descriptors</p>
        </div>
      </CardContent>
    </Card>
  );
}
