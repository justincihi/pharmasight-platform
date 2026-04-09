import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Badge } from '@/components/ui/badge';
import { Download, FileJson, FileText, Loader2 } from 'lucide-react';
import { toast } from 'sonner';

interface DockingResult {
  compoundName: string;
  receptor: string;
  bindingAffinity: number;
  rmsd: number;
  numPoses: number;
  timestamp?: Date;
}

interface DockingResultsExportPanelProps {
  results: DockingResult[];
  title?: string;
}

export const DockingResultsExportPanel: React.FC<DockingResultsExportPanelProps> = ({
  results,
  title = 'Export Docking Results',
}) => {
  const [isExporting, setIsExporting] = useState(false);
  const [exportFormat, setExportFormat] = useState<'csv' | 'json' | 'pdf'>('csv');

  const handleExportCSV = async () => {
    if (results.length === 0) {
      toast.error('No results to export');
      return;
    }

    setIsExporting(true);
    try {
      const headers = [
        'Compound Name',
        'Receptor',
        'Binding Affinity (kcal/mol)',
        'RMSD (Å)',
        'Number of Poses',
        'Timestamp',
      ];

      const rows = results.map((r) => [
        r.compoundName,
        r.receptor,
        r.bindingAffinity.toFixed(2),
        r.rmsd.toFixed(2),
        r.numPoses,
        r.timestamp?.toISOString() || new Date().toISOString(),
      ]);

      const csv = [headers, ...rows.map((row) => row.map((cell) => `"${cell}"`).join(','))].join('\n');

      const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
      const link = document.createElement('a');
      const url = URL.createObjectURL(blob);

      link.setAttribute('href', url);
      link.setAttribute('download', `docking-results-${new Date().toISOString().split('T')[0]}.csv`);
      link.style.visibility = 'hidden';

      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);

      toast.success(`Exported ${results.length} results to CSV`);
    } catch (error) {
      console.error('CSV export error:', error);
      toast.error('Failed to export CSV');
    } finally {
      setIsExporting(false);
    }
  };

  const handleExportJSON = async () => {
    if (results.length === 0) {
      toast.error('No results to export');
      return;
    }

    setIsExporting(true);
    try {
      const json = JSON.stringify(
        {
          exportDate: new Date().toISOString(),
          totalResults: results.length,
          results,
        },
        null,
        2
      );

      const blob = new Blob([json], { type: 'application/json;charset=utf-8;' });
      const link = document.createElement('a');
      const url = URL.createObjectURL(blob);

      link.setAttribute('href', url);
      link.setAttribute('download', `docking-results-${new Date().toISOString().split('T')[0]}.json`);
      link.style.visibility = 'hidden';

      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);

      toast.success(`Exported ${results.length} results to JSON`);
    } catch (error) {
      console.error('JSON export error:', error);
      toast.error('Failed to export JSON');
    } finally {
      setIsExporting(false);
    }
  };

  const handleExportPDF = async () => {
    if (results.length === 0) {
      toast.error('No results to export');
      return;
    }

    setIsExporting(true);
    try {
      // Create a simple HTML table for PDF generation
      const html = `
        <html>
          <head>
            <title>Docking Results Report</title>
            <style>
              body { font-family: Arial, sans-serif; margin: 20px; }
              h1 { color: #333; }
              table { width: 100%; border-collapse: collapse; margin-top: 20px; }
              th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
              th { background-color: #f2f2f2; font-weight: bold; }
              tr:nth-child(even) { background-color: #f9f9f9; }
              .summary { margin-top: 20px; padding: 10px; background-color: #f0f0f0; border-radius: 5px; }
            </style>
          </head>
          <body>
            <h1>Docking Results Report</h1>
            <p><strong>Export Date:</strong> ${new Date().toLocaleString()}</p>
            <p><strong>Total Results:</strong> ${results.length}</p>
            
            <table>
              <thead>
                <tr>
                  <th>Compound Name</th>
                  <th>Receptor</th>
                  <th>Binding Affinity (kcal/mol)</th>
                  <th>RMSD (Å)</th>
                  <th>Poses</th>
                </tr>
              </thead>
              <tbody>
                ${results
                  .map(
                    (r) => `
                  <tr>
                    <td>${r.compoundName}</td>
                    <td>${r.receptor}</td>
                    <td>${r.bindingAffinity.toFixed(2)}</td>
                    <td>${r.rmsd.toFixed(2)}</td>
                    <td>${r.numPoses}</td>
                  </tr>
                `
                  )
                  .join('')}
              </tbody>
            </table>
            
            <div class="summary">
              <h3>Summary Statistics</h3>
              <p><strong>Average Binding Affinity:</strong> ${(results.reduce((a, b) => a + b.bindingAffinity, 0) / results.length).toFixed(2)} kcal/mol</p>
              <p><strong>Best Binding Affinity:</strong> ${Math.min(...results.map((r) => r.bindingAffinity)).toFixed(2)} kcal/mol</p>
              <p><strong>Average RMSD:</strong> ${(results.reduce((a, b) => a + b.rmsd, 0) / results.length).toFixed(2)} Å</p>
            </div>
          </body>
        </html>
      `;

      const blob = new Blob([html], { type: 'text/html;charset=utf-8;' });
      const link = document.createElement('a');
      const url = URL.createObjectURL(blob);

      link.setAttribute('href', url);
      link.setAttribute('download', `docking-results-${new Date().toISOString().split('T')[0]}.html`);
      link.style.visibility = 'hidden';

      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);

      toast.success(`Exported ${results.length} results to HTML`);
    } catch (error) {
      console.error('PDF export error:', error);
      toast.error('Failed to export PDF');
    } finally {
      setIsExporting(false);
    }
  };

  const getAffinityBadgeVariant = (affinity: number) => {
    if (affinity < -8) return 'default';
    if (affinity < -6) return 'secondary';
    return 'outline';
  };

  return (
    <Card className="w-full">
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Download className="w-5 h-5 text-green-500" />
          {title}
        </CardTitle>
        <CardDescription>Export docking results in multiple formats</CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <Tabs defaultValue="export" className="w-full">
          <TabsList className="grid w-full grid-cols-2">
            <TabsTrigger value="export">Export Options</TabsTrigger>
            <TabsTrigger value="preview">Preview ({results.length})</TabsTrigger>
          </TabsList>

          <TabsContent value="export" className="space-y-4 mt-4">
            <div className="grid grid-cols-1 gap-3">
              <Button
                onClick={handleExportCSV}
                disabled={isExporting || results.length === 0}
                className="w-full justify-start"
                variant="outline"
                size="lg"
              >
                {isExporting && exportFormat === 'csv' ? (
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                ) : (
                  <FileText className="w-4 h-4 mr-2" />
                )}
                Export as CSV
              </Button>

              <Button
                onClick={handleExportJSON}
                disabled={isExporting || results.length === 0}
                className="w-full justify-start"
                variant="outline"
                size="lg"
              >
                {isExporting && exportFormat === 'json' ? (
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                ) : (
                  <FileJson className="w-4 h-4 mr-2" />
                )}
                Export as JSON
              </Button>

              <Button
                onClick={handleExportPDF}
                disabled={isExporting || results.length === 0}
                className="w-full justify-start"
                variant="outline"
                size="lg"
              >
                {isExporting && exportFormat === 'pdf' ? (
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                ) : (
                  <Download className="w-4 h-4 mr-2" />
                )}
                Export as HTML/PDF
              </Button>
            </div>

            <div className="p-4 bg-blue-50 rounded-lg border border-blue-200">
              <p className="text-sm text-blue-800">
                <strong>Tip:</strong> CSV is best for spreadsheet analysis, JSON for data integration, and HTML for
                reports.
              </p>
            </div>
          </TabsContent>

          <TabsContent value="preview" className="space-y-4 mt-4">
            {results.length === 0 ? (
              <div className="text-center py-8 text-gray-500">
                <p>No results to preview. Run docking simulations first.</p>
              </div>
            ) : (
              <div className="space-y-2 max-h-96 overflow-y-auto">
                {results.slice(0, 10).map((result, idx) => (
                  <div key={idx} className="p-3 border rounded-lg hover:bg-gray-50">
                    <div className="flex justify-between items-start mb-2">
                      <div>
                        <p className="font-semibold text-sm">{result.compoundName}</p>
                        <p className="text-xs text-gray-600">{result.receptor}</p>
                      </div>
                      <Badge variant={getAffinityBadgeVariant(result.bindingAffinity)}>
                        {result.bindingAffinity.toFixed(2)}
                      </Badge>
                    </div>
                    <div className="flex gap-4 text-xs text-gray-600">
                      <span>RMSD: {result.rmsd.toFixed(2)} Å</span>
                      <span>Poses: {result.numPoses}</span>
                    </div>
                  </div>
                ))}
                {results.length > 10 && (
                  <p className="text-center text-sm text-gray-500 py-2">
                    ... and {results.length - 10} more results
                  </p>
                )}
              </div>
            )}
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  );
};

export default DockingResultsExportPanel;
