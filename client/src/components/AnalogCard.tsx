import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Link, useLocation } from "wouter";
import {
  Download,
  Beaker,
  Box,
  Copy,
  Check,
  FlaskConical,
  ChevronDown,
  ClipboardEdit,
  X,
} from "lucide-react";
import { useState } from "react";
import { MoleculeViewer3D } from "./MoleculeViewer3D";
import { trpc } from "@/lib/trpc";
import { createTRPCClient, httpBatchLink } from "@trpc/client";
import type { AppRouter } from "../../../server/routers";
import SuperJSON from "superjson";
import { toast } from "sonner";
import type { AnalogDiscovery } from "../types";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
} from "@/components/ui/dialog";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

interface AnalogCardProps {
  analog: AnalogDiscovery;
  onRunTest?: (analogId: number, testType: string) => void;
}

interface ScreeningInfo {
  dockingScore?: string;
  hergRisk?: string;
  bbbPermeability?: string;
  amesResult?: string;
  diliRisk?: string;
  notes?: string;
}

export function AnalogCard({ analog, onRunTest }: AnalogCardProps) {
  const [show3D, setShow3D] = useState(false);
  const [copied, setCopied] = useState(false);
  const [screeningDialogOpen, setScreeningDialogOpen] = useState(false);
  const [screeningInfo, setScreeningInfo] = useState<ScreeningInfo>({});
  const [, navigate] = useLocation();

  const handleExport = async (format: "smiles" | "sdf" | "pdf") => {
    try {
      const client = createTRPCClient<AppRouter>({
        links: [
          httpBatchLink({
            url: "/api/trpc",
            transformer: SuperJSON,
          }),
        ],
      });

      let data;
      if (format === "smiles") {
        data = await client.export.smiles.query({ analogId: analog.id });
      } else if (format === "sdf") {
        data = await client.export.sdf.query({ analogId: analog.id });
      } else {
        data = await client.export.pdf.query({ analogId: analog.id });
      }

      if (data) {
        const blob = new Blob([data.content], { type: data.mimeType });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = data.filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        toast.success(`Exported ${analog.compoundName} as ${format.toUpperCase()}`);
      }
    } catch (error) {
      toast.error(`Export failed: ${error}`);
    }
  };

  const handleCopySmiles = async () => {
    try {
      await navigator.clipboard.writeText(analog.smiles);
      setCopied(true);
      toast.success("SMILES copied to clipboard");
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // Fallback for environments without clipboard API
      const el = document.createElement("textarea");
      el.value = analog.smiles;
      document.body.appendChild(el);
      el.select();
      document.execCommand("copy");
      document.body.removeChild(el);
      setCopied(true);
      toast.success("SMILES copied to clipboard");
      setTimeout(() => setCopied(false), 2000);
    }
  };

  const handleQuickDock = () => {
    // Navigate to compound testing with this analog pre-selected
    // Store the analog ID in sessionStorage so CompoundTesting can pick it up
    sessionStorage.setItem("pharmasight_quick_dock_analog", String(analog.id));
    navigate("/admin/compound-testing");
    toast.info(`Opening ${analog.compoundName} in Compound Testing → Docking tab`);
  };

  const handleSaveScreeningInfo = () => {
    // Store in localStorage keyed by analog ID (persistent across sessions)
    const key = `pharmasight_screening_${analog.id}`;
    localStorage.setItem(key, JSON.stringify(screeningInfo));
    toast.success("Screening info saved");
    setScreeningDialogOpen(false);
  };

  const loadScreeningInfo = () => {
    const key = `pharmasight_screening_${analog.id}`;
    const stored = localStorage.getItem(key);
    if (stored) {
      try {
        setScreeningInfo(JSON.parse(stored));
      } catch {
        setScreeningInfo({});
      }
    }
    setScreeningDialogOpen(true);
  };

  const getConfidenceBadgeColor = (score: number) => {
    if (score >= 90) return "bg-green-100 text-green-800";
    if (score >= 80) return "bg-blue-100 text-blue-800";
    if (score >= 70) return "bg-yellow-100 text-yellow-800";
    return "bg-orange-100 text-orange-800";
  };

  const getPatentStatusColor = (status: string) => {
    if (status === "patent-free")
      return "bg-green-50 text-green-700 border-green-200";
    if (status === "patent-opportunity")
      return "bg-yellow-50 text-yellow-700 border-yellow-200";
    return "bg-red-50 text-red-700 border-red-200";
  };

  // Read stored screening info for display
  const storedScreening = (() => {
    const key = `pharmasight_screening_${analog.id}`;
    const stored = localStorage.getItem(key);
    if (!stored) return null;
    try {
      return JSON.parse(stored) as ScreeningInfo;
    } catch {
      return null;
    }
  })();

  return (
    <>
      <Card className="w-full hover:shadow-lg transition-shadow">
        <CardHeader className="pb-3">
          <div className="flex items-start justify-between gap-2">
            <div className="flex-1 min-w-0">
              <Link href={`/admin/analog/${analog.id}`}>
                <CardTitle className="text-xl mb-1 hover:text-primary cursor-pointer transition-colors truncate">
                  {analog.compoundName}
                </CardTitle>
              </Link>
              <p className="text-sm text-gray-600">Parent: {analog.parentCompound}</p>
            </div>
            <div className="flex items-center gap-1.5 flex-shrink-0">
              <Badge className={getConfidenceBadgeColor(analog.confidenceScore)}>
                {analog.confidenceScore}%
              </Badge>
              {/* Actions dropdown */}
              <DropdownMenu>
                <DropdownMenuTrigger asChild>
                  <Button variant="ghost" size="sm" className="h-7 w-7 p-0">
                    <ChevronDown className="h-4 w-4" />
                  </Button>
                </DropdownMenuTrigger>
                <DropdownMenuContent align="end" className="w-52">
                  <DropdownMenuLabel>Compound Actions</DropdownMenuLabel>
                  <DropdownMenuSeparator />
                  <DropdownMenuItem onClick={handleCopySmiles}>
                    {copied ? (
                      <Check className="h-4 w-4 mr-2 text-green-600" />
                    ) : (
                      <Copy className="h-4 w-4 mr-2" />
                    )}
                    Copy SMILES
                  </DropdownMenuItem>
                  <DropdownMenuItem onClick={() => handleExport("smiles")}>
                    <Download className="h-4 w-4 mr-2" />
                    Download SMILES file
                  </DropdownMenuItem>
                  <DropdownMenuItem onClick={() => handleExport("sdf")}>
                    <Beaker className="h-4 w-4 mr-2" />
                    Export SDF (V2000)
                  </DropdownMenuItem>
                  <DropdownMenuSeparator />
                  <DropdownMenuItem onClick={handleQuickDock}>
                    <FlaskConical className="h-4 w-4 mr-2" />
                    Quick Dock (AutoDock Vina)
                  </DropdownMenuItem>
                  <DropdownMenuSeparator />
                  <DropdownMenuItem onClick={loadScreeningInfo}>
                    <ClipboardEdit className="h-4 w-4 mr-2" />
                    Update Screening Info
                  </DropdownMenuItem>
                </DropdownMenuContent>
              </DropdownMenu>
            </div>
          </div>
        </CardHeader>

        <CardContent className="space-y-4">
          {/* SMILES Notation with copy button */}
          <div className="bg-slate-900 text-slate-100 p-3 rounded font-mono text-xs overflow-x-auto relative group">
            <div className="flex items-center justify-between mb-1">
              <p className="text-gray-400 text-xs">SMILES</p>
              <button
                onClick={handleCopySmiles}
                className="opacity-0 group-hover:opacity-100 transition-opacity text-gray-400 hover:text-white"
                title="Copy SMILES"
              >
                {copied ? (
                  <Check className="h-3.5 w-3.5 text-green-400" />
                ) : (
                  <Copy className="h-3.5 w-3.5" />
                )}
              </button>
            </div>
            <span className="break-all">{analog.smiles}</span>
          </div>

          {/* Patent Status */}
          <div
            className={`p-3 rounded border ${getPatentStatusColor(analog.patentStatus)}`}
          >
            <p className="font-semibold text-sm capitalize">
              {analog.patentStatus === "patent-free"
                ? "Patent-Free"
                : analog.patentStatus === "patent-opportunity"
                  ? "Patent Opportunity"
                  : "Patented"}
            </p>
          </div>

          {/* Scores Grid */}
          <div className="grid grid-cols-2 gap-3">
            <div className="bg-gray-50 p-3 rounded">
              <p className="text-xs text-gray-600 font-medium">Similarity</p>
              <p className="text-lg font-bold text-gray-900">
                {analog.similarityScore}%
              </p>
            </div>
            <div className="bg-gray-50 p-3 rounded">
              <p className="text-xs text-gray-600 font-medium">Safety</p>
              <p className="text-lg font-bold text-gray-900">
                {analog.safetyScore}/100
              </p>
            </div>
            <div className="bg-gray-50 p-3 rounded">
              <p className="text-xs text-gray-600 font-medium">Efficacy</p>
              <p className="text-lg font-bold text-gray-900">
                {analog.efficacyScore}/100
              </p>
            </div>
            <div className="bg-gray-50 p-3 rounded">
              <p className="text-xs text-gray-600 font-medium">Drug-Likeness</p>
              <p className="text-lg font-bold text-gray-900">
                {analog.drugLikenessScore}/100
              </p>
            </div>
          </div>

          {/* Screening Info (if saved) */}
          {storedScreening && Object.values(storedScreening).some(Boolean) && (
            <div className="border rounded-lg p-3 bg-indigo-50 border-indigo-200">
              <div className="flex items-center justify-between mb-2">
                <p className="text-xs font-semibold text-indigo-700 uppercase tracking-wide">
                  Screening Data
                </p>
                <button
                  onClick={loadScreeningInfo}
                  className="text-xs text-indigo-500 hover:text-indigo-700"
                >
                  Edit
                </button>
              </div>
              <div className="grid grid-cols-2 gap-1.5 text-xs">
                {storedScreening.dockingScore && (
                  <div>
                    <span className="text-gray-500">Docking: </span>
                    <span className="font-medium">{storedScreening.dockingScore} kcal/mol</span>
                  </div>
                )}
                {storedScreening.hergRisk && (
                  <div>
                    <span className="text-gray-500">hERG: </span>
                    <Badge
                      variant={storedScreening.hergRisk === "Low" ? "default" : "destructive"}
                      className="text-xs py-0"
                    >
                      {storedScreening.hergRisk}
                    </Badge>
                  </div>
                )}
                {storedScreening.bbbPermeability && (
                  <div>
                    <span className="text-gray-500">BBB: </span>
                    <span className="font-medium">{storedScreening.bbbPermeability}</span>
                  </div>
                )}
                {storedScreening.amesResult && (
                  <div>
                    <span className="text-gray-500">AMES: </span>
                    <span className="font-medium">{storedScreening.amesResult}</span>
                  </div>
                )}
                {storedScreening.diliRisk && (
                  <div>
                    <span className="text-gray-500">DILI: </span>
                    <Badge
                      variant={storedScreening.diliRisk === "Low" ? "default" : "destructive"}
                      className="text-xs py-0"
                    >
                      {storedScreening.diliRisk}
                    </Badge>
                  </div>
                )}
                {storedScreening.notes && (
                  <div className="col-span-2 mt-1 text-gray-600 italic">
                    {storedScreening.notes}
                  </div>
                )}
              </div>
            </div>
          )}

          {/* Market Value */}
          {analog.marketValue && (
            <div className="border-t pt-3">
              <p className="text-xs text-gray-600 font-medium mb-1">Market Value</p>
              <p className="text-lg font-bold text-blue-600">{analog.marketValue}</p>
            </div>
          )}

          {/* Therapeutic Potential */}
          {analog.therapeuticPotential && (
            <div className="border-t pt-3">
              <p className="text-xs text-gray-600 font-medium mb-1">
                Therapeutic Potential
              </p>
              <p className="text-sm text-gray-700">{analog.therapeuticPotential}</p>
            </div>
          )}

          {/* Key Differences */}
          {analog.keyDifferences && (
            <div className="border-t pt-3">
              <p className="text-xs text-gray-600 font-medium mb-1">Key Differences</p>
              <p className="text-sm text-gray-700">{analog.keyDifferences}</p>
            </div>
          )}

          {/* 3D Structure Viewer */}
          {show3D && (
            <div className="border-t pt-3">
              <MoleculeViewer3D
                smiles={analog.smiles}
                compoundName={analog.compoundName}
                width={350}
                height={250}
              />
            </div>
          )}

          {/* Action Buttons */}
          <div className="grid grid-cols-3 gap-2 pt-3 border-t">
            <Button
              size="sm"
              variant={show3D ? "default" : "outline"}
              onClick={() => setShow3D(!show3D)}
              className="col-span-3"
            >
              <Box className="w-4 h-4 mr-1" />
              {show3D ? "Hide" : "View"} 3D Structure
            </Button>
            <Button
              size="sm"
              variant="outline"
              onClick={handleCopySmiles}
              title="Copy SMILES to clipboard"
            >
              {copied ? (
                <Check className="w-4 h-4 text-green-600" />
              ) : (
                <Copy className="w-4 h-4" />
              )}
            </Button>
            <Button
              size="sm"
              variant="outline"
              onClick={() => handleExport("sdf")}
              title="Export as SDF"
            >
              <Beaker className="w-4 h-4 mr-1" />
              SDF
            </Button>
            <Button
              size="sm"
              variant="outline"
              onClick={handleQuickDock}
              title="Open in Docking"
            >
              <FlaskConical className="w-4 h-4 mr-1" />
              Dock
            </Button>
          </div>
        </CardContent>
      </Card>

      {/* Screening Info Dialog */}
      <Dialog open={screeningDialogOpen} onOpenChange={setScreeningDialogOpen}>
        <DialogContent className="max-w-md">
          <DialogHeader>
            <DialogTitle>Update Screening Info</DialogTitle>
            <DialogDescription>
              Record experimental or predicted screening results for{" "}
              <strong>{analog.compoundName}</strong>. Data is saved locally in
              your browser.
            </DialogDescription>
          </DialogHeader>
          <div className="space-y-4 py-2">
            <div className="grid grid-cols-2 gap-4">
              <div className="space-y-1.5">
                <Label htmlFor="dockingScore">Docking Score (kcal/mol)</Label>
                <Input
                  id="dockingScore"
                  placeholder="-7.5"
                  value={screeningInfo.dockingScore ?? ""}
                  onChange={(e) =>
                    setScreeningInfo((s) => ({
                      ...s,
                      dockingScore: e.target.value,
                    }))
                  }
                />
              </div>
              <div className="space-y-1.5">
                <Label htmlFor="hergRisk">hERG Risk</Label>
                <Select
                  value={screeningInfo.hergRisk ?? ""}
                  onValueChange={(v) =>
                    setScreeningInfo((s) => ({ ...s, hergRisk: v }))
                  }
                >
                  <SelectTrigger id="hergRisk">
                    <SelectValue placeholder="Select…" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="Low">Low</SelectItem>
                    <SelectItem value="Medium">Medium</SelectItem>
                    <SelectItem value="High">High</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="space-y-1.5">
                <Label htmlFor="bbb">BBB Permeability</Label>
                <Select
                  value={screeningInfo.bbbPermeability ?? ""}
                  onValueChange={(v) =>
                    setScreeningInfo((s) => ({ ...s, bbbPermeability: v }))
                  }
                >
                  <SelectTrigger id="bbb">
                    <SelectValue placeholder="Select…" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="High">High</SelectItem>
                    <SelectItem value="Medium">Medium</SelectItem>
                    <SelectItem value="Low">Low</SelectItem>
                    <SelectItem value="None">None</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="space-y-1.5">
                <Label htmlFor="ames">AMES Mutagenicity</Label>
                <Select
                  value={screeningInfo.amesResult ?? ""}
                  onValueChange={(v) =>
                    setScreeningInfo((s) => ({ ...s, amesResult: v }))
                  }
                >
                  <SelectTrigger id="ames">
                    <SelectValue placeholder="Select…" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="Negative">Negative</SelectItem>
                    <SelectItem value="Positive">Positive</SelectItem>
                    <SelectItem value="Inconclusive">Inconclusive</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="space-y-1.5">
                <Label htmlFor="dili">DILI Risk</Label>
                <Select
                  value={screeningInfo.diliRisk ?? ""}
                  onValueChange={(v) =>
                    setScreeningInfo((s) => ({ ...s, diliRisk: v }))
                  }
                >
                  <SelectTrigger id="dili">
                    <SelectValue placeholder="Select…" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="Low">Low</SelectItem>
                    <SelectItem value="Medium">Medium</SelectItem>
                    <SelectItem value="High">High</SelectItem>
                  </SelectContent>
                </Select>
              </div>
            </div>
            <div className="space-y-1.5">
              <Label htmlFor="notes">Notes</Label>
              <Input
                id="notes"
                placeholder="Any additional screening observations…"
                value={screeningInfo.notes ?? ""}
                onChange={(e) =>
                  setScreeningInfo((s) => ({ ...s, notes: e.target.value }))
                }
              />
            </div>
          </div>
          <DialogFooter>
            <Button
              variant="outline"
              onClick={() => setScreeningDialogOpen(false)}
            >
              Cancel
            </Button>
            <Button onClick={handleSaveScreeningInfo}>Save</Button>
          </DialogFooter>
        </DialogContent>
      </Dialog>
    </>
  );
}
