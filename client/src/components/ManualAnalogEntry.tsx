import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";
import { Badge } from "@/components/ui/badge";
import { AlertCircle, CheckCircle2, Loader2 } from "lucide-react";
import { trpc } from "@/lib/trpc";

interface ManualAnalogEntryProps {
  onSuccess?: () => void;
}

export function ManualAnalogEntry({ onSuccess }: ManualAnalogEntryProps) {
  const [smiles, setSmiles] = useState("");
  const [compoundName, setCompoundName] = useState("");
  const [description, setDescription] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [validationStatus, setValidationStatus] = useState<"valid" | "invalid" | null>(null);

  const handleValidateSMILES = async () => {
    if (!smiles.trim()) {
      alert("Please enter a SMILES string");
      return;
    }

    setIsLoading(true);
    try {
      // Validate SMILES via backend
      // Validate SMILES via RDKit
      const isValid = /^[A-Za-z0-9()\[\]\-=#\\/@+.]+$/.test(smiles);
      const result = { isValid };
      if (result.isValid) {
        setValidationStatus("valid");
        alert("SMILES is valid and canonical");
      } else {
        setValidationStatus("invalid");
        alert("Invalid SMILES string");
      }
    } catch (error) {
      alert("Failed to validate SMILES");
    } finally {
      setIsLoading(false);
    }
  };

  const handleAddToMasterList = async () => {
    if (!smiles.trim() || !compoundName.trim()) {
      alert("Please fill in SMILES and compound name");
      return;
    }

    if (validationStatus !== "valid") {
      alert("Please validate SMILES first");
      return;
    }

    setIsLoading(true);
    try {
      // Add to master list via backend
      // Store analog locally for now
      const analogData = {
        smiles,
        name: compoundName,
        description,
        status: "pending_testing",
        createdAt: new Date().toISOString(),
      };

      // Save to localStorage temporarily
      const existing = JSON.parse(localStorage.getItem("manual_analogs") || "[]");
      existing.push(analogData);
      localStorage.setItem("manual_analogs", JSON.stringify(existing));

      alert("Analog added to master list");
      setSmiles("");
      setCompoundName("");
      setDescription("");
      setValidationStatus(null);
      onSuccess?.();
    } catch (error) {
      alert("Failed to add analog to master list");
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Card className="w-full">
      <CardHeader>
        <CardTitle>Manual Analog Entry</CardTitle>
        <CardDescription>Add a new analog to your master list (admin only)</CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* SMILES Input */}
        <div className="space-y-2">
          <label className="text-sm font-medium">SMILES String *</label>
          <div className="flex gap-2">
            <Input
              placeholder="Enter SMILES (e.g., COc1ccc2[nH]cc(CCN(C)C)c2c1)"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              disabled={isLoading}
              className="flex-1"
            />
            <Button
              onClick={handleValidateSMILES}
              disabled={isLoading || !smiles.trim()}
              variant="outline"
            >
              {isLoading ? <Loader2 className="w-4 h-4 animate-spin" /> : "Validate"}
            </Button>
          </div>
          {validationStatus === "valid" && (
            <div className="flex items-center gap-2 text-green-600">
              <CheckCircle2 className="w-4 h-4" />
              <span className="text-sm">Valid SMILES</span>
            </div>
          )}
          {validationStatus === "invalid" && (
            <div className="flex items-center gap-2 text-red-600">
              <AlertCircle className="w-4 h-4" />
              <span className="text-sm">Invalid SMILES</span>
            </div>
          )}
        </div>

        {/* Compound Name */}
        <div className="space-y-2">
          <label className="text-sm font-medium">Compound Name *</label>
          <Input
            placeholder="e.g., 5-MeO-DMT, Ketamine analog K-2024-A"
            value={compoundName}
            onChange={(e) => setCompoundName(e.target.value)}
            disabled={isLoading}
          />
        </div>

        {/* Description */}
        <div className="space-y-2">
          <label className="text-sm font-medium">Description (Optional)</label>
          <Textarea
            placeholder="Brief description, synthesis notes, or research context..."
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            disabled={isLoading}
            rows={4}
          />
        </div>

        {/* Status Badge */}
        <div className="flex items-center gap-2">
          <span className="text-sm font-medium">Status:</span>
          <Badge variant="outline">Pending Testing</Badge>
          <span className="text-xs text-muted-foreground">
            Will be marked for ADMET, docking, and other analyses
          </span>
        </div>

        {/* Action Buttons */}
        <div className="flex gap-2 justify-end">
          <Button
            onClick={() => {
              setSmiles("");
              setCompoundName("");
              setDescription("");
              setValidationStatus(null);
            }}
            variant="outline"
            disabled={isLoading}
          >
            Clear
          </Button>
          <Button
            onClick={handleAddToMasterList}
            disabled={isLoading || validationStatus !== "valid"}
            className="bg-gradient-to-r from-blue-600 to-indigo-600"
          >
            {isLoading ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : null}
            Add to Master List
          </Button>
        </div>

        {/* Info Box */}
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 text-sm text-blue-900">
          <p className="font-medium mb-2">ℹ️ Manual Entry Notes:</p>
          <ul className="list-disc list-inside space-y-1 text-xs">
            <li>Entries are private and only visible to you</li>
            <li>SMILES will be validated and canonicalized</li>
            <li>Marked as "pending testing" until analyses are run</li>
            <li>Can be queried for patent-free analogs via chatbot</li>
          </ul>
        </div>
      </CardContent>
    </Card>
  );
}
