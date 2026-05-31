import { useState } from "react";
import { Button } from "@/components/ui/button";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import { Download, FileText, Table2, Loader2 } from "lucide-react";
import { toast } from "sonner";

// ── Types ──────────────────────────────────────────────────────────────────────

interface ToxicityEntry {
  risk_score?: number;
  risk_level?: string;
  prediction?: string;
  recommendation?: string;
  logP?: number;
  molecular_weight?: number;
  tpsa?: number;
}

interface ToxicityProfile {
  hERG?: ToxicityEntry;
  hepatotoxicity?: ToxicityEntry;
  mutagenicity?: ToxicityEntry;
  carcinogenicity?: ToxicityEntry;
  [key: string]: ToxicityEntry | undefined;
}

interface SyntheticAccessibility {
  sa_score?: number;
  difficulty?: string;
  estimated_steps?: string;
  recommendation?: string;
}

interface OptimizationSuggestion {
  modification?: string;
  rationale?: string;
  expected_improvement?: string;
}

export interface AdmetResultData {
  smiles?: string;
  compound_name?: string;
  analyzed_at?: string;
  toxicity_profile?: ToxicityProfile;
  synthetic_accessibility?: SyntheticAccessibility;
  optimization_suggestions?: OptimizationSuggestion[];
  // ADMET-AI flat properties
  herg?: number;
  bbb?: number;
  ames?: number;
  dili?: number;
  caco2?: number;
  half_life?: number;
  clearance?: number;
  [key: string]: unknown;
}

interface ExportResultsButtonProps {
  data: AdmetResultData;
  compoundName?: string;
  filename?: string;
  /** "admet" = full ADMET-AI + toxicity panel, "toxicity" = toxicity-only */
  mode?: "admet" | "toxicity";
}

// ── CSV helpers ────────────────────────────────────────────────────────────────

function flattenAdmetToCsv(data: AdmetResultData, mode: "admet" | "toxicity"): string {
  const rows: [string, string | number | undefined][] = [];

  const add = (label: string, value: string | number | undefined) => rows.push([label, value]);

  add("Compound Name", data.compound_name ?? "Unknown");
  add("SMILES", data.smiles ?? "");
  add("Analyzed At", data.analyzed_at ?? new Date().toISOString());

  if (mode === "admet") {
    // ADMET-AI flat properties
    if (data.herg !== undefined) add("hERG (ML)", data.herg);
    if (data.bbb !== undefined) add("BBB Penetration (ML)", data.bbb);
    if (data.ames !== undefined) add("AMES Mutagenicity (ML)", data.ames);
    if (data.dili !== undefined) add("DILI Risk (ML)", data.dili);
    if (data.caco2 !== undefined) add("Caco-2 Permeability (ML)", data.caco2);
    if (data.half_life !== undefined) add("Half-life (ML)", data.half_life);
    if (data.clearance !== undefined) add("Clearance (ML)", data.clearance);
  }

  // Toxicity profile (both modes)
  const tox = data.toxicity_profile ?? {};
  const toxKeys: string[] = ["hERG", "hepatotoxicity", "mutagenicity", "carcinogenicity"];
  for (const key of toxKeys) {
    const entry = tox[key];
    if (!entry) continue;
    const label = key.charAt(0).toUpperCase() + key.slice(1);
    if (entry.risk_score !== undefined) add(`${label} Risk Score`, entry.risk_score);
    if (entry.risk_level) add(`${label} Risk Level`, entry.risk_level);
    if (entry.prediction) add(`${label} Prediction`, entry.prediction);
    if (entry.recommendation) add(`${label} Recommendation`, entry.recommendation);
    if (entry.logP !== undefined) add(`${label} LogP`, entry.logP);
    if (entry.molecular_weight !== undefined) add(`${label} MW`, entry.molecular_weight);
    if (entry.tpsa !== undefined) add(`${label} TPSA`, entry.tpsa);
  }

  // Synthetic accessibility
  const sa = data.synthetic_accessibility;
  if (sa) {
    if (sa.sa_score !== undefined) add("SA Score", sa.sa_score);
    if (sa.difficulty) add("Synthesis Difficulty", sa.difficulty);
    if (sa.estimated_steps) add("Estimated Steps", sa.estimated_steps);
    if (sa.recommendation) add("SA Recommendation", sa.recommendation);
  }

  // Optimization suggestions
  const opts = data.optimization_suggestions ?? [];
  opts.forEach((s, i) => {
    if (s.modification) add(`Optimization ${i + 1} - Modification`, s.modification);
    if (s.rationale) add(`Optimization ${i + 1} - Rationale`, s.rationale);
    if (s.expected_improvement) add(`Optimization ${i + 1} - Expected Improvement`, s.expected_improvement);
  });

  const header = "Property,Value";
  const body = rows.map(([k, v]) => `"${String(k).replace(/"/g, '""')}","${String(v ?? "").replace(/"/g, '""')}"`).join("\n");
  return `${header}\n${body}`;
}

// ── PDF helpers (pure DOM, no server) ─────────────────────────────────────────

function buildPdfHtml(data: AdmetResultData, mode: "admet" | "toxicity", compoundName: string): string {
  const tox = data.toxicity_profile ?? {};
  const sa = data.synthetic_accessibility;
  const opts = data.optimization_suggestions ?? [];
  const now = new Date().toLocaleString();

  const riskBadge = (level?: string) => {
    const color = level === "Low" ? "#16a34a" : level === "Medium" ? "#d97706" : "#dc2626";
    return `<span style="background:${color};color:#fff;padding:2px 8px;border-radius:4px;font-size:12px">${level ?? "N/A"}</span>`;
  };

  const toxRows = (["hERG", "hepatotoxicity", "mutagenicity", "carcinogenicity"] as const)
    .filter(k => tox[k])
    .map(k => {
      const e = tox[k]!;
      return `<tr>
        <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;font-weight:500">${k.charAt(0).toUpperCase() + k.slice(1)}</td>
        <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${e.risk_score?.toFixed(1) ?? e.prediction ?? "—"}</td>
        <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${riskBadge(e.risk_level)}</td>
        <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;font-size:12px;color:#6b7280">${e.recommendation ?? ""}</td>
      </tr>`;
    }).join("");

  const admetRows = mode === "admet" ? [
    ["hERG (ML)", data.herg?.toFixed(4)],
    ["BBB Penetration (ML)", data.bbb?.toFixed(4)],
    ["AMES Mutagenicity (ML)", data.ames?.toFixed(4)],
    ["DILI Risk (ML)", data.dili?.toFixed(4)],
    ["Caco-2 Permeability (ML)", data.caco2?.toFixed(4)],
    ["Half-life (ML)", data.half_life?.toFixed(4)],
    ["Clearance (ML)", data.clearance?.toFixed(4)],
  ].filter(([, v]) => v !== undefined).map(([k, v]) => `<tr>
    <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;font-weight:500">${k}</td>
    <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${v}</td>
  </tr>`).join("") : "";

  const optRows = opts.map((s, i) => `<tr>
    <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${i + 1}</td>
    <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${s.modification ?? ""}</td>
    <td style="padding:6px 12px;border-bottom:1px solid #e5e7eb;font-size:12px;color:#6b7280">${s.rationale ?? ""}</td>
  </tr>`).join("");

  return `<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>${compoundName} — ADMET Report</title>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 40px; color: #111; }
    h1 { font-size: 22px; color: #1e3a5f; margin-bottom: 4px; }
    .meta { font-size: 13px; color: #6b7280; margin-bottom: 24px; }
    h2 { font-size: 15px; color: #1e3a5f; margin: 20px 0 8px; border-bottom: 2px solid #e5e7eb; padding-bottom: 4px; }
    table { width: 100%; border-collapse: collapse; font-size: 13px; }
    th { text-align: left; padding: 8px 12px; background: #f3f4f6; font-weight: 600; }
    .smiles { font-family: monospace; font-size: 12px; background: #f9fafb; padding: 8px 12px; border-radius: 6px; word-break: break-all; }
    .footer { margin-top: 32px; font-size: 11px; color: #9ca3af; border-top: 1px solid #e5e7eb; padding-top: 12px; }
  </style>
</head>
<body>
  <h1>PharmaSight™ — ADMET &amp; Toxicity Report</h1>
  <div class="meta">
    <strong>Compound:</strong> ${compoundName} &nbsp;|&nbsp;
    <strong>Generated:</strong> ${now}
  </div>
  ${data.smiles ? `<div class="smiles"><strong>SMILES:</strong> ${data.smiles}</div>` : ""}

  ${admetRows ? `<h2>ML-Predicted ADMET Properties (ADMET-AI / Chemprop)</h2>
  <table><thead><tr><th>Property</th><th>Value</th></tr></thead><tbody>${admetRows}</tbody></table>` : ""}

  ${toxRows ? `<h2>Toxicity Profile</h2>
  <table><thead><tr><th>Category</th><th>Score / Prediction</th><th>Risk Level</th><th>Recommendation</th></tr></thead><tbody>${toxRows}</tbody></table>` : ""}

  ${sa ? `<h2>Synthetic Accessibility</h2>
  <table><thead><tr><th>Property</th><th>Value</th></tr></thead><tbody>
    ${sa.sa_score !== undefined ? `<tr><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">SA Score</td><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${sa.sa_score.toFixed(2)}</td></tr>` : ""}
    ${sa.difficulty ? `<tr><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">Difficulty</td><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${sa.difficulty}</td></tr>` : ""}
    ${sa.estimated_steps ? `<tr><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">Estimated Steps</td><td style="padding:6px 12px;border-bottom:1px solid #e5e7eb">${sa.estimated_steps}</td></tr>` : ""}
  </tbody></table>` : ""}

  ${optRows ? `<h2>Optimization Suggestions</h2>
  <table><thead><tr><th>#</th><th>Modification</th><th>Rationale</th></tr></thead><tbody>${optRows}</tbody></table>` : ""}

  <div class="footer">Generated by PharmaSight™ Admin Dashboard · ${now}</div>
</body>
</html>`;
}

function downloadCsv(content: string, filename: string) {
  const blob = new Blob([content], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function downloadPdf(html: string, filename: string) {
  const win = window.open("", "_blank");
  if (!win) {
    toast.error("Pop-up blocked. Please allow pop-ups to download PDF.");
    return;
  }
  win.document.write(html);
  win.document.close();
  win.focus();
  // Give browser time to render, then print
  setTimeout(() => {
    win.print();
    win.close();
  }, 600);
}

// ── Component ──────────────────────────────────────────────────────────────────

export function ExportResultsButton({
  data,
  compoundName = "compound",
  filename,
  mode = "admet",
}: ExportResultsButtonProps) {
  const [exporting, setExporting] = useState(false);
  const base = filename ?? `${compoundName.replace(/\s+/g, "_")}_${mode}_${new Date().toISOString().slice(0, 10)}`;

  const handleCsv = () => {
    setExporting(true);
    try {
      const csv = flattenAdmetToCsv(data, mode);
      downloadCsv(csv, `${base}.csv`);
      toast.success("CSV downloaded");
    } catch (e) {
      toast.error("CSV export failed");
    } finally {
      setExporting(false);
    }
  };

  const handlePdf = () => {
    setExporting(true);
    try {
      const html = buildPdfHtml(data, mode, compoundName);
      downloadPdf(html, `${base}.pdf`);
      toast.success("PDF report opened — use your browser's Print → Save as PDF");
    } catch (e) {
      toast.error("PDF export failed");
    } finally {
      setExporting(false);
    }
  };

  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="outline" size="sm" disabled={exporting} className="gap-2">
          {exporting ? (
            <Loader2 className="h-4 w-4 animate-spin" />
          ) : (
            <Download className="h-4 w-4" />
          )}
          Export
        </Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent align="end">
        <DropdownMenuItem onClick={handleCsv} className="gap-2 cursor-pointer">
          <Table2 className="h-4 w-4 text-green-600" />
          Download CSV
        </DropdownMenuItem>
        <DropdownMenuItem onClick={handlePdf} className="gap-2 cursor-pointer">
          <FileText className="h-4 w-4 text-red-500" />
          Export as PDF
        </DropdownMenuItem>
      </DropdownMenuContent>
    </DropdownMenu>
  );
}
