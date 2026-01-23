import type { SynthesisRoute } from "./retrosynthesis";

/**
 * Generate CSV content from synthesis routes
 */
export function generateRoutesCSV(routes: SynthesisRoute[]): string {
  const headers = [
    "Route ID",
    "Target Name",
    "Target SMILES",
    "Total Steps",
    "Overall Yield",
    "Total Cost ($)",
    "Feasibility Score",
    "Difficulty",
    "Estimated Time",
    "Summary"
  ];

  const rows = routes.map(route => [
    route.routeId,
    route.targetName,
    route.targetSmiles,
    route.totalSteps.toString(),
    route.overallYield,
    route.totalCost.toFixed(2),
    route.feasibilityScore.toString(),
    route.difficulty,
    route.estimatedTime,
    `"${route.summary.replace(/"/g, '""')}"` // Escape quotes
  ]);

  return [headers.join(","), ...rows.map(row => row.join(","))].join("\n");
}

/**
 * Generate detailed CSV with step-by-step breakdown
 */
export function generateDetailedRoutesCSV(routes: SynthesisRoute[]): string {
  const headers = [
    "Route ID",
    "Target Name",
    "Step Number",
    "Reaction",
    "Reagents",
    "Conditions",
    "Yield",
    "Difficulty",
    "Cost ($)",
    "Notes"
  ];

  const rows: string[] = [];
  
  routes.forEach(route => {
    route.steps.forEach(step => {
      rows.push([
        route.routeId,
        route.targetName,
        step.stepNumber.toString(),
        `"${step.reaction.replace(/"/g, '""')}"`,
        `"${step.reagents.join("; ").replace(/"/g, '""')}"`,
        `"${step.conditions.replace(/"/g, '""')}"`,
        step.yield,
        step.difficulty,
        step.estimatedCost.toFixed(2),
        `"${step.notes.replace(/"/g, '""')}"`
      ].join(","));
    });
  });

  return [headers.join(","), ...rows].join("\n");
}

/**
 * Generate comparison CSV for multiple routes
 */
export function generateComparisonCSV(routes: SynthesisRoute[]): string {
  if (routes.length === 0) return "";

  const headers = ["Metric", ...routes.map((_, i) => `Route ${i + 1}`)];
  
  const metrics = [
    ["Feasibility Score", ...routes.map(r => r.feasibilityScore.toString())],
    ["Total Steps", ...routes.map(r => r.totalSteps.toString())],
    ["Total Cost ($)", ...routes.map(r => r.totalCost.toFixed(2))],
    ["Overall Yield", ...routes.map(r => r.overallYield)],
    ["Difficulty", ...routes.map(r => r.difficulty)],
    ["Estimated Time", ...routes.map(r => r.estimatedTime)],
    ["Starting Materials Count", ...routes.map(r => r.startingMaterials.length.toString())],
  ];

  return [headers.join(","), ...metrics.map(row => row.join(","))].join("\n");
}

/**
 * Generate simple markdown report for synthesis routes
 */
export function generateRoutesMarkdown(routes: SynthesisRoute[]): string {
  let md = `# Synthesis Routes Report\n\n`;
  md += `Generated: ${new Date().toLocaleString()}\n\n`;
  md += `Total Routes: ${routes.length}\n\n`;
  md += `---\n\n`;

  routes.forEach((route, index) => {
    md += `## Route ${index + 1}: ${route.targetName}\n\n`;
    md += `**Route ID:** ${route.routeId}\n\n`;
    md += `**Target SMILES:** \`${route.targetSmiles}\`\n\n`;
    md += `### Overview\n\n`;
    md += `- **Feasibility Score:** ${route.feasibilityScore}/100\n`;
    md += `- **Total Steps:** ${route.totalSteps}\n`;
    md += `- **Total Cost:** $${route.totalCost.toFixed(2)}\n`;
    md += `- **Overall Yield:** ${route.overallYield}\n`;
    md += `- **Difficulty:** ${route.difficulty}\n`;
    md += `- **Estimated Time:** ${route.estimatedTime}\n\n`;
    
    md += `### Summary\n\n${route.summary}\n\n`;
    
    md += `### Starting Materials\n\n`;
    route.startingMaterials.forEach(material => {
      md += `- **${material.name}** (${material.availability}): $${material.cost}\n`;
      md += `  - SMILES: \`${material.smiles}\`\n`;
    });
    md += `\n`;
    
    md += `### Synthesis Steps\n\n`;
    route.steps.forEach(step => {
      md += `#### Step ${step.stepNumber}: ${step.reaction}\n\n`;
      md += `- **Difficulty:** ${step.difficulty}\n`;
      md += `- **Expected Yield:** ${step.yield}\n`;
      md += `- **Estimated Cost:** $${step.estimatedCost}\n`;
      md += `- **Reagents:** ${step.reagents.join(", ")}\n`;
      md += `- **Conditions:** ${step.conditions}\n`;
      md += `- **Notes:** ${step.notes}\n\n`;
    });
    
    md += `---\n\n`;
  });

  return md;
}
