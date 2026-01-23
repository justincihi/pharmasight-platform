/**
 * Export utilities for analog discoveries
 * Supports SMILES, SDF (3D structure), and PDF reports
 */

interface AnalogData {
  compoundName: string;
  smiles: string;
  parentCompound: string;
  confidenceScore: number;
  similarityScore: number;
  safetyScore: number;
  efficacyScore: number;
  marketValue: number;
  patentStatus: string;
  therapeuticPotential?: string;
  keyDifferences?: string;
}

/**
 * Export SMILES notation as text file
 */
export function exportSMILES(analog: AnalogData): string {
  return `# ${analog.compoundName}
# Parent: ${analog.parentCompound}
# Confidence: ${analog.confidenceScore}%
${analog.smiles}
`;
}

/**
 * Export as SDF (Structure-Data File) format
 * This is a simplified SDF - in production, use RDKit to generate proper 3D coordinates
 */
export function exportSDF(analog: AnalogData): string {
  const sdf = `${analog.compoundName}
  PharmaSight Generated
  
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
> <COMPOUND_NAME>
${analog.compoundName}

> <PARENT_COMPOUND>
${analog.parentCompound}

> <SMILES>
${analog.smiles}

> <CONFIDENCE_SCORE>
${analog.confidenceScore}

> <SIMILARITY_SCORE>
${analog.similarityScore}

> <SAFETY_SCORE>
${analog.safetyScore}

> <EFFICACY_SCORE>
${analog.efficacyScore}

> <MARKET_VALUE>
${analog.marketValue}

> <PATENT_STATUS>
${analog.patentStatus}

${analog.therapeuticPotential ? `> <THERAPEUTIC_POTENTIAL>\n${analog.therapeuticPotential}\n\n` : ''}
${analog.keyDifferences ? `> <KEY_DIFFERENCES>\n${analog.keyDifferences}\n\n` : ''}
$$$$
`;

  return sdf;
}

/**
 * Generate HTML for PDF export
 * This HTML can be converted to PDF using a headless browser or PDF library
 */
export function generatePDFHTML(analog: AnalogData): string {
  return `<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>${analog.compoundName} - PharmaSight Report</title>
  <style>
    body {
      font-family: 'Arial', sans-serif;
      max-width: 800px;
      margin: 40px auto;
      padding: 20px;
      color: #333;
    }
    .header {
      text-align: center;
      border-bottom: 3px solid #2563eb;
      padding-bottom: 20px;
      margin-bottom: 30px;
    }
    .header h1 {
      color: #1e40af;
      margin: 0;
      font-size: 28px;
    }
    .header p {
      color: #64748b;
      margin: 5px 0 0 0;
    }
    .section {
      margin: 25px 0;
    }
    .section-title {
      font-size: 18px;
      font-weight: bold;
      color: #1e40af;
      border-bottom: 2px solid #e2e8f0;
      padding-bottom: 8px;
      margin-bottom: 15px;
    }
    .info-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 15px;
    }
    .info-item {
      padding: 12px;
      background: #f8fafc;
      border-radius: 6px;
    }
    .info-label {
      font-weight: bold;
      color: #475569;
      font-size: 12px;
      text-transform: uppercase;
      margin-bottom: 4px;
    }
    .info-value {
      font-size: 16px;
      color: #0f172a;
    }
    .smiles-box {
      background: #1e293b;
      color: #e2e8f0;
      padding: 15px;
      border-radius: 6px;
      font-family: 'Courier New', monospace;
      word-break: break-all;
      font-size: 14px;
    }
    .badge {
      display: inline-block;
      padding: 4px 12px;
      border-radius: 12px;
      font-size: 14px;
      font-weight: 600;
    }
    .badge-success {
      background: #dcfce7;
      color: #166534;
    }
    .badge-warning {
      background: #fef3c7;
      color: #92400e;
    }
    .badge-info {
      background: #dbeafe;
      color: #1e40af;
    }
    .footer {
      margin-top: 40px;
      padding-top: 20px;
      border-top: 2px solid #e2e8f0;
      text-align: center;
      color: #64748b;
      font-size: 12px;
    }
  </style>
</head>
<body>
  <div class="header">
    <h1>${analog.compoundName}</h1>
    <p>PharmaSight™ Analog Discovery Report</p>
  </div>

  <div class="section">
    <div class="section-title">Compound Information</div>
    <div class="info-grid">
      <div class="info-item">
        <div class="info-label">Parent Compound</div>
        <div class="info-value">${analog.parentCompound}</div>
      </div>
      <div class="info-item">
        <div class="info-label">Confidence Score</div>
        <div class="info-value">
          <span class="badge ${analog.confidenceScore >= 85 ? 'badge-success' : 'badge-warning'}">
            ${analog.confidenceScore}%
          </span>
        </div>
      </div>
      <div class="info-item">
        <div class="info-label">Similarity Score</div>
        <div class="info-value">${analog.similarityScore}%</div>
      </div>
      <div class="info-item">
        <div class="info-label">Patent Status</div>
        <div class="info-value">
          <span class="badge ${analog.patentStatus === 'patent_free' ? 'badge-success' : 'badge-info'}">
            ${analog.patentStatus === 'patent_free' ? 'Patent-Free' : analog.patentStatus === 'patent_opportunity' ? 'Patent Opportunity' : 'Patented'}
          </span>
        </div>
      </div>
    </div>
  </div>

  <div class="section">
    <div class="section-title">Chemical Structure (SMILES)</div>
    <div class="smiles-box">${analog.smiles}</div>
  </div>

  <div class="section">
    <div class="section-title">Predicted Properties</div>
    <div class="info-grid">
      <div class="info-item">
        <div class="info-label">Safety Score</div>
        <div class="info-value">${analog.safetyScore}/100</div>
      </div>
      <div class="info-item">
        <div class="info-label">Efficacy Score</div>
        <div class="info-value">${analog.efficacyScore}/100</div>
      </div>
      <div class="info-item">
        <div class="info-label">Market Value</div>
        <div class="info-value">$${analog.marketValue.toLocaleString()}M</div>
      </div>
    </div>
  </div>

  ${analog.therapeuticPotential ? `
  <div class="section">
    <div class="section-title">Therapeutic Potential</div>
    <p>${analog.therapeuticPotential}</p>
  </div>
  ` : ''}

  ${analog.keyDifferences ? `
  <div class="section">
    <div class="section-title">Key Differences from Parent</div>
    <p>${analog.keyDifferences}</p>
  </div>
  ` : ''}

  <div class="footer">
    <p>Generated by PharmaSight™ Autonomous Discovery Platform</p>
    <p>Report generated on ${new Date().toLocaleDateString()}</p>
  </div>
</body>
</html>`;
}

/**
 * Convert HTML to PDF using Puppeteer (if available)
 * Falls back to returning HTML if Puppeteer is not installed
 */
export async function generatePDF(analog: AnalogData): Promise<Buffer | string> {
  const html = generatePDFHTML(analog);

  // For now, return HTML - PDF generation can be added later with Puppeteer
  // or by using a PDF service
  return html;
}
