import { describe, it, expect } from 'vitest';

// ── Helpers extracted from the chat.send procedure ──────────────────────────

/** Extract [ID:N] citations from LLM response text */
function extractCitedIds(content: string): number[] {
  const citedIds: number[] = [];
  const idPattern = /\[ID:(\d+)\]/g;
  let match;
  while ((match = idPattern.exec(content)) !== null) {
    const id = parseInt(match[1]);
    if (!citedIds.includes(id)) citedIds.push(id);
  }
  return citedIds;
}

/** Build cited compounds from fetched analog list */
function buildCitedCompounds(
  citedIds: number[],
  analogs: Array<{
    id: number;
    compoundId: string;
    compoundName: string;
    confidenceScore: number;
    patentStatus: string;
    safetyScore: number;
    efficacyScore: number;
  }>
) {
  return citedIds
    .map(id => analogs.find(a => a.id === id))
    .filter(Boolean)
    .map(a => ({
      id: a!.id,
      compoundId: a!.compoundId,
      compoundName: a!.compoundName,
      confidenceScore: a!.confidenceScore,
      patentStatus: a!.patentStatus,
      safetyScore: a!.safetyScore,
      efficacyScore: a!.efficacyScore,
    }));
}

/** Format analog context line for system prompt */
function formatAnalogContextLine(a: {
  id: number;
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  smiles?: string | null;
  safetyScore: number;
  efficacyScore: number;
  confidenceScore: number;
  patentStatus: string;
  therapeuticArea?: string | null;
  createdAt?: Date | null;
}): string {
  return `[ID:${a.id}] ${a.compoundId}: ${a.compoundName} (parent: ${a.parentCompound}) | SMILES: ${a.smiles || 'N/A'} | Safety: ${a.safetyScore}/100 | Efficacy: ${a.efficacyScore}/100 | Confidence: ${a.confidenceScore}% | Patent: ${a.patentStatus} | Therapeutic: ${a.therapeuticArea || 'N/A'} | Discovered: ${a.createdAt ? new Date(a.createdAt).toLocaleDateString() : 'N/A'}`;
}

// ── Test data ────────────────────────────────────────────────────────────────

const MOCK_ANALOGS = [
  {
    id: 1,
    compoundId: 'PST-001',
    compoundName: 'Ketamine Analog A',
    parentCompound: 'Ketamine',
    smiles: 'CC1(NC(=O)c2ccccc2Cl)CCCCC1=O',
    safetyScore: 82,
    efficacyScore: 91,
    confidenceScore: 94,
    patentStatus: 'Patent-free',
    therapeuticArea: 'Anesthesia',
    createdAt: new Date('2025-01-15'),
  },
  {
    id: 2,
    compoundId: 'PST-002',
    compoundName: 'Ketamine Analog B',
    parentCompound: 'Ketamine',
    smiles: 'CC1(NC(=O)c2ccccc2F)CCCCC1=O',
    safetyScore: 75,
    efficacyScore: 88,
    confidenceScore: 87,
    patentStatus: 'Patented',
    therapeuticArea: 'Depression',
    createdAt: new Date('2025-02-20'),
  },
  {
    id: 42,
    compoundId: 'PST-042',
    compoundName: 'Novel NMDA Blocker',
    parentCompound: 'MK-801',
    smiles: 'C1CC2=CC=CC=C2CC1NC(=O)c1ccccc1',
    safetyScore: 90,
    efficacyScore: 95,
    confidenceScore: 97,
    patentStatus: 'Patent-free',
    therapeuticArea: 'Neuroprotection',
    createdAt: new Date('2025-03-10'),
  },
];

// ── Tests ────────────────────────────────────────────────────────────────────

describe('Chat context: extractCitedIds', () => {
  it('extracts a single ID from response', () => {
    const content = 'The best candidate is **Ketamine Analog A** [ID:1] with 94% confidence.';
    expect(extractCitedIds(content)).toEqual([1]);
  });

  it('extracts multiple unique IDs', () => {
    const content = 'Compare [ID:1] and [ID:42] — both are patent-free. Also see [ID:2].';
    expect(extractCitedIds(content)).toEqual([1, 42, 2]);
  });

  it('deduplicates repeated IDs', () => {
    const content = '[ID:1] has high safety. [ID:1] also has high efficacy.';
    expect(extractCitedIds(content)).toEqual([1]);
  });

  it('returns empty array when no IDs present', () => {
    const content = 'There are no compounds matching your query in the database.';
    expect(extractCitedIds(content)).toEqual([]);
  });

  it('handles large IDs correctly', () => {
    const content = 'Compound [ID:12345] was discovered recently.';
    expect(extractCitedIds(content)).toEqual([12345]);
  });
});

describe('Chat context: buildCitedCompounds', () => {
  it('returns matching compound objects for cited IDs', () => {
    const ids = [1, 42];
    const result = buildCitedCompounds(ids, MOCK_ANALOGS);
    expect(result).toHaveLength(2);
    expect(result[0].compoundName).toBe('Ketamine Analog A');
    expect(result[1].compoundName).toBe('Novel NMDA Blocker');
  });

  it('silently skips IDs not found in analog list', () => {
    const ids = [1, 999]; // 999 doesn't exist
    const result = buildCitedCompounds(ids, MOCK_ANALOGS);
    expect(result).toHaveLength(1);
    expect(result[0].id).toBe(1);
  });

  it('returns empty array when no IDs match', () => {
    const result = buildCitedCompounds([100, 200], MOCK_ANALOGS);
    expect(result).toHaveLength(0);
  });

  it('returns all required fields in each cited compound', () => {
    const result = buildCitedCompounds([42], MOCK_ANALOGS);
    const c = result[0];
    expect(c).toHaveProperty('id', 42);
    expect(c).toHaveProperty('compoundId', 'PST-042');
    expect(c).toHaveProperty('compoundName', 'Novel NMDA Blocker');
    expect(c).toHaveProperty('confidenceScore', 97);
    expect(c).toHaveProperty('patentStatus', 'Patent-free');
    expect(c).toHaveProperty('safetyScore', 90);
    expect(c).toHaveProperty('efficacyScore', 95);
  });
});

describe('Chat context: formatAnalogContextLine', () => {
  it('includes [ID:N] prefix for compound ID citation', () => {
    const line = formatAnalogContextLine(MOCK_ANALOGS[0]);
    expect(line).toMatch(/^\[ID:1\]/);
  });

  it('includes compound name and parent', () => {
    const line = formatAnalogContextLine(MOCK_ANALOGS[0]);
    expect(line).toContain('Ketamine Analog A');
    expect(line).toContain('parent: Ketamine');
  });

  it('uses N/A for missing SMILES', () => {
    const analog = { ...MOCK_ANALOGS[0], smiles: null };
    const line = formatAnalogContextLine(analog);
    expect(line).toContain('SMILES: N/A');
  });

  it('uses N/A for missing therapeutic area', () => {
    const analog = { ...MOCK_ANALOGS[0], therapeuticArea: null };
    const line = formatAnalogContextLine(analog);
    expect(line).toContain('Therapeutic: N/A');
  });

  it('includes safety, efficacy, and confidence scores', () => {
    const line = formatAnalogContextLine(MOCK_ANALOGS[2]);
    expect(line).toContain('Safety: 90/100');
    expect(line).toContain('Efficacy: 95/100');
    expect(line).toContain('Confidence: 97%');
  });
});

describe('Chat context: end-to-end citation flow', () => {
  it('extracts IDs from response and maps them to compound objects', () => {
    const llmResponse = `Based on your database, the top candidates are:
1. **Ketamine Analog A** [ID:1] — 94% confidence, patent-free
2. **Novel NMDA Blocker** [ID:42] — 97% confidence, excellent safety profile

I recommend prioritizing [ID:42] for further docking studies.`;

    const ids = extractCitedIds(llmResponse);
    expect(ids).toEqual([1, 42]);

    const cited = buildCitedCompounds(ids, MOCK_ANALOGS);
    expect(cited).toHaveLength(2);
    expect(cited.map(c => c.id)).toEqual([1, 42]);
  });
});
