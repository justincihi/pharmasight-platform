# BioTransformer Integration API

## Overview

BioTransformer 3.0 integration for metabolite prediction in the PharmaSight platform. Predicts Phase I, Phase II, and gut microbiome metabolism pathways.

## Features

- **Phase I Metabolism**: CYP450-mediated oxidation, reduction, hydrolysis
- **Phase II Metabolism**: Glucuronidation, sulfation, acetylation, methylation
- **Gut Microbiome**: Bacterial transformation pathways
- **Batch Processing**: Multiple compounds in a single request
- **Mock Mode**: Fallback with realistic sample data when JAR unavailable

## API Endpoints

### 1. Predict Metabolites

Predict metabolites for a single compound.

**Endpoint**: `biotransformer.predictMetabolites`

**Input**:
```typescript
{
  smiles: string,           // SMILES string of parent compound
  metabolismType: string,   // 'human' | 'cyp450' | 'phase2' | 'gut' | 'ecbased' | 'environmental' | 'superbio'
  steps: number             // Number of transformation steps (1-3)
}
```

**Output**:
```typescript
{
  success: boolean,
  parent_smiles: string,
  metabolism_type: string,
  steps: number,
  num_metabolites: number,
  metabolites: [
    {
      smiles: string,           // SMILES of metabolite
      name: string,             // Descriptive name
      reaction: string,         // Type of transformation
      enzyme: string,           // Enzyme responsible
      molecular_weight: string, // MW in g/mol
      generation: string        // Transformation step (1, 2, 3)
    }
  ],
  mock_mode: boolean,         // true if using mock data
  note?: string               // Optional note about mock mode
}
```

**Example Usage**:
```typescript
const result = await trpc.biotransformer.predictMetabolites.mutate({
  smiles: "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",  // Ketamine
  metabolismType: 'human',
  steps: 2
});

console.log(`Found ${result.num_metabolites} metabolites`);
result.metabolites.forEach(m => {
  console.log(`${m.name}: ${m.reaction} via ${m.enzyme}`);
});
```

### 2. Batch Predict

Process multiple compounds in a single request.

**Endpoint**: `biotransformer.batchPredict`

**Input**:
```typescript
{
  compounds: [
    { smiles: string, id?: string }
  ],
  metabolismType: string,
  steps: number
}
```

**Output**:
```typescript
[
  {
    id: string,
    success: boolean,
    parent_smiles: string,
    metabolites: [...],
    // ... same as single prediction
  }
]
```

**Example Usage**:
```typescript
const compounds = [
  { smiles: "CCN(C1CCCCC1=O)c2cccc(F)c2Cl", id: "ketamine" },
  { smiles: "CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3", id: "atropine" }
];

const results = await trpc.biotransformer.batchPredict.mutate({
  compounds,
  metabolismType: 'human',
  steps: 1
});
```

### 3. Get Metabolism Types

Retrieve available metabolism type options.

**Endpoint**: `biotransformer.getMetabolismTypes`

**Output**:
```typescript
{
  human: "allHuman",
  ecbased: "ecbased",
  cyp450: "cyp450",
  phase2: "phaseII",
  gut: "gut",
  environmental: "env",
  superbio: "superbio"
}
```

### 4. Save Results

Store metabolite prediction results in the database.

**Endpoint**: `biotransformer.saveResults`

**Input**:
```typescript
{
  analogId: number,
  results: object  // Full prediction result object
}
```

**Output**:
```typescript
{
  id: number,
  analogId: number,
  testType: "metabolite_prediction",
  testStatus: "completed",
  results: string,
  runBy: number,
  createdAt: Date
}
```

## Metabolism Types

| Type | Description | Use Case |
|------|-------------|----------|
| `human` | Complete human metabolism (Phase I + II) | General drug metabolism |
| `cyp450` | CYP450 enzyme metabolism only | Phase I oxidation focus |
| `phase2` | Conjugation reactions only | Glucuronidation, sulfation |
| `gut` | Gut microbiome transformations | Oral bioavailability |
| `ecbased` | EC number-based predictions | Enzyme-specific |
| `environmental` | Environmental degradation | Stability assessment |
| `superbio` | Super bio transformation | Comprehensive analysis |

## Common Reactions

### Phase I
- **Hydroxylation**: Addition of -OH group (CYP3A4, CYP2D6)
- **N-oxidation**: Nitrogen oxidation (FMO3)
- **O-dealkylation**: Ether cleavage (CYP3A4)
- **Reduction**: Carbonyl → alcohol (AKR, CBR)

### Phase II
- **Glucuronidation**: Addition of glucuronic acid (UGT1A1, UGT2B7)
- **Sulfation**: Addition of sulfate (SULT1A1)
- **Acetylation**: Addition of acetyl group (NAT2)
- **Methylation**: Addition of methyl group (COMT)

### Gut Microbiome
- **Dehalogenation**: Removal of halogens
- **Reduction**: Nitro → amine
- **Hydrolysis**: Ester/amide cleavage

## Installation

### BioTransformer JAR (Optional)

For production use with actual BioTransformer predictions:

```bash
# Download BioTransformer 3.0 JAR
wget https://bitbucket.org/djoumbou/biotransformer/downloads/BioTransformer3.0.jar \
  -O server/python_modules/BioTransformer3.0.jar

# Verify Java is installed
java -version  # Requires Java 11+
```

### Python Dependencies

```bash
pip install -r server/python_modules/requirements.txt
```

Required packages:
- No additional packages beyond standard library for mock mode
- BioTransformer JAR requires Java 11+ runtime

## Mock Mode

When the BioTransformer JAR is not available, the module runs in **mock mode**:

- Returns realistic sample metabolites based on common transformation patterns
- Includes Phase I (hydroxylation, N-oxidation) and Phase II (glucuronidation, sulfation)
- Useful for development, testing, and demonstration
- Results include `mock_mode: true` flag and explanatory note

**Mock Transformations**:
- Hydroxylation (CYP3A4)
- N-oxidation (FMO3)
- Glucuronidation (UGT1A1) - when steps ≥ 2
- Sulfation (SULT1A1) - when steps ≥ 2
- Dehalogenation (gut microbiome) - when metabolism_type = 'gut'

## Testing

Run BioTransformer integration tests:

```bash
pnpm vitest run server/biotransformer.test.ts
```

Tests cover:
- Metabolism type retrieval
- Single compound prediction (Phase I)
- Multi-step prediction (Phase I + II)
- Batch processing
- Gut microbiome metabolism
- Mock mode verification

## Performance

### Mock Mode
- **Response Time**: <100ms per compound
- **Throughput**: ~500 compounds/minute
- **Memory**: Minimal (<50MB)

### JAR Mode (Estimated)
- **Response Time**: 1-5 seconds per compound
- **Throughput**: ~20-60 compounds/minute
- **Memory**: 2-4GB Java heap recommended
- **Timeout**: 5 minutes per prediction

## Example Workflow

```typescript
// 1. Get available metabolism types
const types = await trpc.biotransformer.getMetabolismTypes.query();

// 2. Predict metabolites for a new analog
const prediction = await trpc.biotransformer.predictMetabolites.mutate({
  smiles: analog.smiles,
  metabolismType: 'human',
  steps: 2
});

// 3. Save results to database
if (prediction.success) {
  await trpc.biotransformer.saveResults.mutate({
    analogId: analog.id,
    results: prediction
  });

  // 4. Display results in UI
  console.log(`Predicted ${prediction.num_metabolites} metabolites:`);
  prediction.metabolites.forEach(m => {
    console.log(`- ${m.name} (${m.enzyme})`);
  });
}
```

## Troubleshooting

### "BioTransformer JAR not found"
- Running in mock mode - download JAR for production use
- Check `server/python_modules/BioTransformer3.0.jar` exists
- Verify file permissions (readable)

### "Python process exited with code 1"
- Check Python dependencies installed: `pip install -r requirements.txt`
- Verify Python 3.8+ available: `python3 --version`
- Check error logs for import failures

### "Prediction timeout (>5 minutes)"
- Reduce complexity or number of steps
- Consider batch processing for multiple compounds
- Increase timeout in `biotransformer_wrapper.py` if needed

### Invalid SMILES
- Validate SMILES format before submission
- Use RDKit for SMILES canonicalization
- Check for special characters or encoding issues

## References

- **BioTransformer**: Djoumbou-Feunang et al. (2019) Journal of Cheminformatics
- **Bitbucket Repository**: https://bitbucket.org/djoumbou/biotransformer
- **Documentation**: https://bitbucket.org/djoumbou/biotransformer/wiki/Home
- **Publication**: https://doi.org/10.1186/s13321-018-0324-5

## License

BioTransformer is developed by the Wishart Lab (University of Alberta) and is available for academic and commercial use. Please cite the original publication when using this tool.
