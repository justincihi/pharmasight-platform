# BioTransformer Service

REST API wrapper for BioTransformer 3.0 - Metabolite prediction tool for pharmaceutical compounds.

## Overview

BioTransformer predicts small molecule metabolism in:
- **Human** (Phase I & II)
- **CYP450** enzymes
- **Phase II** conjugation
- **Gut microbiome**
- **Environmental** degradation

## API Endpoints

### Health Check
```bash
GET /health

Response:
{
  "status": "healthy",
  "service": "biotransformer",
  "mock_mode": false,
  "jar_exists": true
}
```

### Predict Metabolites
```bash
POST /predict

Request:
{
  "smiles": "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",
  "metabolism_type": "human",
  "steps": 2
}

Response:
{
  "success": true,
  "parent_smiles": "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",
  "metabolism_type": "human",
  "steps": 2,
  "num_metabolites": 15,
  "metabolites": [
    {
      "smiles": "...",
      "name": "Hydroxylated metabolite",
      "reaction": "Hydroxylation",
      "enzyme": "CYP3A4",
      "molecular_weight": "285.73",
      "generation": "1"
    }
  ]
}
```

### Metabolism Types
```bash
GET /metabolism-types

Response:
{
  "metabolism_types": [
    {
      "id": "human",
      "name": "allHuman",
      "description": "Comprehensive human metabolism (Phase I & II)"
    },
    ...
  ]
}
```

### Batch Prediction
```bash
POST /batch

Request:
{
  "compounds": [
    {"id": "COMP-001", "smiles": "CCO"},
    {"id": "COMP-002", "smiles": "CC(C)O"}
  ],
  "metabolism_type": "cyp450",
  "steps": 1
}

Response:
{
  "success": true,
  "total_compounds": 2,
  "results": [...]
}
```

## Docker Usage

### Build
```bash
docker build -t pharmasight/biotransformer:latest .
```

### Run
```bash
docker run -p 8007:8000 \
  -e JAVA_OPTS="-Xmx4g" \
  -e DEBUG=False \
  pharmasight/biotransformer:latest
```

### With docker-compose
```yaml
biotransformer:
  build: ./services/biotransformer
  ports:
    - "8007:8000"
  environment:
    - JAVA_OPTS=-Xmx4g
    - PREDICTION_TIMEOUT=300
  healthcheck:
    test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
    interval: 30s
    timeout: 10s
    retries: 3
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PORT` | 8000 | Service port |
| `BIOTRANSFORMER_JAR` | /opt/BioTransformer3.0.jar | Path to JAR file |
| `JAVA_OPTS` | -Xmx4g | Java heap size |
| `PREDICTION_TIMEOUT` | 300 | Timeout in seconds |
| `MAX_BATCH_SIZE` | 100 | Max compounds in batch |
| `DEBUG` | False | Enable debug mode |

## Mock Mode

If BioTransformer JAR is not available, the service runs in **mock mode**:
- Returns realistic sample metabolites
- Useful for testing integration
- Clearly labeled in responses: `"mock_mode": true`

## Integration with PharmaSight

### Python Integration
```python
import requests

def predict_metabolites(smiles: str) -> dict:
    response = requests.post(
        'http://biotransformer:8000/predict',
        json={
            'smiles': smiles,
            'metabolism_type': 'human',
            'steps': 2
        }
    )
    return response.json()

# Example
result = predict_metabolites("CCN(C1CCCCC1=O)c2cccc(F)c2Cl")
print(f"Found {result['num_metabolites']} metabolites")
```

### TypeScript Integration (tRPC)
```typescript
// admin-dashboard/server/metaboliteRouter.ts
export const metaboliteRouter = router({
  predict: publicProcedure
    .input(z.object({
      smiles: z.string(),
      metabolismType: z.enum(['human', 'cyp450', 'phase2', 'gut']),
      steps: z.number().min(1).max(3)
    }))
    .mutation(async ({ input }) => {
      const response = await fetch('http://biotransformer:8000/predict', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          smiles: input.smiles,
          metabolism_type: input.metabolismType,
          steps: input.steps
        })
      });

      return response.json();
    })
});
```

## Performance

- **Single compound**: 10-60 seconds depending on complexity
- **Batch mode**: ~30-90 seconds per compound
- **Memory**: 2-4GB heap recommended
- **Steps**: Each step multiplies prediction time

## Metabolism Types Details

| Type | BioTransformer Code | Description |
|------|---------------------|-------------|
| `human` | `allHuman` | Phase I + Phase II (comprehensive) |
| `cyp450` | `cyp450` | Cytochrome P450 only |
| `phase2` | `phaseII` | Conjugation reactions only |
| `gut` | `gut` | Gut microbiome metabolism |
| `environmental` | `env` | Environmental degradation |
| `superbio` | `superbio` | All pathways combined |
| `ecbased` | `ecbased` | EC number-based prediction |

## Troubleshooting

### JAR Not Found
- Download manually from: https://bitbucket.org/djoumbou/biotransformer/downloads/
- Place in `/opt/BioTransformer3.0.jar`
- Or set `BIOTRANSFORMER_JAR` environment variable

### Out of Memory
- Increase Java heap: `JAVA_OPTS=-Xmx8g`
- Reduce batch size
- Use fewer prediction steps

### Slow Predictions
- Normal for complex molecules with 2-3 steps
- Consider using `cyp450` or `phase2` instead of `human`
- Use batch mode for multiple compounds

## References

- BioTransformer: https://bitbucket.org/djoumbou/biotransformer/
- Publication: https://doi.org/10.1186/s13321-018-0324-5
- Documentation: https://bitbucket.org/djoumbou/biotransformer/wiki/Home

## License

This wrapper is part of PharmaSight Platform.
BioTransformer is licensed under GPL-3.0 (see original repository).
