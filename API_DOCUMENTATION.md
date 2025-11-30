# PharmaSight™ API Documentation

## Overview

The PharmaSight platform provides a comprehensive RESTful API for pharmaceutical research and drug discovery. All requests are routed through the API Gateway at `http://localhost:8080`.

## Authentication

Most endpoints require authentication using JWT tokens.

### Obtain Access Token

**Endpoint:** `POST /auth-service/token`

**Request:**
```bash
curl -X POST http://localhost:8080/auth-service/token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=researcher&password=password"
```

**Response:**
```json
{
  "access_token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "token_type": "bearer"
}
```

**Default Users:**
- **Username:** `researcher`, **Password:** `password`, **Role:** researcher
- **Username:** `admin`, **Password:** `password`, **Role:** admin

### Use Access Token

Include the token in the Authorization header:
```bash
curl -H "Authorization: Bearer YOUR_TOKEN_HERE" http://localhost:8080/auth-service/users/me
```

---

## Health Checks

### Gateway Health Check

**Endpoint:** `GET /health`

**Description:** Returns the health status of the API Gateway and all downstream services.

**Request:**
```bash
curl http://localhost:8080/health
```

**Response:**
```json
{
  "gateway": "healthy",
  "services": {
    "compound-analysis": {"status": "healthy", "status_code": 200},
    "analog-generation": {"status": "healthy", "status_code": 200},
    "ml-models": {"status": "healthy", "status_code": 200},
    "quantum-calculator": {"status": "healthy", "status_code": 200},
    "auth-service": {"status": "healthy", "status_code": 200}
  },
  "overall_status": "healthy"
}
```

---

## Compound Analysis Service

### Analyze Compound

**Endpoint:** `POST /compound-analysis/analyze`

**Description:** Analyzes a chemical compound and returns its molecular properties, including RDKit calculations and 2D structure visualization.

**Request:**
```bash
curl -X POST http://localhost:8080/compound-analysis/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "psilocybin",
    "analysis_type": "full"
  }'
```

**Request Body:**
- `compound` (string, required): Compound name or SMILES string
- `analysis_type` (string, optional): Type of analysis (default: "full")

**Response:**
```json
{
  "base_properties": {
    "name": "Psilocybin",
    "smiles": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
    "molecular_weight": 284.25,
    "therapeutic_area": "Psychedelic Therapy",
    "status": "Phase II Clinical Trials"
  },
  "rdkit_analysis": {
    "molecular_weight": 284.25,
    "logp": 1.74,
    "tpsa": 86.15,
    "num_h_donors": 2,
    "num_h_acceptors": 6,
    "num_rotatable_bonds": 5,
    "formal_charge": 0
  },
  "svg_image": "<svg>...</svg>"
}
```

---

## Analog Generation Service

### Generate Analogs

**Endpoint:** `POST /analog-generation/generate`

**Description:** Generates structural analogs for a parent compound with similarity scoring and patent information.

**Request:**
```bash
curl -X POST http://localhost:8080/analog-generation/generate \
  -H "Content-Type: application/json" \
  -d '{
    "parent_compound": "psilocybin",
    "target_properties": "all"
  }'
```

**Request Body:**
- `parent_compound` (string, required): Parent compound name or SMILES
- `target_properties` (string, optional): Filter criteria
  - `"all"` - Return all analogs (default)
  - `"patent-free"` - Return only patent-free analogs
  - `"high-similarity"` - Return analogs with similarity > 0.9

**Response:**
```json
{
  "parent_compound": "psilocybin",
  "analogs": [
    {
      "name": "4-HO-MET",
      "smiles": "CCN(C)CCc1c[nH]c2ccc(O)cc12",
      "similarity": 0.92,
      "patent_status": "Expired"
    },
    {
      "name": "4-AcO-DMT",
      "smiles": "CC(=O)Oc1ccc2c(c1)c(CCN(C)C)c[nH]2",
      "similarity": 0.88,
      "patent_status": "Expired"
    }
  ],
  "count": 2
}
```

---

## ML Models Service

### Predict ADMET Properties

**Endpoint:** `POST /ml-models/predict/admet`

**Description:** Predicts Absorption, Distribution, Metabolism, Excretion, and Toxicity properties for multiple molecules.

**Request:**
```bash
curl -X POST http://localhost:8080/ml-models/predict/admet \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": [
      "CC(=O)Oc1ccccc1C(=O)O",
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ]
  }'
```

**Request Body:**
- `smiles_list` (array of strings, required): List of SMILES strings

**Response:**
```json
{
  "warning": "ADMET model not loaded. Using heuristic predictions.",
  "predictions": {
    "CC(=O)Oc1ccccc1C(=O)O": {
      "absorption": 0.85,
      "toxicity": 0.15,
      "lipinski_violations": 0,
      "bioavailability_score": 0.72
    },
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": {
      "absorption": 0.92,
      "toxicity": 0.08,
      "lipinski_violations": 0,
      "bioavailability_score": 0.85
    }
  }
}
```

### Predict Molecular Properties

**Endpoint:** `POST /ml-models/predict/properties`

**Description:** Predicts various molecular properties for a single molecule.

**Request:**
```bash
curl -X POST http://localhost:8080/ml-models/predict/properties \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O"
  }'
```

**Request Body:**
- `smiles` (string, required): SMILES string

**Response:**
```json
{
  "warning": "Property model not loaded. Using RDKit descriptors.",
  "molecular_weight": 180.16,
  "logp": 1.31,
  "tpsa": 63.6,
  "hbd": 1,
  "hba": 4
}
```

---

## Quantum Calculator Service

### Calculate Quantum Properties

**Endpoint:** `POST /quantum-calculator/calculate`

**Description:** Performs DFT calculations to determine electronic properties of molecules.

**Request:**
```bash
curl -X POST http://localhost:8080/quantum-calculator/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "mol_geometry": "O 0 0 0; H 0 1 0; H 0 0 1",
    "basis": "sto-3g",
    "xc_functional": "b3lyp"
  }'
```

**Request Body:**
- `mol_geometry` (string, required): Molecular geometry in XYZ format
- `basis` (string, optional): Basis set (default: "sto-3g")
- `xc_functional` (string, optional): Exchange-correlation functional (default: "b3lyp")

**Response:**
```json
{
  "total_energy_hartree": -75.98,
  "homo_lumo_gap_ev": 12.5,
  "dipole_moment_debye": 1.85,
  "calculation_details": {
    "basis_set": "sto-3g",
    "xc_functional": "b3lyp",
    "num_basis_functions": 7
  }
}
```

---

## Error Responses

All services return standardized error responses:

### 400 Bad Request
```json
{
  "detail": "Invalid SMILES string"
}
```

### 401 Unauthorized
```json
{
  "detail": "Could not validate credentials"
}
```

### 404 Not Found
```json
{
  "detail": "Service 'unknown-service' not found"
}
```

### 503 Service Unavailable
```json
{
  "detail": "Service unavailable: compound-analysis"
}
```

### 504 Gateway Timeout
```json
{
  "detail": "Service timeout: quantum-calculator"
}
```

---

## Rate Limiting

Currently, there are no rate limits enforced. In production, consider implementing rate limiting using:
- API Gateway middleware
- Redis-based rate limiting
- Token bucket algorithm

---

## Best Practices

1. **Always check service health** before making requests
2. **Cache JWT tokens** for the duration of their validity (30 minutes)
3. **Use batch endpoints** when processing multiple molecules
4. **Handle timeouts gracefully**, especially for quantum calculations
5. **Validate SMILES strings** on the client side when possible
6. **Use appropriate timeout values** for long-running calculations

---

## SDK Examples

### Python Example

```python
import requests

# Authenticate
token_response = requests.post(
    "http://localhost:8080/auth-service/token",
    data={"username": "researcher", "password": "password"}
)
token = token_response.json()["access_token"]
headers = {"Authorization": f"Bearer {token}"}

# Analyze compound
response = requests.post(
    "http://localhost:8080/compound-analysis/analyze",
    headers=headers,
    json={"compound": "psilocybin"}
)
print(response.json())
```

### JavaScript Example

```javascript
// Authenticate
const tokenResponse = await fetch('http://localhost:8080/auth-service/token', {
  method: 'POST',
  headers: {'Content-Type': 'application/x-www-form-urlencoded'},
  body: 'username=researcher&password=password'
});
const {access_token} = await tokenResponse.json();

// Analyze compound
const response = await fetch('http://localhost:8080/compound-analysis/analyze', {
  method: 'POST',
  headers: {
    'Authorization': `Bearer ${access_token}`,
    'Content-Type': 'application/json'
  },
  body: JSON.stringify({compound: 'psilocybin'})
});
const data = await response.json();
console.log(data);
```

---

## OpenAPI/Swagger Documentation

Each service exposes its own OpenAPI documentation:

- **API Gateway:** http://localhost:8080/docs
- **Compound Analysis:** http://localhost:8001/docs
- **Analog Generation:** http://localhost:8002/docs
- **ML Models:** http://localhost:8003/docs
- **Quantum Calculator:** http://localhost:8004/docs
- **Auth Service:** http://localhost:8005/docs

---

## Support

For issues or questions:
1. Check service health endpoints
2. Review logs: `docker-compose logs -f [service-name]`
3. Consult the repository documentation
4. Open an issue on GitHub

---

**Last Updated:** November 2024  
**Version:** 1.0.0
