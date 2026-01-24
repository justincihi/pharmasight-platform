# PharmaSight Platform API Documentation

REST API endpoints for integrating the Python research engine with the Admin Dashboard.

## Authentication

All endpoints require an API key passed in the request:
- **Query parameter**: `?apiKey=YOUR_API_KEY`
- **Request body**: `{ "apiKey": "YOUR_API_KEY", ... }`

Set the API key in environment variable: `PLATFORM_API_KEY`

## Endpoints

### Health Check

**GET** `/api/platform/health`

Check if the API is running.

**Response:**
```json
{
  "status": "ok",
  "service": "pharmasight-admin-dashboard",
  "timestamp": "2025-12-24T19:00:00.000Z"
}
```

---

### Import Analog Discoveries

**POST** `/api/platform/discoveries/import`

Import analog discoveries from the research engine into the dashboard database.

**Request Body:**
```json
{
  "apiKey": "YOUR_API_KEY",
  "discoveries": [
    {
      "compoundId": "KETAMINE-20251224-A001",
      "compoundName": "Fluoroketamine Analog 1",
      "smiles": "CNC1(c2cccc(F)c2Cl)CCCCC1=O",
      "parentCompound": "Ketamine",
      "mechanismOfAction": "Enhanced NMDA receptor antagonism",
      "keyDifferences": "Fluorinated derivative with improved binding",
      "confidence": 92,
      "similarity": 88,
      "safetyScore": 85,
      "efficacyScore": 90,
      "drugLikenessScore": 95,
      "patentStatus": "patent-free",
      "marketValue": "$65M",
      "discoveryMethod": "autonomous-engine"
    }
  ]
}
```

**Response:**
```json
{
  "success": true,
  "imported": 1,
  "skipped": 0,
  "errors": 0,
  "details": {
    "imported": ["KETAMINE-20251224-A001"],
    "errors": []
  }
}
```

---

### Get Recent Discoveries

**GET** `/api/platform/discoveries/recent?apiKey=YOUR_API_KEY&limit=10`

Retrieve recent analog discoveries from the dashboard.

**Query Parameters:**
- `apiKey` (required): API key for authentication
- `limit` (optional): Number of discoveries to return (default: 10)

**Response:**
```json
{
  "success": true,
  "count": 10,
  "discoveries": [
    {
      "id": 1,
      "compoundId": "KETAMINE-20251224-A001",
      "compoundName": "Fluoroketamine Analog 1",
      "smiles": "CNC1(c2cccc(F)c2Cl)CCCCC1=O",
      "parentCompound": "Ketamine",
      "confidenceScore": 92,
      "similarityScore": 88,
      "discoveredAt": "2025-12-24T19:00:00.000Z",
      ...
    }
  ]
}
```

---

### Get Analog by Compound ID

**GET** `/api/platform/analogs/:compoundId?apiKey=YOUR_API_KEY`

Retrieve a specific analog by its compound ID.

**Example:**
```
GET /api/platform/analogs/KETAMINE-20251224-A001?apiKey=YOUR_API_KEY
```

**Response:**
```json
{
  "success": true,
  "analog": {
    "id": 1,
    "compoundId": "KETAMINE-20251224-A001",
    "compoundName": "Fluoroketamine Analog 1",
    "smiles": "CNC1(c2cccc(F)c2Cl)CCCCC1=O",
    ...
  }
}
```

---

### Update Analog Data

**PUT** `/api/platform/analogs/:compoundId`

Update an existing analog's data.

**Request Body:**
```json
{
  "apiKey": "YOUR_API_KEY",
  "safetyScore": 90,
  "efficacyScore": 92,
  "patentStatus": "patent-opportunity"
}
```

**Response:**
```json
{
  "success": true,
  "message": "Analog KETAMINE-20251224-A001 updated successfully"
}
```

---

## Python Integration Example

```python
import requests
import json

# Configuration
API_BASE_URL = "https://your-dashboard-url.com"
API_KEY = "your-api-key-here"

# Import discoveries
def import_discoveries(discoveries):
    url = f"{API_BASE_URL}/api/platform/discoveries/import"
    payload = {
        "apiKey": API_KEY,
        "discoveries": discoveries
    }
    response = requests.post(url, json=payload)
    return response.json()

# Get recent discoveries
def get_recent_discoveries(limit=10):
    url = f"{API_BASE_URL}/api/platform/discoveries/recent"
    params = {"apiKey": API_KEY, "limit": limit}
    response = requests.get(url, params=params)
    return response.json()

# Example usage
if __name__ == "__main__":
    # Import a discovery
    new_discoveries = [{
        "compoundId": "KETAMINE-20251224-A001",
        "compoundName": "Fluoroketamine Analog 1",
        "smiles": "CNC1(c2cccc(F)c2Cl)CCCCC1=O",
        "parentCompound": "Ketamine",
        "confidence": 92,
        "similarity": 88,
        "safetyScore": 85,
        "efficacyScore": 90,
        "drugLikenessScore": 95,
        "patentStatus": "patent-free",
        "marketValue": "$65M"
    }]
    
    result = import_discoveries(new_discoveries)
    print(f"Imported: {result['imported']} discoveries")
    
    # Get recent discoveries
    recent = get_recent_discoveries(limit=5)
    print(f"Found {recent['count']} recent discoveries")
```

---

## Error Responses

All endpoints return standard error responses:

```json
{
  "error": "Error message",
  "message": "Detailed error information"
}
```

**Common HTTP Status Codes:**
- `200`: Success
- `400`: Bad Request (invalid input)
- `401`: Unauthorized (invalid API key)
- `404`: Not Found
- `500`: Internal Server Error

---

## Security Notes

1. **API Key**: Store the API key securely in environment variables
2. **HTTPS**: Always use HTTPS in production
3. **Rate Limiting**: Consider implementing rate limiting for production use
4. **Validation**: All input is validated before processing
5. **Logging**: All API calls are logged for audit purposes
