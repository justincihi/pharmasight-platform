# PharmaSight Dashboard Python Integration

Python client library for integrating the autonomous research engine with the PharmaSight Admin Dashboard.

## Installation

```bash
# Install required dependency
pip install requests

# Set environment variables
export PHARMASIGHT_DASHBOARD_URL="https://your-dashboard.manus.space"
export PLATFORM_API_KEY="your-api-key-here"
```

## Quick Start

```python
from pharmasight_dashboard_client import PharmaSightDashboard

# Initialize client
dashboard = PharmaSightDashboard(
    api_url="https://your-dashboard.manus.space",
    api_key="your-api-key-here"
)

# Check dashboard health
health = dashboard.health_check()
print(f"Status: {health['status']}")

# Import discoveries
discoveries = [
    {
        "compoundId": "KETAMINE-NEW-001",
        "compoundName": "Novel Ketamine Analog",
        "smiles": "CNC1(c2ccccc2F)CCCCC1=O",
        "parentCompound": "Ketamine",
        "mechanismOfAction": "NMDA antagonist",
        "keyDifferences": "Fluorinated derivative",
        "confidence": 90,
        "similarity": 88,
        "safetyScore": 85,
        "efficacyScore": 90,
        "drugLikenessScore": 95,
        "patentStatus": "patent-free",
        "marketValue": "$50M",
        "discoveryMethod": "ai-generation"
    }
]

result = dashboard.import_discoveries(discoveries)
print(f"Imported: {result['imported']}, Skipped: {result['skipped']}")
```

## Integration with Autonomous Research Engine

Add this to your autonomous discovery script:

```python
# At the top of your script
from pharmasight_dashboard_client import PharmaSightDashboard
import os

# Initialize dashboard client
dashboard = PharmaSightDashboard(
    api_url=os.getenv("PHARMASIGHT_DASHBOARD_URL"),
    api_key=os.getenv("PLATFORM_API_KEY")
)

# After discovering new analogs
def process_discoveries(new_analogs):
    """Process and import new analog discoveries"""
    
    # Format discoveries for dashboard
    formatted_discoveries = []
    for analog in new_analogs:
        formatted_discoveries.append({
            "compoundId": analog['id'],
            "compoundName": analog['name'],
            "smiles": analog['smiles'],
            "parentCompound": analog['parent'],
            "mechanismOfAction": analog['mechanism'],
            "keyDifferences": analog['differences'],
            "confidence": analog['confidence_score'],
            "similarity": analog['similarity_score'],
            "safetyScore": analog['safety'],
            "efficacyScore": analog['efficacy'],
            "drugLikenessScore": analog['drug_likeness'],
            "patentStatus": analog['patent_status'],
            "marketValue": analog['market_value'],
            "discoveryMethod": "autonomous-engine"
        })
    
    # Import to dashboard
    result = dashboard.import_discoveries(formatted_discoveries)
    print(f"✅ Imported {result['imported']} new analogs to dashboard")
    
    return result
```

## API Methods

### `health_check()`
Check if the dashboard API is available.

**Returns:** `Dict` with status information

### `import_discoveries(discoveries: List[Dict])`
Import analog discoveries to the dashboard.

**Parameters:**
- `discoveries`: List of discovery dictionaries

**Required fields for each discovery:**
- `compoundId`: Unique identifier
- `compoundName`: Human-readable name
- `smiles`: SMILES notation
- `parentCompound`: Parent compound name
- `mechanismOfAction`: Description of mechanism
- `keyDifferences`: What makes this analog unique
- `confidence`: Confidence score (0-100)
- `similarity`: Similarity to parent (0-100)
- `safetyScore`: Safety score (0-100)
- `efficacyScore`: Efficacy score (0-100)
- `drugLikenessScore`: Drug-likeness score (0-100)
- `patentStatus`: "patent-free", "patent-opportunity", etc.
- `marketValue`: Estimated market value (e.g., "$50M")
- `discoveryMethod`: How it was discovered

**Returns:** `Dict` with import results

### `get_recent_discoveries(limit: int = 10)`
Get recent analog discoveries from the dashboard.

**Parameters:**
- `limit`: Maximum number of discoveries to return

**Returns:** `Dict` with discoveries list

### `get_analog(compound_id: str)`
Get specific analog by compound ID.

**Parameters:**
- `compound_id`: Compound identifier

**Returns:** Analog data or `None` if not found

### `update_analog(compound_id: str, update_data: Dict)`
Update analog data.

**Parameters:**
- `compound_id`: Compound identifier
- `update_data`: Fields to update

**Returns:** `bool` indicating success

## Testing

Run the example script to test the integration:

```bash
python pharmasight_dashboard_client.py
```

## Master Analog File

All discoveries are automatically synced to `master_analogs.json` in the project root. This file serves as a single source of truth for all analog data.

## Troubleshooting

**Connection Error:**
- Verify `PHARMASIGHT_DASHBOARD_URL` is correct
- Ensure dashboard is running and accessible

**Authentication Error:**
- Check `PLATFORM_API_KEY` matches the dashboard configuration
- Verify API key is set in dashboard environment variables

**Import Errors:**
- Check discovery data format matches required fields
- Review dashboard logs for specific error messages
