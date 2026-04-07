# 47 Articles Issue - Root Cause Analysis

## Problem Description

The platform shows "only reading 47 articles each time" with:
- Different results on each run
- Results not saving when navigating to analog creation/research module
- Inconsistent analog counts displayed

## Root Cause Analysis

### 1. API Result Limits

Most research paper APIs have default result limits:
- **PubMed API:** Default 20, max 10,000 per request
- **Semantic Scholar:** Default 100, max 1,000
- **arXiv:** Default 10, max 30,000
- **CrossRef:** Default 20, max 1,000

**Likely Culprit:** The research engine is using an API with a 47-50 result limit.

### 2. Where the Limit Might Be Set

Check these files in the `pharmasight-platform` repo:

```python
# Look for:
- max_results=47
- limit=47  
- top_k=47
- n_results=47
```

Common locations:
- `services/compound-analysis/main.py`
- `src/research_engine.py`
- `src/literature_search.py`
- Any file importing `requests` or API clients

### 3. Why Results Aren't Saving

**Issue:** Results are generated in-memory but not persisted to database.

**Solution:** Ensure the research module calls the dashboard's Platform API:

```python
import requests

# After discovering analogs
response = requests.post(
    "https://your-dashboard-url/api/platform/import",
    headers={"Authorization": f"Bearer {PLATFORM_API_KEY}"},
    json={"discoveries": discovered_analogs}
)
```

### 4. Why Results Differ Each Time

**Cause:** Random sampling or API pagination without cursor tracking.

**Fix:** Implement deterministic search with:
- Fixed seed for random operations
- Cursor-based pagination
- Date-range filtering to avoid duplicates

## Solutions

### Solution 1: Increase API Limit

Find and update the research module:

```python
# BEFORE
results = api.search(query, max_results=47)

# AFTER  
results = api.search(query, max_results=500)
```

### Solution 2: Implement Pagination

```python
all_results = []
offset = 0
batch_size = 100

while offset < 500:  # Get 500 total results
    batch = api.search(query, limit=batch_size, offset=offset)
    all_results.extend(batch)
    offset += batch_size
    if len(batch) < batch_size:
        break  # No more results
```

### Solution 3: Save Results to Database

Add this to the research engine:

```python
def save_discoveries_to_dashboard(discoveries):
    """Save discovered analogs to dashboard database"""
    import requests
    import os
    
    api_key = os.getenv('PLATFORM_API_KEY')
    dashboard_url = os.getenv('DASHBOARD_URL', 'http://localhost:3000')
    
    response = requests.post(
        f"{dashboard_url}/api/platform/import",
        headers={"Authorization": f"Bearer {api_key}"},
        json={"discoveries": discoveries}
    )
    
    return response.json()
```

### Solution 4: Fix Inconsistent Counts

The dashboard shows different counts because:
- **Database count:** Total analogs in `analog_discoveries` table
- **Display count:** Filtered results based on current view/filters
- **Research count:** Number from last research run

**Fix:** Standardize count queries:

```typescript
// In server/routers.ts
const totalCount = await db.select({ count: sql`count(*)` })
  .from(analogDiscoveries);

const displayCount = await db.select({ count: sql`count(*)` })
  .from(analogDiscoveries)
  .where(filters);
```

## Testing Plan

1. **Find the 47 limit:**
   ```bash
   cd /path/to/pharmasight-platform
   grep -r "47\|max_results\|limit" --include="*.py" | grep -v venv
   ```

2. **Increase limit and test:**
   ```python
   # Update the limit to 500
   # Run research module
   # Verify 500 results returned
   ```

3. **Verify persistence:**
   ```bash
   # Check database after run
   mysql> SELECT COUNT(*) FROM analog_discoveries WHERE discovered_at > NOW() - INTERVAL 1 HOUR;
   ```

4. **Test consistency:**
   ```bash
   # Run research 3 times
   # Compare result counts
   # Should be identical if using fixed seed
   ```

## Quick Fixes

### Immediate (5 minutes)

Add to research module:

```python
# At top of file
MAX_RESULTS = int(os.getenv('MAX_RESEARCH_RESULTS', 500))

# In search function
results = api.search(query, max_results=MAX_RESULTS)
```

### Short-term (30 minutes)

1. Find all API calls with result limits
2. Increase limits to 500-1000
3. Add pagination for APIs that require it
4. Implement result caching to avoid re-fetching

### Long-term (2 hours)

1. Build comprehensive research pipeline:
   - Multiple API sources
   - Deduplication logic
   - Incremental updates
   - Result versioning

2. Add research run tracking:
   - Store run metadata (date, query, count)
   - Track which analogs came from which run
   - Enable run comparison

3. Implement smart filtering:
   - Novelty detection (avoid duplicates)
   - Quality scoring (prioritize high-confidence)
   - Relevance ranking (match research goals)

## Files to Check

In `pharmasight-platform` repo:

```
services/
├── analog-generation/main.py
├── compound-analysis/main.py
└── ml-models/main.py

src/
├── research_engine.py (if exists)
├── literature_search.py (if exists)
└── daily_discovery_engine.py (if exists)
```

In `pharmasight-admin-dashboard` repo:

```
server/
├── runAutonomousResearch.ts
├── importDiscoveries.ts
└── autonomousScheduler.ts
```

## Environment Variables to Add

```env
# Research Configuration
MAX_RESEARCH_RESULTS=500
ENABLE_PAGINATION=true
RESEARCH_BATCH_SIZE=100
DEDUPLICATE_RESULTS=true

# API Keys (if using external research APIs)
PUBMED_API_KEY=your-key
SEMANTIC_SCHOLAR_API_KEY=your-key
```

## Monitoring

Add logging to track:
- Number of results fetched per run
- API response times
- Deduplication stats
- Save success/failure rates

```python
logger.info(f"Fetched {len(results)} results from API")
logger.info(f"Deduplicated to {len(unique_results)} unique analogs")
logger.info(f"Saved {saved_count} new analogs to database")
```
