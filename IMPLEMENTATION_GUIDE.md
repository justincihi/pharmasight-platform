# PharmaSight Implementation Guide

## Current Status

### ✅ Completed
1. **Python Service Gateway** - HTTP-based architecture to replace spawn ENOENT errors
2. **Metabolite Router** - Biotransformer integration for metabolite prediction
3. **Service Router Updates** - Docking and toxicity now use gateway
4. **Result Caching** - SHA256-based SMILES hashing with TTL
5. **Analysis Test Suite** - 21 comprehensive tests
6. **Lead Optimization Router** - tRPC procedures for dragonfly_gen integration

### ⚠️ Critical - Requires Setup

#### 1. Python Microservices Deployment
The application now expects Python services running on these URLs (configurable via env vars):

```
PYTHON_DOCKING_SERVICE_URL=http://localhost:5001
PYTHON_ADMET_SERVICE_URL=http://localhost:5002
PYTHON_TOXICITY_SERVICE_URL=http://localhost:5003
PYTHON_METABOLITES_SERVICE_URL=http://localhost:5004
PYTHON_LEAD_OPT_SERVICE_URL=http://localhost:5005
PYTHON_BIONEMO_SERVICE_URL=http://localhost:5006
```

**Setup Options:**
- **Option A (Recommended)**: Deploy to Manus Cloud Computer ($10/month) with persistent Python environment
- **Option B**: Use local machine with Manus Desktop client
- **Option C**: Docker containers on any cloud provider

**Required Python Services:**
Each service needs a Flask/FastAPI endpoint with `/health` and specific POST endpoints:
- `/dock` - Docking simulation
- `/predict` - ADMET/Toxicity prediction
- `/optimize` - Lead optimization
- `/process` - BioNemo protein processing

#### 2. Autonomous Research Engine Activation
**Current Issues:**
- `research_goals.json` file path error
- LLM API not configured for trend analysis
- Manual trigger button not implemented

**To Fix:**
1. Create `/home/ubuntu/pharmasight-admin-dashboard/data/research_goals.json`:
```json
{
  "goals": [
    {
      "id": "psychedelics",
      "name": "Psychedelics Research",
      "description": "Discover novel psychedelic compounds",
      "therapeutic_areas": ["psychiatry", "neurology"],
      "enabled": true
    }
  ],
  "lastUpdated": "2026-05-29T00:00:00Z"
}
```

2. Set LLM API key for trend analysis:
```bash
PERPLEXITY_API_KEY=your_key_here  # For medical trend analysis
# OR
GEMINI_API_KEY=your_key_here      # Alternative
```

3. Create ResearchEnginePanel component with:
   - Manual "Run Research" button
   - Progress indicator
   - Results viewer
   - History log

#### 3. Biotransformer UI Integration
**What's Done:**
- `metaboliteRouter.ts` with 4 procedures (predict, getHistory, analyzeStability, compareProfiles)
- Service gateway integration

**What's Needed:**
- Add Metabolite tab to CompoundTesting page
- Create MetaboliteViewer component showing:
  - Phase I/II/III metabolites
  - Metabolite structures
  - Stability scores
  - Toxicity predictions
- Add metabolite pathway visualization

#### 4. Dragonfly Lead Optimization UI
**What's Done:**
- `leadOptimizationRouter.ts` with 6 procedures
- `LeadOptimizationPanel.tsx` component

**What's Needed:**
- Wire LeadOptimizationPanel into CompoundTesting page
- Add parameter controls for:
  - Number of analogs (5-100)
  - Optimization objectives (potency, selectivity, ADMET)
  - Constraint values (MW, logP, HBD, HBA)
- Create lead ranking visualization
- Add structure comparison viewer

#### 5. BioNemo Protein Analysis
**What's Done:**
- Service gateway support
- Mock responses

**What's Needed:**
1. Install BioNemo:
```bash
pip install bionemo
```

2. Create BioNemoRouter with procedures for:
   - Protein embedding generation
   - Target prediction
   - Binding site identification
   - Protein-ligand interaction prediction

3. Create ProteinAnalysisPanel component showing:
   - Protein sequence input
   - Embedding visualization
   - Predicted binding sites
   - Interaction heatmaps

## Environment Variables Required

```bash
# Python Services (set these to point to your persistent compute)
PYTHON_DOCKING_SERVICE_URL=http://your-service:5001
PYTHON_ADMET_SERVICE_URL=http://your-service:5002
PYTHON_TOXICITY_SERVICE_URL=http://your-service:5003
PYTHON_METABOLITES_SERVICE_URL=http://your-service:5004
PYTHON_LEAD_OPT_SERVICE_URL=http://your-service:5005
PYTHON_BIONEMO_SERVICE_URL=http://your-service:5006

# LLM APIs for Autonomous Research
PERPLEXITY_API_KEY=your_key_here
GEMINI_API_KEY=your_key_here

# Database
DATABASE_URL=mysql://user:pass@host/db

# OAuth
VITE_APP_ID=your_app_id
OAUTH_SERVER_URL=https://api.manus.im
JWT_SECRET=your_secret
```

## Database Migrations Required

```bash
# Add metabolite_prediction analysis type
ALTER TABLE analysis_results MODIFY analysis_type ENUM('docking', 'toxicity', 'admet', 'pkpd', 'metabolite_prediction', 'lead_optimization', 'bionemo');

# Create research results table
CREATE TABLE research_results (
  id INT AUTO_INCREMENT PRIMARY KEY,
  goal_id VARCHAR(128) NOT NULL,
  research_type VARCHAR(64) NOT NULL,
  findings JSON NOT NULL,
  confidence_score INT,
  source VARCHAR(128),
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
);

# Create research execution log
CREATE TABLE research_execution_log (
  id INT AUTO_INCREMENT PRIMARY KEY,
  goal_id VARCHAR(128),
  status ENUM('pending', 'running', 'completed', 'failed'),
  started_at TIMESTAMP,
  completed_at TIMESTAMP,
  error_message TEXT,
  results_count INT
);
```

## Testing Checklist

- [ ] Python services health check endpoints responding
- [ ] Docking service returning valid binding affinity scores
- [ ] ADMET service returning property predictions
- [ ] Toxicity service returning toxicity profiles
- [ ] Metabolite service returning metabolite structures
- [ ] Lead optimization service returning optimized analogs
- [ ] BioNemo service returning protein embeddings
- [ ] Research goals loading from JSON file
- [ ] Autonomous research engine triggering manually
- [ ] Results persisting to database
- [ ] Cache hit rates tracking correctly

## Next Steps (Priority Order)

1. **Deploy Python Services** - Use persistent compute to run microservices
2. **Configure Environment Variables** - Point to Python service URLs
3. **Activate Research Engine** - Fix file paths, set LLM API keys
4. **Build UI Components** - Add metabolite, lead optimization, protein analysis panels
5. **Run Integration Tests** - Verify all services working end-to-end
6. **Deploy to Production** - Publish with all features enabled

## Support

For issues with:
- **Python service spawning**: Check service gateway logs, verify HTTP endpoints
- **Research engine**: Verify research_goals.json exists, LLM API keys set
- **UI components**: Check browser console for tRPC errors
- **Database**: Run migrations, verify schema matches types

## Architecture Diagram

```
┌─────────────────────────────────────────┐
│     PharmaSight Admin Dashboard         │
│     (Node.js + React + tRPC)            │
└──────────────┬──────────────────────────┘
               │
               ├─→ Python Service Gateway
               │   (HTTP client with fallback)
               │
       ┌───────┴────────┬────────┬────────┬────────┬────────┐
       │                │        │        │        │        │
   Docking         ADMET    Toxicity  Metabolite  Lead    BioNemo
   Service        Service   Service   Service    Optim    Service
   (5001)         (5002)    (5003)    (5004)     (5005)   (5006)
   
   All running on persistent compute (Cloud Computer or local machine)
```
