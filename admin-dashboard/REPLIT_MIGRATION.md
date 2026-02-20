# PharmaSight Admin Dashboard - Replit Migration Guide

## Quick Start Commands for Replit

```bash
# 1. Clone the dashboard branch
git clone -b dashboard-manus https://github.com/justincihi/pharmasight-platform.git pharmasight-dashboard
cd pharmasight-dashboard

# 2. Install dependencies
pnpm install

# 3. Set up environment variables in Replit Secrets
# Add these secrets in Replit's Secrets tab:
# - DATABASE_URL (your MySQL/TiDB connection string)
# - JWT_SECRET (generate with: openssl rand -base64 32)
# - GEMINI_API_KEY
# - ANTHROPIC_API_KEY  
# - SONAR_API_KEY (Perplexity)
# - PLATFORM_API_KEY (generate a secure random string)

# 4. Push database schema
pnpm db:push

# 5. Start development server
pnpm dev
```

## Dashboard File Structure

The dashboard is a **complete standalone application** located in:

```
pharmasight-admin-dashboard/
├── client/                    # React frontend
│   ├── src/
│   │   ├── pages/            # Main pages (Home, Dashboard, Analytics, etc.)
│   │   ├── components/       # Reusable UI components
│   │   └── lib/trpc.ts       # API client
├── server/                    # Express + tRPC backend
│   ├── routers.ts            # API endpoints
│   ├── db.ts                 # Database queries
│   ├── multiLLM.ts           # Chat integration
│   ├── retrosynthesis.ts     # Synthesis route generation
│   └── autonomousScheduler.ts # Research engine scheduler
├── drizzle/                   # Database schema
│   └── schema.ts
└── master_analogs.json        # Master analog list (150 compounds)
```

## Key Components

### Main Dashboard Pages
- `client/src/pages/Home.tsx` - Landing page with chatbot
- `client/src/pages/AdminDashboard.tsx` - Main analog discovery dashboard
- `client/src/pages/AnalogDetail.tsx` - Individual analog details with 3D viewer
- `client/src/pages/BatchOperations.tsx` - Bulk analog management
- `client/src/pages/SchedulerDashboard.tsx` - Autonomous engine monitoring

### Backend Services
- `server/routers.ts` - All tRPC API endpoints
- `server/autonomousScheduler.ts` - Scheduled research runs
- `server/platformAPI.ts` - REST API for Python integration

## Replit Configuration

Create `.replit` file:

```toml
run = "pnpm dev"
language = "nodejs"

[nix]
channel = "stable-22_11"

[deployment]
run = ["pnpm", "build", "&&", "pnpm", "start"]
deploymentTarget = "cloudrun"
```

## Environment Variables Required

```env
# Database
DATABASE_URL=mysql://user:pass@host:port/database

# Auth
JWT_SECRET=your-secret-key
OAUTH_SERVER_URL=https://api.manus.im
VITE_OAUTH_PORTAL_URL=https://portal.manus.im

# AI APIs
GEMINI_API_KEY=your-gemini-key
ANTHROPIC_API_KEY=your-anthropic-key
SONAR_API_KEY=your-perplexity-key

# Platform Integration
PLATFORM_API_KEY=your-secure-api-key
BUILT_IN_FORGE_API_KEY=your-manus-api-key
BUILT_IN_FORGE_API_URL=https://forge.manus.im
```

## Connecting to Python Research Engine

The dashboard exposes REST API endpoints at `/api/platform/*`:

```python
import requests

# Import discovered analogs
response = requests.post(
    "https://your-replit-url.repl.co/api/platform/import",
    headers={"Authorization": f"Bearer {PLATFORM_API_KEY}"},
    json={
        "discoveries": [
            {
                "compoundId": "KETAMINE-20251227-A001",
                "compoundName": "N-Ethylketamine",
                "parentCompound": "Ketamine",
                "smiles": "CCC(=O)C1(NC)CCCCC1",
                # ... other fields
            }
        ]
    }
)
```

See `python_integration/pharmasight_dashboard_client.py` for full Python client.

## Troubleshooting

### Port Issues
If port 3000 is in use, change in `server/_core/index.ts`:
```typescript
const PORT = process.env.PORT || 3001;
```

### Database Connection
Ensure DATABASE_URL uses SSL for production:
```
mysql://user:pass@host:port/db?ssl={"rejectUnauthorized":true}
```

### Build Errors
Clear cache and rebuild:
```bash
rm -rf node_modules .next
pnpm install
pnpm build
```

## Working Frontend URL

Current Manus deployment: https://3000-ipncy26as6oda8zafsq4p-991d17ea.us2.manus.computer

After Replit deployment, your URL will be: `https://your-repl-name.your-username.repl.co`
