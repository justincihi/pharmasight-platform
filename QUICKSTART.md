# PharmaSight Combined - Quick Start Guide

Get PharmaSight Platform running in 5 minutes!

## Prerequisites

- Docker and Docker Compose installed
- Git installed
- 4GB RAM available

## Local Development (Fastest)

```bash
# 1. Clone and checkout Combined branch
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
git checkout Combined

# 2. Start all services
docker-compose up --build

# 3. Access applications
# Main app: http://localhost:5000
# Admin dashboard: http://localhost:3000
```

That's it! The platform is now running locally.

## Deploy to Render (Production)

### Step 1: Create PostgreSQL Database

1. Go to Render Dashboard → New → PostgreSQL
2. Name: `pharmasight-db`
3. Copy the **Internal Database URL**

### Step 2: Create Redis Instance

1. Go to Render Dashboard → New → Redis
2. Name: `pharmasight-redis`
3. Copy the **Internal Redis URL**

### Step 3: Deploy Flask Application

1. Go to Render Dashboard → New → Web Service
2. Connect your GitHub repository
3. Configure:
   - **Name:** `pharmasight-app`
   - **Branch:** `Combined`
   - **Root Directory:** `/`
   - **Environment:** `Docker`
   - **Docker Command:** (leave default)

4. Add Environment Variables:
   ```
   DATABASE_URL=<paste-internal-database-url>
   REDIS_URL=<paste-internal-redis-url>
   FLASK_ENV=production
   PORT=5000
   ```

5. Click **Create Web Service**

### Step 4: Deploy Admin Dashboard

1. Go to Render Dashboard → New → Web Service
2. Connect your GitHub repository
3. Configure:
   - **Name:** `pharmasight-admin`
   - **Branch:** `Combined`
   - **Root Directory:** `/admin-dashboard`
   - **Environment:** `Docker`
   - **Dockerfile Path:** `admin-dashboard/Dockerfile`

4. Add Environment Variables:
   ```
   DATABASE_URL=<same-database-url-as-flask>
   FLASK_API_URL=https://pharmasight-app.onrender.com
   JWT_SECRET=<generate-random-secret>
   PORT=3000
   NODE_ENV=production
   ```

5. Click **Create Web Service**

### Step 5: Configure Custom Domains (Optional)

1. In Flask app settings → Custom Domains → Add `app.pharmasight.com`
2. In Admin app settings → Custom Domains → Add `admin.pharmasight.com`
3. Update your DNS records as instructed by Render

## First Time Setup

### Create Admin User

Access the admin dashboard and create your first user account.

### Configure LLM APIs (Optional)

In the admin dashboard, go to Settings → API Configuration and add your API keys:
- OpenAI API Key
- Perplexity API Key
- Google Gemini API Key

## Testing the Deployment

### Test Flask Application

```bash
curl https://pharmasight-app.onrender.com/health
```

Expected response:
```json
{
  "status": "healthy",
  "timestamp": "2024-12-24T...",
  "version": "4.0.0-combined"
}
```

### Test Admin Dashboard

```bash
curl https://pharmasight-admin.onrender.com/api/health
```

### Test Compound Analysis

```bash
curl -X POST https://pharmasight-app.onrender.com/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"compound": "aspirin"}'
```

## Troubleshooting

**Services won't start:**
- Check logs in Render dashboard
- Verify environment variables are set correctly
- Ensure DATABASE_URL and REDIS_URL are correct

**Database connection error:**
- Use **Internal Database URL** not External
- Check database is running in Render dashboard

**Admin dashboard can't reach Flask app:**
- Verify FLASK_API_URL is correct
- Use the Render service URL (e.g., `https://pharmasight-app.onrender.com`)
- Check CORS settings

## Next Steps

1. Import compound data
2. Configure autonomous research schedule
3. Set up user accounts and roles
4. Test all drug discovery features
5. Configure backups

## Support

- Documentation: See `COMBINED_DEPLOYMENT_GUIDE.md`
- Integration Report: See `COMBINED_INTEGRATION_REPORT.md`
- Admin Dashboard: See `admin-dashboard/README.md`

---

**You're ready to go!** 🚀
