# PharmaSight Platform - Render Deployment Guide

This guide will walk you through deploying the PharmaSight Combined platform to Render.

## Prerequisites

- Render account (sign up at https://render.com)
- GitHub account with access to pharmasight-platform repository
- Combined branch pushed to GitHub

## Deployment Options

### Option 1: Automatic Deployment with render.yaml (Recommended)

The Combined branch includes a `render.yaml` file for Infrastructure as Code deployment.

#### Steps:

1. **Go to Render Dashboard**
   - Visit https://dashboard.render.com

2. **Create New Blueprint**
   - Click "New" → "Blueprint"
   - Connect your GitHub account if not already connected
   - Select the `pharmasight-platform` repository
   - Select the `Combined` branch
   - Render will automatically detect `render.yaml`

3. **Review Services**
   Render will create:
   - PostgreSQL database (`pharmasight-db`)
   - Redis instance (`pharmasight-redis`)
   - Flask web service (`pharmasight-app`)
   - Admin dashboard web service (`pharmasight-admin`)

4. **Configure Environment Variables**
   - LLM API keys (optional, can be added later):
     - `OPENAI_API_KEY`
     - `PERPLEXITY_API_KEY`
     - `GEMINI_API_KEY`
     - `ANTHROPIC_API_KEY`

5. **Deploy**
   - Click "Apply" to create all services
   - Wait for deployment to complete (5-10 minutes)

6. **Access Your Applications**
   - Flask app: `https://pharmasight-app.onrender.com`
   - Admin dashboard: `https://pharmasight-admin.onrender.com`

### Option 2: Manual Deployment

If you prefer to create services manually:

#### Step 1: Create PostgreSQL Database

1. Go to Render Dashboard → New → PostgreSQL
2. Configure:
   - **Name:** `pharmasight-db`
   - **Database:** `pharmasight`
   - **User:** `pharmasight`
   - **Region:** Oregon (or your preferred region)
   - **Plan:** Starter ($7/month)
3. Click "Create Database"
4. **Copy the Internal Database URL** (you'll need this)

#### Step 2: Create Redis Instance

1. Go to Render Dashboard → New → Redis
2. Configure:
   - **Name:** `pharmasight-redis`
   - **Region:** Oregon (same as database)
   - **Plan:** Starter ($10/month)
   - **Maxmemory Policy:** allkeys-lru
3. Click "Create Redis"
4. **Copy the Internal Redis URL**

#### Step 3: Deploy Flask Application

1. Go to Render Dashboard → New → Web Service
2. Connect your GitHub repository
3. Configure:
   - **Name:** `pharmasight-app`
   - **Region:** Oregon
   - **Branch:** `Combined`
   - **Root Directory:** `/` (leave empty)
   - **Environment:** Docker
   - **Plan:** Starter ($7/month)
   - **Health Check Path:** `/health`

4. **Environment Variables:**
   ```
   DATABASE_URL=<paste-internal-database-url>
   REDIS_URL=<paste-internal-redis-url>
   FLASK_ENV=production
   PORT=5000
   SECRET_KEY=<generate-random-string>
   ```

   Optional LLM API keys:
   ```
   OPENAI_API_KEY=<your-key>
   PERPLEXITY_API_KEY=<your-key>
   GEMINI_API_KEY=<your-key>
   ANTHROPIC_API_KEY=<your-key>
   ```

5. Click "Create Web Service"

#### Step 4: Deploy Admin Dashboard

1. Go to Render Dashboard → New → Web Service
2. Connect your GitHub repository
3. Configure:
   - **Name:** `pharmasight-admin`
   - **Region:** Oregon
   - **Branch:** `Combined`
   - **Root Directory:** `admin-dashboard`
   - **Environment:** Docker
   - **Dockerfile Path:** `admin-dashboard/Dockerfile`
   - **Plan:** Starter ($7/month)
   - **Health Check Path:** `/api/health`

4. **Environment Variables:**
   ```
   DATABASE_URL=<same-database-url-as-flask>
   FLASK_API_URL=https://pharmasight-app.onrender.com
   JWT_SECRET=<generate-random-string>
   PORT=3000
   NODE_ENV=production
   ALLOWED_ORIGINS=https://pharmasight-app.onrender.com,https://pharmasight-admin.onrender.com
   ```

5. Click "Create Web Service"

## Post-Deployment Configuration

### 1. Configure Custom Domains (Optional)

#### For Flask Application:
1. Go to `pharmasight-app` service → Settings → Custom Domain
2. Add domain: `app.pharmasight.com`
3. Update your DNS records as instructed by Render:
   - Type: CNAME
   - Name: app
   - Value: pharmasight-app.onrender.com

#### For Admin Dashboard:
1. Go to `pharmasight-admin` service → Settings → Custom Domain
2. Add domain: `admin.pharmasight.com`
3. Update your DNS records:
   - Type: CNAME
   - Name: admin
   - Value: pharmasight-admin.onrender.com

**Note:** SSL certificates are automatically provisioned by Render.

### 2. Update CORS Settings

After setting up custom domains, update the `ALLOWED_ORIGINS` environment variable in the admin dashboard:

```
ALLOWED_ORIGINS=https://app.pharmasight.com,https://admin.pharmasight.com
```

### 3. Configure LLM API Keys

You can add LLM API keys in two ways:

**Option A: Via Render Dashboard**
1. Go to each service → Environment
2. Add environment variables for API keys
3. Save changes (service will redeploy)

**Option B: Via Admin Dashboard UI**
1. Access admin dashboard
2. Go to Settings → API Configuration
3. Enter API keys through the UI
4. Keys are stored in database

### 4. Create Admin User

1. Access the admin dashboard
2. Click "Sign Up" or use the registration endpoint
3. Create your first admin account
4. Set role to "admin" in the database if needed

## Verification

### Test Flask Application

```bash
# Health check
curl https://pharmasight-app.onrender.com/health

# Expected response:
# {"status": "healthy", "timestamp": "...", "version": "4.0.0-combined"}

# Test compound analysis
curl -X POST https://pharmasight-app.onrender.com/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"compound": "aspirin"}'
```

### Test Admin Dashboard

```bash
# Health check
curl https://pharmasight-admin.onrender.com/api/health

# Expected response:
# {"status": "ok", "timestamp": "..."}
```

### Test Integration

1. Log in to admin dashboard
2. Navigate to compound management
3. Create or view a compound
4. Verify data appears in Flask application
5. Test data export features

## Monitoring

### View Logs

**Flask Application:**
1. Go to `pharmasight-app` service
2. Click "Logs" tab
3. Monitor application logs in real-time

**Admin Dashboard:**
1. Go to `pharmasight-admin` service
2. Click "Logs" tab
3. Monitor application logs

### Metrics

Render provides built-in metrics:
- CPU usage
- Memory usage
- Request count
- Response times

Access metrics from each service's "Metrics" tab.

### Health Checks

Render automatically monitors health check endpoints:
- Flask: `/health`
- Admin: `/api/health`

If health checks fail, Render will attempt to restart the service.

## Troubleshooting

### Flask App Won't Start

**Check logs for errors:**
```bash
# Common issues:
# - Missing environment variables
# - Database connection failed
# - RDKit import errors
```

**Solutions:**
1. Verify all environment variables are set
2. Check DATABASE_URL is the Internal URL (not External)
3. Ensure Dockerfile includes all dependencies
4. Review build logs for Python package installation errors

### Admin Dashboard Won't Start

**Check logs for errors:**
```bash
# Common issues:
# - Node.js build failures
# - Missing dependencies
# - Database connection failed
```

**Solutions:**
1. Verify `pnpm-lock.yaml` is committed
2. Check DATABASE_URL matches Flask app
3. Verify FLASK_API_URL is correct
4. Ensure Node.js version is 18+

### Database Connection Issues

**Symptoms:**
- Services fail to start
- "Connection refused" errors
- Timeout errors

**Solutions:**
1. Use **Internal Database URL** not External
2. Verify database is running (check Render dashboard)
3. Check database credentials
4. Ensure services are in the same region

### CORS Errors

**Symptoms:**
- Admin dashboard can't reach Flask API
- Browser console shows CORS errors

**Solutions:**
1. Update `ALLOWED_ORIGINS` in admin dashboard
2. Include both Render URLs and custom domains
3. Restart admin dashboard after updating

### Health Check Failures

**Symptoms:**
- Service shows as unhealthy
- Automatic restarts

**Solutions:**
1. Verify health check paths are correct
2. Test health endpoints manually
3. Check application is listening on correct port
4. Review application logs for startup errors

## Scaling

### Vertical Scaling (Upgrade Plan)

Upgrade to higher-tier plans for more resources:
- **Starter:** 512 MB RAM, 0.5 CPU
- **Standard:** 2 GB RAM, 1 CPU
- **Pro:** 4 GB RAM, 2 CPU

### Horizontal Scaling

For high-traffic scenarios:
1. Enable auto-scaling in service settings
2. Set min/max instance counts
3. Configure scaling triggers

### Database Scaling

1. Upgrade PostgreSQL plan for more storage/connections
2. Enable connection pooling
3. Consider read replicas for read-heavy workloads

## Cost Estimation

### Minimum Configuration (Starter Plans)

| Service | Plan | Cost/Month |
|---------|------|------------|
| PostgreSQL | Starter | $7 |
| Redis | Starter | $10 |
| Flask App | Starter | $7 |
| Admin Dashboard | Starter | $7 |
| **Total** | | **$31/month** |

### Recommended Configuration (Standard Plans)

| Service | Plan | Cost/Month |
|---------|------|------------|
| PostgreSQL | Standard | $20 |
| Redis | Standard | $25 |
| Flask App | Standard | $25 |
| Admin Dashboard | Standard | $25 |
| **Total** | | **$95/month** |

**Note:** First month may be free with Render's trial credits.

## Maintenance

### Updating the Application

1. Push changes to Combined branch on GitHub
2. Render automatically detects changes
3. Services redeploy automatically
4. Monitor deployment logs

### Database Backups

Render automatically backs up PostgreSQL:
- **Starter:** Daily backups, 7-day retention
- **Standard:** Daily backups, 14-day retention
- **Pro:** Daily backups, 30-day retention

**Manual backup:**
```bash
# Download backup from Render dashboard
# Or use pg_dump via shell access
```

### Monitoring and Alerts

Set up alerts in Render dashboard:
1. Go to service → Settings → Notifications
2. Configure email/Slack notifications
3. Set up alerts for:
   - Service failures
   - High CPU/memory usage
   - Health check failures

## Security Best Practices

1. **Use Environment Variables**
   - Never commit secrets to Git
   - Use Render's environment variable management

2. **Enable HTTPS**
   - Render provides free SSL certificates
   - Enforce HTTPS in application

3. **Restrict Database Access**
   - Use Internal Database URL
   - Don't expose database publicly

4. **Regular Updates**
   - Keep dependencies updated
   - Monitor security advisories

5. **Access Control**
   - Use strong JWT secrets
   - Implement RBAC in admin dashboard
   - Regular audit of user permissions

## Support

### Render Support
- Documentation: https://render.com/docs
- Community: https://community.render.com
- Support: support@render.com

### PharmaSight Support
- Check logs first
- Review troubleshooting section
- Consult COMBINED_DEPLOYMENT_GUIDE.md

## Next Steps

After successful deployment:

1. ✅ Verify both applications are running
2. ✅ Configure custom domains
3. ✅ Add LLM API keys
4. ✅ Create admin user accounts
5. ✅ Import compound data
6. ✅ Test all features
7. ✅ Set up monitoring and alerts
8. ✅ Configure automated backups

---

**Deployment complete!** Your PharmaSight Platform is now live on Render. 🚀
