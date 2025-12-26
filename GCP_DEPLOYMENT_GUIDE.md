# PharmaSight Platform - Google Cloud Run Deployment Guide

Complete guide for deploying PharmaSight Platform to Google Cloud Platform using Cloud Run.

---

## Overview

This guide will help you deploy the PharmaSight Combined platform to Google Cloud Platform (GCP) using Cloud Run, which provides automatic scaling, managed infrastructure, and enterprise-grade reliability.

### What You'll Deploy

**Architecture:**
- **Flask Application** (Drug Discovery Platform) → Cloud Run
- **Admin Dashboard** (Management Interface) → Cloud Run
- **PostgreSQL Database** → Cloud SQL
- **Redis Cache** → Memorystore
- **Container Registry** → Google Container Registry (GCR)

**Estimated Cost:** $25-75/month (depending on usage)

---

## Prerequisites

### 1. Google Cloud Platform Account

- Sign up at https://cloud.google.com
- Enable billing (required for Cloud Run)
- **New users get $300 free credit** (valid for 90 days)

### 2. Install Google Cloud SDK

**macOS:**
```bash
brew install google-cloud-sdk
```

**Linux:**
```bash
curl https://sdk.cloud.google.com | bash
exec -l $SHELL
```

**Windows:**
Download from https://cloud.google.com/sdk/docs/install

### 3. Install Docker

- Download from https://www.docker.com/products/docker-desktop
- Ensure Docker is running before deployment

### 4. Authenticate with GCP

```bash
gcloud auth login
gcloud auth configure-docker
```

---

## Deployment Methods

### Method 1: Automated Deployment (Recommended)

Use the provided deployment script for one-command deployment.

#### Steps:

1. **Clone the Repository**
```bash
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
git checkout Combined
```

2. **Set Your GCP Project ID**
```bash
export GCP_PROJECT_ID="your-project-id"
```

3. **Run Deployment Script**
```bash
./deploy-gcp.sh
```

The script will:
- Enable required GCP APIs
- Create Cloud SQL (PostgreSQL) instance
- Create Memorystore (Redis) instance
- Build and push Docker images
- Deploy both applications to Cloud Run
- Configure environment variables

**Deployment Time:** 15-20 minutes

---

### Method 2: Manual Deployment

For more control over the deployment process.

#### Step 1: Enable Required APIs

```bash
gcloud services enable \
    cloudbuild.googleapis.com \
    run.googleapis.com \
    sql-component.googleapis.com \
    sqladmin.googleapis.com \
    redis.googleapis.com \
    containerregistry.googleapis.com \
    secretmanager.googleapis.com
```

#### Step 2: Create Cloud SQL Instance

```bash
# Create PostgreSQL instance
gcloud sql instances create pharmasight-db \
    --database-version=POSTGRES_15 \
    --tier=db-f1-micro \
    --region=us-central1 \
    --root-password="YOUR_SECURE_PASSWORD" \
    --backup \
    --backup-start-time=03:00

# Create database
gcloud sql databases create pharmasight \
    --instance=pharmasight-db

# Create user
gcloud sql users create pharmasight \
    --instance=pharmasight-db \
    --password="YOUR_DB_PASSWORD"
```

#### Step 3: Create Memorystore (Redis) Instance

```bash
gcloud redis instances create pharmasight-redis \
    --size=1 \
    --region=us-central1 \
    --redis-version=redis_7_0 \
    --tier=basic
```

#### Step 4: Build and Push Docker Images

```bash
# Set project ID
export GCP_PROJECT_ID="your-project-id"

# Build Flask application
docker build -t gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest \
    -f Dockerfile .

# Push Flask image
docker push gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest

# Build Admin dashboard
docker build -t gcr.io/$GCP_PROJECT_ID/pharmasight-admin:latest \
    -f admin-dashboard/Dockerfile ./admin-dashboard

# Push Admin image
docker push gcr.io/$GCP_PROJECT_ID/pharmasight-admin:latest
```

#### Step 5: Create Secrets

```bash
# Database URL
echo -n "postgresql://pharmasight:YOUR_DB_PASSWORD@/pharmasight?host=/cloudsql/YOUR_CONNECTION_NAME" | \
    gcloud secrets create pharmasight-db-url --data-file=-

# Redis URL
REDIS_HOST=$(gcloud redis instances describe pharmasight-redis --region=us-central1 --format='value(host)')
echo -n "redis://$REDIS_HOST:6379" | \
    gcloud secrets create pharmasight-redis-url --data-file=-

# JWT Secret
openssl rand -base64 32 | gcloud secrets create pharmasight-jwt-secret --data-file=-
```

#### Step 6: Deploy Flask Application

```bash
# Get Cloud SQL connection name
DB_CONNECTION=$(gcloud sql instances describe pharmasight-db --format='value(connectionName)')

# Deploy Flask app
gcloud run deploy pharmasight-app \
    --image=gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest \
    --platform=managed \
    --region=us-central1 \
    --allow-unauthenticated \
    --memory=2Gi \
    --cpu=2 \
    --timeout=300 \
    --max-instances=10 \
    --add-cloudsql-instances=$DB_CONNECTION \
    --set-env-vars="FLASK_ENV=production,PORT=5000" \
    --set-secrets="DATABASE_URL=pharmasight-db-url:latest,REDIS_URL=pharmasight-redis-url:latest"
```

#### Step 7: Deploy Admin Dashboard

```bash
# Get Flask app URL
FLASK_URL=$(gcloud run services describe pharmasight-app --region=us-central1 --format='value(status.url)')

# Deploy Admin dashboard
gcloud run deploy pharmasight-admin \
    --image=gcr.io/$GCP_PROJECT_ID/pharmasight-admin:latest \
    --platform=managed \
    --region=us-central1 \
    --allow-unauthenticated \
    --memory=1Gi \
    --cpu=1 \
    --timeout=60 \
    --max-instances=5 \
    --add-cloudsql-instances=$DB_CONNECTION \
    --set-env-vars="NODE_ENV=production,PORT=3000,FLASK_API_URL=$FLASK_URL" \
    --set-secrets="DATABASE_URL=pharmasight-db-url:latest,JWT_SECRET=pharmasight-jwt-secret:latest"
```

---

### Method 3: Cloud Build (CI/CD)

Use Cloud Build for automated deployments from GitHub.

#### Setup:

1. **Connect GitHub Repository**
```bash
gcloud alpha builds triggers create github \
    --repo-name=pharmasight-platform \
    --repo-owner=justincihi \
    --branch-pattern="^Combined$" \
    --build-config=cloudbuild.yaml
```

2. **Push to GitHub**
```bash
git push origin Combined
```

Cloud Build will automatically build and deploy both applications.

---

## Post-Deployment Configuration

### 1. Get Application URLs

```bash
# Flask application URL
gcloud run services describe pharmasight-app --region=us-central1 --format='value(status.url)'

# Admin dashboard URL
gcloud run services describe pharmasight-admin --region=us-central1 --format='value(status.url)'
```

### 2. Configure Custom Domains (Optional)

#### For Flask Application:

```bash
gcloud run domain-mappings create \
    --service=pharmasight-app \
    --domain=app.pharmasight.com \
    --region=us-central1
```

#### For Admin Dashboard:

```bash
gcloud run domain-mappings create \
    --service=pharmasight-admin \
    --domain=admin.pharmasight.com \
    --region=us-central1
```

Update your DNS records with the provided values.

### 3. Add LLM API Keys

```bash
# OpenAI
echo -n "your-openai-key" | gcloud secrets create openai-api-key --data-file=-

# Perplexity
echo -n "your-perplexity-key" | gcloud secrets create perplexity-api-key --data-file=-

# Gemini
echo -n "your-gemini-key" | gcloud secrets create gemini-api-key --data-file=-

# Anthropic
echo -n "your-anthropic-key" | gcloud secrets create anthropic-api-key --data-file=-

# Update Flask service
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --update-secrets=OPENAI_API_KEY=openai-api-key:latest,PERPLEXITY_API_KEY=perplexity-api-key:latest,GEMINI_API_KEY=gemini-api-key:latest,ANTHROPIC_API_KEY=anthropic-api-key:latest
```

### 4. Create Admin User

Access the admin dashboard and create your first user account.

---

## Verification

### Test Flask Application

```bash
FLASK_URL=$(gcloud run services describe pharmasight-app --region=us-central1 --format='value(status.url)')

# Health check
curl $FLASK_URL/health

# Test compound analysis
curl -X POST $FLASK_URL/api/analyze \
    -H "Content-Type: application/json" \
    -d '{"compound": "aspirin"}'
```

### Test Admin Dashboard

```bash
ADMIN_URL=$(gcloud run services describe pharmasight-admin --region=us-central1 --format='value(status.url)')

# Health check
curl $ADMIN_URL/api/health
```

---

## Monitoring and Logging

### View Logs

```bash
# Flask application logs
gcloud run services logs read pharmasight-app --region=us-central1 --limit=50

# Admin dashboard logs
gcloud run services logs read pharmasight-admin --region=us-central1 --limit=50

# Follow logs in real-time
gcloud run services logs tail pharmasight-app --region=us-central1
```

### View Metrics

Access Cloud Console:
1. Go to https://console.cloud.google.com
2. Navigate to Cloud Run
3. Select your service
4. Click "Metrics" tab

Metrics include:
- Request count
- Request latency
- Container CPU utilization
- Container memory utilization
- Billable container time

### Set Up Alerts

```bash
# Create alert for high error rate
gcloud alpha monitoring policies create \
    --notification-channels=YOUR_CHANNEL_ID \
    --display-name="PharmaSight High Error Rate" \
    --condition-display-name="Error rate > 5%" \
    --condition-threshold-value=0.05 \
    --condition-threshold-duration=300s
```

---

## Scaling Configuration

### Automatic Scaling

Cloud Run automatically scales based on traffic. Configure limits:

```bash
# Update Flask app scaling
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --min-instances=0 \
    --max-instances=20 \
    --concurrency=80

# Update Admin dashboard scaling
gcloud run services update pharmasight-admin \
    --region=us-central1 \
    --min-instances=0 \
    --max-instances=10 \
    --concurrency=80
```

### Resource Allocation

Adjust CPU and memory based on usage:

```bash
# Upgrade Flask app resources
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --memory=4Gi \
    --cpu=4

# Upgrade Admin dashboard resources
gcloud run services update pharmasight-admin \
    --region=us-central1 \
    --memory=2Gi \
    --cpu=2
```

---

## Cost Optimization

### Estimated Monthly Costs

| Service | Configuration | Estimated Cost |
|---------|--------------|----------------|
| Cloud Run (Flask) | 2 vCPU, 2GB RAM | $10-30 |
| Cloud Run (Admin) | 1 vCPU, 1GB RAM | $5-15 |
| Cloud SQL | db-f1-micro | $10-15 |
| Memorystore | 1GB Basic | $10-15 |
| Container Registry | 10GB storage | $1 |
| **Total** | | **$36-76/month** |

### Cost Reduction Tips

1. **Use Minimum Instances = 0**
   - Scales to zero when not in use
   - Saves money during low traffic periods

2. **Optimize Container Size**
   - Remove unnecessary dependencies
   - Use multi-stage Docker builds

3. **Use Cloud SQL Proxy**
   - Reduces connection overhead
   - Improves performance

4. **Enable Request Compression**
   - Reduces bandwidth costs
   - Improves response times

5. **Set Appropriate Timeouts**
   - Prevents long-running requests
   - Reduces billable time

---

## Troubleshooting

### Common Issues

#### 1. Deployment Fails

**Error:** "Permission denied"

**Solution:**
```bash
# Grant necessary permissions
gcloud projects add-iam-policy-binding $GCP_PROJECT_ID \
    --member="serviceAccount:YOUR_SERVICE_ACCOUNT" \
    --role="roles/run.admin"
```

#### 2. Database Connection Fails

**Error:** "Could not connect to Cloud SQL"

**Solution:**
- Verify Cloud SQL instance is running
- Check connection name is correct
- Ensure Cloud SQL Admin API is enabled
- Verify secrets are created correctly

#### 3. Container Build Fails

**Error:** "Error building Docker image"

**Solution:**
- Check Dockerfile syntax
- Verify all dependencies are listed
- Increase Cloud Build timeout
- Check build logs for specific errors

#### 4. Out of Memory Errors

**Error:** "Container exceeded memory limit"

**Solution:**
```bash
# Increase memory allocation
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --memory=4Gi
```

#### 5. Timeout Errors

**Error:** "Request timeout"

**Solution:**
```bash
# Increase timeout
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --timeout=600
```

---

## Maintenance

### Update Applications

```bash
# Rebuild and redeploy
./deploy-gcp.sh
```

Or manually:

```bash
# Build new images
docker build -t gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest -f Dockerfile .
docker push gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest

# Deploy update
gcloud run deploy pharmasight-app \
    --image=gcr.io/$GCP_PROJECT_ID/pharmasight-flask:latest \
    --region=us-central1
```

### Database Backups

Cloud SQL automatically backs up your database daily. To create manual backup:

```bash
gcloud sql backups create \
    --instance=pharmasight-db \
    --description="Manual backup before update"
```

### Restore from Backup

```bash
# List backups
gcloud sql backups list --instance=pharmasight-db

# Restore
gcloud sql backups restore BACKUP_ID \
    --backup-instance=pharmasight-db \
    --backup-id=BACKUP_ID
```

---

## Security Best Practices

### 1. Use Secret Manager

Never hardcode secrets. Always use Secret Manager:

```bash
echo -n "your-secret" | gcloud secrets create secret-name --data-file=-
```

### 2. Enable VPC Connector (Optional)

For enhanced security, use VPC connector:

```bash
gcloud compute networks vpc-access connectors create pharmasight-connector \
    --region=us-central1 \
    --range=10.8.0.0/28

gcloud run services update pharmasight-app \
    --region=us-central1 \
    --vpc-connector=pharmasight-connector
```

### 3. Restrict Access

For production, consider requiring authentication:

```bash
gcloud run services update pharmasight-app \
    --region=us-central1 \
    --no-allow-unauthenticated
```

### 4. Enable Cloud Armor

Protect against DDoS attacks:

```bash
gcloud compute security-policies create pharmasight-policy \
    --description="DDoS protection for PharmaSight"
```

---

## Support and Resources

### GCP Documentation
- Cloud Run: https://cloud.google.com/run/docs
- Cloud SQL: https://cloud.google.com/sql/docs
- Memorystore: https://cloud.google.com/memorystore/docs

### PharmaSight Documentation
- Main README: `README.md`
- Quick Start: `QUICKSTART.md`
- Deployment Guide: `COMBINED_DEPLOYMENT_GUIDE.md`

### Getting Help
- GCP Support: https://cloud.google.com/support
- Community: https://stackoverflow.com/questions/tagged/google-cloud-run
- GitHub Issues: https://github.com/justincihi/pharmasight-platform/issues

---

## Summary

You now have PharmaSight Platform deployed on Google Cloud Run with:

✅ Automatic scaling and load balancing  
✅ Managed PostgreSQL database  
✅ Redis caching for performance  
✅ Enterprise-grade security  
✅ Automatic SSL certificates  
✅ Built-in monitoring and logging  
✅ Cost-effective pay-per-use pricing  

**Your applications are live and ready to use!**

---

**Deployment Guide Version:** 1.0  
**Last Updated:** December 26, 2024  
**Author:** Manus AI
