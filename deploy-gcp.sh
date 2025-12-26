#!/bin/bash
# PharmaSight Platform - Google Cloud Run Deployment Script
# This script deploys both Flask and Admin applications to Google Cloud Run

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}PharmaSight Platform - GCP Deployment${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    echo -e "${RED}Error: gcloud CLI is not installed${NC}"
    echo "Please install it from: https://cloud.google.com/sdk/docs/install"
    exit 1
fi

# Get project ID
if [ -z "$GCP_PROJECT_ID" ]; then
    echo -e "${YELLOW}Enter your GCP Project ID:${NC}"
    read -r GCP_PROJECT_ID
fi

echo -e "${GREEN}Using GCP Project: $GCP_PROJECT_ID${NC}"
echo ""

# Set project
gcloud config set project "$GCP_PROJECT_ID"

# Enable required APIs
echo -e "${BLUE}[1/8] Enabling required GCP APIs...${NC}"
gcloud services enable \
    cloudbuild.googleapis.com \
    run.googleapis.com \
    sql-component.googleapis.com \
    sqladmin.googleapis.com \
    redis.googleapis.com \
    containerregistry.googleapis.com \
    secretmanager.googleapis.com

echo -e "${GREEN}✓ APIs enabled${NC}"
echo ""

# Create Cloud SQL instance (PostgreSQL)
echo -e "${BLUE}[2/8] Creating Cloud SQL (PostgreSQL) instance...${NC}"
if gcloud sql instances describe pharmasight-db --project="$GCP_PROJECT_ID" &> /dev/null; then
    echo -e "${YELLOW}Cloud SQL instance already exists${NC}"
else
    gcloud sql instances create pharmasight-db \
        --database-version=POSTGRES_15 \
        --tier=db-f1-micro \
        --region=us-central1 \
        --root-password="$(openssl rand -base64 32)" \
        --backup \
        --backup-start-time=03:00
    
    # Create database
    gcloud sql databases create pharmasight \
        --instance=pharmasight-db
    
    echo -e "${GREEN}✓ Cloud SQL instance created${NC}"
fi
echo ""

# Create Memorystore (Redis) instance
echo -e "${BLUE}[3/8] Creating Memorystore (Redis) instance...${NC}"
if gcloud redis instances describe pharmasight-redis --region=us-central1 --project="$GCP_PROJECT_ID" &> /dev/null; then
    echo -e "${YELLOW}Redis instance already exists${NC}"
else
    gcloud redis instances create pharmasight-redis \
        --size=1 \
        --region=us-central1 \
        --redis-version=redis_7_0 \
        --tier=basic
    
    echo -e "${GREEN}✓ Redis instance created${NC}"
fi
echo ""

# Build and push Docker images
echo -e "${BLUE}[4/8] Building Flask application Docker image...${NC}"
docker build -t gcr.io/"$GCP_PROJECT_ID"/pharmasight-flask:latest -f Dockerfile .
docker push gcr.io/"$GCP_PROJECT_ID"/pharmasight-flask:latest
echo -e "${GREEN}✓ Flask image built and pushed${NC}"
echo ""

echo -e "${BLUE}[5/8] Building Admin dashboard Docker image...${NC}"
docker build -t gcr.io/"$GCP_PROJECT_ID"/pharmasight-admin:latest -f admin-dashboard/Dockerfile ./admin-dashboard
docker push gcr.io/"$GCP_PROJECT_ID"/pharmasight-admin:latest
echo -e "${GREEN}✓ Admin image built and pushed${NC}"
echo ""

# Get database connection string
echo -e "${BLUE}[6/8] Getting database connection details...${NC}"
DB_INSTANCE_CONNECTION_NAME=$(gcloud sql instances describe pharmasight-db --format='value(connectionName)')
echo -e "${GREEN}✓ Database connection name: $DB_INSTANCE_CONNECTION_NAME${NC}"
echo ""

# Deploy Flask application
echo -e "${BLUE}[7/8] Deploying Flask application to Cloud Run...${NC}"
gcloud run deploy pharmasight-app \
    --image=gcr.io/"$GCP_PROJECT_ID"/pharmasight-flask:latest \
    --platform=managed \
    --region=us-central1 \
    --allow-unauthenticated \
    --memory=2Gi \
    --cpu=2 \
    --timeout=300 \
    --max-instances=10 \
    --add-cloudsql-instances="$DB_INSTANCE_CONNECTION_NAME" \
    --set-env-vars="FLASK_ENV=production,PORT=5000" \
    --set-secrets="DATABASE_URL=pharmasight-db-url:latest,REDIS_URL=pharmasight-redis-url:latest"

FLASK_URL=$(gcloud run services describe pharmasight-app --region=us-central1 --format='value(status.url)')
echo -e "${GREEN}✓ Flask app deployed at: $FLASK_URL${NC}"
echo ""

# Deploy Admin dashboard
echo -e "${BLUE}[8/8] Deploying Admin dashboard to Cloud Run...${NC}"
gcloud run deploy pharmasight-admin \
    --image=gcr.io/"$GCP_PROJECT_ID"/pharmasight-admin:latest \
    --platform=managed \
    --region=us-central1 \
    --allow-unauthenticated \
    --memory=1Gi \
    --cpu=1 \
    --timeout=60 \
    --max-instances=5 \
    --add-cloudsql-instances="$DB_INSTANCE_CONNECTION_NAME" \
    --set-env-vars="NODE_ENV=production,PORT=3000,FLASK_API_URL=$FLASK_URL" \
    --set-secrets="DATABASE_URL=pharmasight-db-url:latest,JWT_SECRET=pharmasight-jwt-secret:latest"

ADMIN_URL=$(gcloud run services describe pharmasight-admin --region=us-central1 --format='value(status.url)')
echo -e "${GREEN}✓ Admin dashboard deployed at: $ADMIN_URL${NC}"
echo ""

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Deployment Complete!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${GREEN}Flask Application:${NC}"
echo -e "  URL: $FLASK_URL"
echo -e "  Health: $FLASK_URL/health"
echo ""
echo -e "${GREEN}Admin Dashboard:${NC}"
echo -e "  URL: $ADMIN_URL"
echo -e "  Health: $ADMIN_URL/api/health"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "1. Configure custom domains (optional)"
echo "2. Add LLM API keys via Secret Manager"
echo "3. Create admin user account"
echo "4. Import compound data"
echo ""
echo -e "${BLUE}To view logs:${NC}"
echo "  gcloud run services logs read pharmasight-app --region=us-central1"
echo "  gcloud run services logs read pharmasight-admin --region=us-central1"
echo ""
