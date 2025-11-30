# PharmaSight™ Microservices Deployment Guide

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Quick Start](#quick-start)
3. [Configuration](#configuration)
4. [Service Architecture](#service-architecture)
5. [Development Workflow](#development-workflow)
6. [Production Deployment](#production-deployment)
7. [Monitoring & Debugging](#monitoring--debugging)
8. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Software
- **Docker** (version 20.10+)
- **Docker Compose** (version 1.29+)

### Verify Installation
```bash
docker --version
docker-compose --version
```

### System Requirements
- **Minimum:** 4GB RAM, 2 CPU cores, 10GB disk space
- **Recommended:** 8GB RAM, 4 CPU cores, 20GB disk space
- **For Quantum Calculations:** 16GB RAM, 8 CPU cores

---

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
```

### 2. Configure Environment
```bash
# Copy example environment file
cp .env.example .env

# Edit .env with your preferred settings (optional)
nano .env
```

### 3. Build and Start Services
```bash
# Build all images and start services in detached mode
docker-compose up --build -d

# View startup logs
docker-compose logs -f
```

### 4. Verify Deployment
```bash
# Check all services are healthy
curl http://localhost:8080/health

# Should return:
# {"gateway": "healthy", "services": {...}, "overall_status": "healthy"}
```

### 5. Access Services
- **API Gateway:** http://localhost:8080
- **API Documentation:** http://localhost:8080/docs
- **Compound Analysis:** http://localhost:8001/docs
- **Analog Generation:** http://localhost:8002/docs
- **ML Models:** http://localhost:8003/docs
- **Quantum Calculator:** http://localhost:8004/docs
- **Auth Service:** http://localhost:8005/docs

---

## Configuration

### Environment Variables

Create a `.env` file in the project root:

```bash
# Database Configuration
POSTGRES_USER=pharmasight_user
POSTGRES_PASSWORD=strong_password_here
POSTGRES_DB=pharmasight_db

# Authentication Secret (CHANGE THIS!)
SECRET_KEY=your-secret-key-minimum-32-characters-long

# Service Configuration (optional)
LOG_LEVEL=INFO
REQUEST_TIMEOUT=30
```

**Security Note:** Always use strong passwords and secrets in production!

### Port Configuration

Default ports defined in `docker-compose.yml`:
- **8080** - API Gateway
- **8001** - Compound Analysis
- **8002** - Analog Generation
- **8003** - ML Models
- **8004** - Quantum Calculator
- **8005** - Authentication
- **5432** - PostgreSQL
- **6379** - Redis

To change ports, edit `docker-compose.yml` or use environment variables.

---

## Service Architecture

### Service Dependency Graph
```
                    ┌─────────────┐
                    │ API Gateway │ :8080
                    └──────┬──────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
   ┌────▼────┐      ┌──────▼──────┐   ┌──────▼──────┐
   │Compound │      │   Analog    │   │  ML Models  │
   │Analysis │      │ Generation  │   │             │
   └────┬────┘      └──────┬──────┘   └──────┬──────┘
        │                  │                  │
   ┌────▼────┐      ┌──────▼──────┐   ┌──────▼──────┐
   │Quantum  │      │    Auth     │   │   Redis     │
   │Calc     │      │  Service    │   │   Cache     │
   └─────────┘      └─────────────┘   └──────┬──────┘
                                              │
                                       ┌──────▼──────┐
                                       │ PostgreSQL  │
                                       │  Database   │
                                       └─────────────┘
```

### Service Responsibilities

| Service | Purpose | Dependencies |
|---------|---------|--------------|
| **API Gateway** | Routes requests, load balancing | All services |
| **Compound Analysis** | RDKit molecular analysis | Redis, PostgreSQL |
| **Analog Generation** | Generate chemical analogs | Redis |
| **ML Models** | ADMET predictions | None |
| **Quantum Calculator** | DFT calculations | None |
| **Auth Service** | JWT authentication, RBAC | None |
| **PostgreSQL** | Persistent data storage | None |
| **Redis** | Caching, session storage | None |

---

## Development Workflow

### Running Individual Services

Start only specific services for focused development:

```bash
# Start just the compound analysis service with dependencies
docker-compose up -d db redis compound-service

# View logs
docker-compose logs -f compound-service
```

### Making Code Changes

1. **Edit service code** in `services/[service-name]/main.py`
2. **Rebuild the service:**
   ```bash
   docker-compose up -d --build compound-service
   ```
3. **View logs:**
   ```bash
   docker-compose logs -f compound-service
   ```

### Testing Services

#### Unit Tests (if available)
```bash
# Run tests in container
docker-compose exec compound-service python -m pytest
```

#### Manual API Testing
```bash
# Get authentication token
TOKEN=$(curl -X POST http://localhost:8080/auth-service/token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=researcher&password=password" | jq -r '.access_token')

# Test compound analysis
curl -X POST http://localhost:8080/compound-analysis/analyze \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"compound": "psilocybin"}'
```

### Adding New Services

1. **Create service directory:**
   ```bash
   mkdir -p services/new-service
   cd services/new-service
   ```

2. **Create files:**
   - `main.py` - FastAPI application
   - `requirements.txt` - Python dependencies
   - `Dockerfile` - Container definition

3. **Add to docker-compose.yml:**
   ```yaml
   new-service:
     build: ./services/new-service
     ports:
       - "8006:8000"
     restart: unless-stopped
     healthcheck:
       test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
       interval: 30s
       timeout: 10s
       retries: 3
     networks:
       - pharmasight-network
   ```

4. **Update API Gateway** to route to new service in `services/api-gateway/main.py`

---

## Production Deployment

### Security Checklist

- [ ] Change all default passwords
- [ ] Use strong SECRET_KEY (minimum 32 characters)
- [ ] Enable HTTPS/TLS
- [ ] Implement rate limiting
- [ ] Configure firewall rules
- [ ] Enable service-to-service authentication
- [ ] Regular security updates
- [ ] Backup strategy for database

### Performance Optimization

#### 1. Resource Limits
Add to docker-compose.yml:
```yaml
services:
  compound-service:
    deploy:
      resources:
        limits:
          cpus: '2.0'
          memory: 2G
        reservations:
          cpus: '1.0'
          memory: 1G
```

#### 2. Enable GPU (for ML service)
Uncomment GPU section in docker-compose.yml:
```yaml
ml-service:
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            count: 1
            capabilities: [gpu]
```

#### 3. Redis Persistence
Add to docker-compose.yml:
```yaml
redis:
  volumes:
    - redis_data:/data
  command: redis-server --appendonly yes
```

### Scaling Services

#### Horizontal Scaling
```bash
# Scale compound service to 3 instances
docker-compose up -d --scale compound-service=3
```

**Note:** Requires load balancer configuration (e.g., nginx)

#### Load Balancing with Nginx
Create `nginx.conf`:
```nginx
upstream compound_backend {
    server compound-service-1:8000;
    server compound-service-2:8000;
    server compound-service-3:8000;
}

server {
    listen 80;
    location /compound-analysis/ {
        proxy_pass http://compound_backend/;
    }
}
```

---

## Monitoring & Debugging

### View Logs

```bash
# All services
docker-compose logs -f

# Specific service
docker-compose logs -f compound-service

# Last 100 lines
docker-compose logs --tail=100 compound-service

# Follow logs with timestamps
docker-compose logs -f -t compound-service
```

### Health Checks

```bash
# Gateway health (includes all services)
curl http://localhost:8080/health | jq

# Individual service health
curl http://localhost:8001/health
curl http://localhost:8002/health
curl http://localhost:8003/health
curl http://localhost:8004/health
curl http://localhost:8005/health
```

### Service Status

```bash
# List all services
docker-compose ps

# Inspect specific service
docker-compose ps compound-service

# View resource usage
docker stats
```

### Debugging Inside Containers

```bash
# Execute bash in running container
docker-compose exec compound-service /bin/bash

# Run Python REPL
docker-compose exec compound-service python

# Check installed packages
docker-compose exec compound-service pip list
```

### Database Access

```bash
# Connect to PostgreSQL
docker-compose exec db psql -U pharmasight_user -d pharmasight_db

# Common queries
\dt              # List tables
\d+ compounds    # Describe compounds table
SELECT * FROM compounds LIMIT 10;
```

### Redis Debugging

```bash
# Connect to Redis CLI
docker-compose exec redis redis-cli

# Common commands
KEYS *           # List all keys
GET key_name     # Get value
FLUSHALL         # Clear all data (careful!)
INFO             # Server info
```

---

## Troubleshooting

### Common Issues

#### 1. Services Won't Start
```bash
# Check logs for errors
docker-compose logs

# Remove old containers and volumes
docker-compose down -v
docker-compose up --build -d
```

#### 2. Port Already in Use
```bash
# Find process using port 8080
lsof -i :8080

# Kill process
kill -9 <PID>

# Or change port in docker-compose.yml
```

#### 3. Database Connection Errors
```bash
# Check database is running
docker-compose ps db

# Verify environment variables
docker-compose exec compound-service env | grep DATABASE_URL

# Reset database
docker-compose down -v
docker-compose up -d db
```

#### 4. RDKit Import Errors
```bash
# Verify RDKit installation
docker-compose exec compound-service python -c "from rdkit import Chem; print(Chem.rdBase.rdkitVersion)"

# Rebuild with --no-cache
docker-compose build --no-cache compound-service
```

#### 5. Health Check Failures
```bash
# Check if curl is installed in container
docker-compose exec compound-service which curl

# Test health endpoint directly
docker-compose exec compound-service curl http://localhost:8000/health
```

#### 6. Authentication Issues
```bash
# Test token generation
curl -X POST http://localhost:8080/auth-service/token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=researcher&password=password"

# Verify SECRET_KEY is set
docker-compose exec auth-service env | grep SECRET_KEY
```

### Performance Issues

#### High Memory Usage
```bash
# Check container resource usage
docker stats

# Limit container memory
# Add to docker-compose.yml:
mem_limit: 2g
```

#### Slow Response Times
```bash
# Check Redis cache
docker-compose exec redis redis-cli DBSIZE

# Clear cache if needed
docker-compose exec redis redis-cli FLUSHALL

# Check database connections
docker-compose exec db psql -U pharmasight_user -d pharmasight_db -c "SELECT count(*) FROM pg_stat_activity;"
```

### Clean Slate

```bash
# Stop all services
docker-compose down

# Remove all data (WARNING: destroys data!)
docker-compose down -v

# Remove all images
docker-compose down --rmi all

# Start fresh
docker-compose up --build -d
```

---

## Backup & Recovery

### Database Backup
```bash
# Create backup
docker-compose exec db pg_dump -U pharmasight_user pharmasight_db > backup.sql

# Restore backup
cat backup.sql | docker-compose exec -T db psql -U pharmasight_user pharmasight_db
```

### Redis Backup
```bash
# Trigger save
docker-compose exec redis redis-cli BGSAVE

# Copy RDB file
docker cp $(docker-compose ps -q redis):/data/dump.rdb ./redis-backup.rdb
```

---

## Maintenance

### Update Services
```bash
# Pull latest images
docker-compose pull

# Rebuild and restart
docker-compose up -d --build
```

### Clean Up
```bash
# Remove unused containers
docker container prune

# Remove unused images
docker image prune

# Remove unused volumes (careful!)
docker volume prune
```

---

## Additional Resources

- **API Documentation:** See `API_DOCUMENTATION.md`
- **README:** See `README.md`
- **Docker Compose Reference:** https://docs.docker.com/compose/
- **FastAPI Documentation:** https://fastapi.tiangolo.com/

---

**Last Updated:** November 2024  
**Version:** 1.0.0
