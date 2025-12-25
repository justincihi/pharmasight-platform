# PharmaSight Combined Platform - Deployment Guide

## Overview

The Combined branch integrates the best features from multiple development branches:

- **Replit-December**: Complete drug discovery platform with 40+ modules, BioTransformer integration, molecular editor, and data export
- **manus-December**: Modern TypeScript/React admin dashboard with database management
- **Infrastructure**: Docker containerization and multi-subdomain architecture

## Architecture

```
PharmaSight Combined Platform
├── Flask Application (Main Drug Discovery Platform)
│   ├── Port: 5000
│   ├── Subdomain: app.pharmasight.com
│   └── Features: All drug discovery modules, RDKit, BioTransformer
│
├── Admin Dashboard (Management Interface)
│   ├── Port: 3000
│   ├── Subdomain: admin.pharmasight.com
│   └── Features: User management, analytics, LLM API configuration
│
├── PostgreSQL Database (Shared)
│   ├── Port: 5432
│   └── Shared by both applications
│
├── Redis Cache
│   ├── Port: 6379
│   └── Used by Flask application
│
└── Nginx Reverse Proxy
    ├── Port: 80/443
    └── Routes traffic to appropriate application
```

## Quick Start

### Prerequisites

- Docker and Docker Compose
- Domain name configured with subdomains (optional for local testing)
- At least 4GB RAM

### Local Development

1. **Clone and checkout Combined branch:**
   ```bash
   git clone https://github.com/justincihi/pharmasight-platform.git
   cd pharmasight-platform
   git checkout Combined
   ```

2. **Configure environment variables:**
   ```bash
   cp .env.example .env
   # Edit .env with your settings
   ```

3. **Start all services:**
   ```bash
   docker-compose up --build
   ```

4. **Access the applications:**
   - Main app: http://localhost:5000
   - Admin dashboard: http://localhost:3000
   - With nginx: http://localhost (routes based on Host header)

### Production Deployment

## Option 1: Docker Compose (Recommended for VPS)

1. **Set up your server:**
   ```bash
   # Install Docker and Docker Compose
   curl -fsSL https://get.docker.com -o get-docker.sh
   sh get-docker.sh
   ```

2. **Configure DNS:**
   - Point `app.pharmasight.com` to your server IP
   - Point `admin.pharmasight.com` to your server IP

3. **Set environment variables:**
   ```bash
   export DB_PASSWORD="secure-password"
   export JWT_SECRET="secure-jwt-secret"
   ```

4. **Deploy:**
   ```bash
   docker-compose up -d --build
   ```

5. **Set up SSL (Let's Encrypt):**
   ```bash
   # Install certbot
   apt-get install certbot python3-certbot-nginx
   
   # Get certificates
   certbot --nginx -d app.pharmasight.com -d admin.pharmasight.com
   ```

## Option 2: Render Deployment

### Deploy Flask Application

1. **Create new Web Service on Render:**
   - Name: `pharmasight-app`
   - Environment: `Docker`
   - Branch: `Combined`
   - Root Directory: `/`

2. **Configure environment variables:**
   ```
   DATABASE_URL=<your-postgres-url>
   REDIS_URL=<your-redis-url>
   FLASK_ENV=production
   PORT=5000
   ```

3. **Add PostgreSQL database:**
   - Create new PostgreSQL instance
   - Copy connection URL to DATABASE_URL

4. **Add Redis cache:**
   - Create new Redis instance
   - Copy connection URL to REDIS_URL

### Deploy Admin Dashboard

1. **Create new Web Service on Render:**
   - Name: `pharmasight-admin`
   - Environment: `Docker`
   - Branch: `Combined`
   - Root Directory: `/admin-dashboard`
   - Dockerfile Path: `admin-dashboard/Dockerfile`

2. **Configure environment variables:**
   ```
   DATABASE_URL=<same-postgres-url-as-flask>
   FLASK_API_URL=https://pharmasight-app.onrender.com
   JWT_SECRET=<secure-secret>
   PORT=3000
   NODE_ENV=production
   ```

### Configure Custom Domain

1. **Add custom domain to Flask app:**
   - Domain: `app.pharmasight.com`
   - Render will provide DNS instructions

2. **Add custom domain to Admin dashboard:**
   - Domain: `admin.pharmasight.com`
   - Follow Render's DNS instructions

## Option 3: Manual Deployment (Without Docker)

### Deploy Flask Application

```bash
# Install Python dependencies
pip install -r requirements.txt
pip install -r requirements-rdkit.txt

# Set environment variables
export DATABASE_URL="postgresql://..."
export REDIS_URL="redis://..."

# Run with gunicorn
gunicorn --bind 0.0.0.0:5000 --workers 4 app:app
```

### Deploy Admin Dashboard

```bash
# Navigate to admin dashboard
cd admin-dashboard

# Install dependencies
pnpm install

# Build application
pnpm build

# Set environment variables
export DATABASE_URL="postgresql://..."
export FLASK_API_URL="http://localhost:5000"
export JWT_SECRET="your-secret"

# Run production server
pnpm start
```

## Environment Variables

### Flask Application (.env)

```env
# Database
DATABASE_URL=postgresql://user:password@host:5432/pharmasight

# Redis
REDIS_URL=redis://localhost:6379/0

# Flask
FLASK_ENV=production
FLASK_APP=app.py
SECRET_KEY=your-flask-secret-key

# Server
PORT=5000

# Optional: LLM API Keys
OPENAI_API_KEY=
PERPLEXITY_API_KEY=
GEMINI_API_KEY=
```

### Admin Dashboard (admin-dashboard/.env)

```env
# Database (same as Flask)
DATABASE_URL=postgresql://user:password@host:5432/pharmasight

# Flask API
FLASK_API_URL=http://localhost:5000

# JWT
JWT_SECRET=your-jwt-secret

# Server
PORT=3000
NODE_ENV=production

# CORS
ALLOWED_ORIGINS=http://localhost:5000,https://app.pharmasight.com
```

## Database Setup

The applications share a PostgreSQL database. Run migrations for both:

### Flask Migrations

```bash
# Create tables (if using SQLAlchemy migrations)
python -c "from src.pharmasight_complete import app; app.app_context().push(); db.create_all()"
```

### Admin Dashboard Migrations

```bash
cd admin-dashboard
pnpm drizzle-kit push:pg
```

## Testing the Deployment

### Health Checks

```bash
# Flask application
curl http://localhost:5000/health

# Admin dashboard
curl http://localhost:3000/api/health
```

### Integration Test

```bash
# Test Flask API
curl -X POST http://localhost:5000/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"compound": "aspirin"}'

# Test admin dashboard API
curl http://localhost:3000/api/users
```

## Monitoring and Logs

### Docker Compose Logs

```bash
# All services
docker-compose logs -f

# Specific service
docker-compose logs -f flask-app
docker-compose logs -f admin-dashboard
```

### Application Logs

- Flask logs: Check stdout or configure logging to file
- Admin dashboard: Check stdout or PM2 logs if using PM2

## Backup and Maintenance

### Database Backup

```bash
# Backup PostgreSQL
docker exec pharmasight-db pg_dump -U pharmasight pharmasight > backup.sql

# Restore
docker exec -i pharmasight-db psql -U pharmasight pharmasight < backup.sql
```

### Update Deployment

```bash
# Pull latest changes
git pull origin Combined

# Rebuild and restart
docker-compose up -d --build
```

## Troubleshooting

### Flask App Won't Start

1. Check logs: `docker-compose logs flask-app`
2. Verify DATABASE_URL is correct
3. Ensure RDKit dependencies are installed
4. Check port 5000 is not in use

### Admin Dashboard Won't Start

1. Check logs: `docker-compose logs admin-dashboard`
2. Verify DATABASE_URL matches Flask app
3. Ensure FLASK_API_URL is reachable
4. Check Node.js version (should be 18+)

### Database Connection Issues

1. Verify PostgreSQL is running: `docker-compose ps`
2. Test connection: `psql $DATABASE_URL`
3. Check firewall rules
4. Verify credentials

### Nginx Routing Issues

1. Check nginx logs: `docker-compose logs nginx`
2. Verify subdomain DNS is configured
3. Test with Host header: `curl -H "Host: app.pharmasight.com" http://localhost`
4. Check nginx.conf syntax: `nginx -t`

### CORS Issues

1. Add domains to ALLOWED_ORIGINS in admin dashboard
2. Configure Flask CORS settings
3. Check browser console for CORS errors

## Security Checklist

- [ ] Change all default passwords
- [ ] Set strong JWT_SECRET
- [ ] Enable HTTPS with valid SSL certificates
- [ ] Configure firewall (only allow 80, 443)
- [ ] Set up database backups
- [ ] Enable rate limiting
- [ ] Configure CORS properly
- [ ] Use environment variables for secrets
- [ ] Enable PostgreSQL SSL
- [ ] Set up monitoring and alerts

## Performance Optimization

### Flask Application

- Increase gunicorn workers: `--workers 8`
- Enable Redis caching
- Use connection pooling for database
- Optimize RDKit calculations (cache results)

### Admin Dashboard

- Enable gzip compression in nginx
- Use CDN for static assets
- Optimize database queries
- Enable browser caching

### Database

- Add indexes on frequently queried columns
- Use connection pooling
- Regular VACUUM and ANALYZE
- Monitor query performance

## Scaling

### Horizontal Scaling

1. **Multiple Flask instances:**
   ```yaml
   flask-app:
     deploy:
       replicas: 3
   ```

2. **Load balancer:**
   - Use nginx upstream with multiple backends
   - Or use cloud load balancer (AWS ALB, etc.)

3. **Database replication:**
   - Set up PostgreSQL read replicas
   - Route read queries to replicas

### Vertical Scaling

- Increase server resources (CPU, RAM)
- Optimize Docker resource limits
- Use faster storage (SSD)

## Support

For issues or questions:
- Check logs first
- Review this guide
- Check GitHub issues
- Contact development team

## Next Steps

After successful deployment:

1. Configure LLM API keys in admin dashboard
2. Create admin user accounts
3. Import initial compound data
4. Set up monitoring and alerts
5. Configure automated backups
6. Test all features thoroughly
7. Set up CI/CD pipeline for updates

---

**PharmaSight Combined Platform** - Integrated drug discovery and administration system.
