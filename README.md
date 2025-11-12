# PharmaSight™ - Microservices Architecture with Autonomous Research

**Version:** 2.0  
**Last Updated:** 2025-11-11  
**Status:** Production Ready ✅

---

## 🎯 Overview

PharmaSight™ is an enterprise drug discovery platform featuring a modern microservices architecture with integrated autonomous research capabilities. The platform combines molecular analysis, analog generation, machine learning predictions, and automated literature scanning to accelerate pharmaceutical research and IP discovery.

---

## 🏗️ Architecture

### Microservices

| Service | Port | Description |
|---------|------|-------------|
| **API Gateway** | 8080 | Single entry point for all client requests |
| **Compound Analysis** | 8001 | Molecular analysis, RDKit calculations, property prediction |
| **Analog Generation** | 8002 | Chemical analog generation, similarity scoring, patent analysis |
| **ML Models** | 8003 | ADMET predictions, toxicity models, ML-based property prediction |
| **Quantum Calculator** | 8004 | DFT calculations using PySCF for advanced properties |
| **Auth Service** | 8005 | JWT authentication, authorization, RBAC |
| **Research Engine** | 8006 | **NEW**: Autonomous research, article database, analog discovery |

### Infrastructure

- **PostgreSQL 13**: Primary relational database
- **Redis**: In-memory caching and session management
- **Docker Compose**: Container orchestration
- **Health Checks**: Automated service monitoring

---

## ✨ Key Features

### Core Capabilities
- ✅ **500+ Compounds** in database with complete profiles
- ✅ **70+ Receptor Subtypes** with pharmacology data
- ✅ **RDKit Integration** for molecular calculations
- ✅ **ADMET Predictions** using machine learning
- ✅ **Quantum Chemistry** calculations (DFT)
- ✅ **PKPD/PBPK Simulations** for population modeling
- ✅ **Drug-Drug Interactions** analysis

### Autonomous Research System
- ✅ **Daily Literature Scanning** via PubMed API
- ✅ **26 Research Articles** in database (1985-2024)
- ✅ **119 Novel Analogs** discovered and documented
- ✅ **IP Protection** with public timestamped registry
- ✅ **Rate Limiting** (50 API calls/day for cost control)
- ✅ **Automatic Analog Generation** from research goals

### Enterprise Features
- ✅ **Microservices Architecture** for scalability
- ✅ **JWT Authentication** with role-based access
- ✅ **Health Monitoring** for all services
- ✅ **Audit Logging** for compliance
- ✅ **Data Export** (PDF, CSV, Excel)

---

## 🚀 Quick Start

### Prerequisites

- Docker 20.10+
- Docker Compose 1.29+
- Python 3.11+ (for testing)
- 4GB RAM minimum (8GB recommended)

### Installation

1. **Clone Repository**
```bash
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
git checkout microservices-with-research
```

2. **Configure Environment (Optional)**
```bash
cp .env.example .env
# Edit .env with your settings
```

3. **Start All Services**
```bash
docker-compose up --build -d
```

Or use the Makefile:
```bash
make up
```

4. **Verify Services**
```bash
curl http://localhost:8080/health
```

Or:
```bash
make health
```

5. **Run Tests**
```bash
python3 test_microservices.py
```

Or:
```bash
make test
```

### Accessing the Platform

- **API Gateway:** http://localhost:8080
- **Interactive API Docs:** http://localhost:8080/docs
- **Research Engine:** http://localhost:8006
- **Research Engine Docs:** http://localhost:8006/docs

---

## 📊 Research Engine

### Features

The Research Engine microservice provides:

1. **Autonomous Literature Scanning**
   - Automated PubMed searches
   - Article metadata extraction
   - DOI, PMID, authors, journals
   - Session logging

2. **Article Database**
   - 26 research articles
   - Searchable by keyword, year, journal
   - Statistics and analytics
   - CSV export

3. **Analog Discovery**
   - 119 novel analogs
   - Complete SMILES strings
   - Patent status tracking
   - IP-protected registry

4. **RDKit Integration**
   - Automatic analog generation
   - Research goal synchronization
   - Property calculations
   - Master discovery log

### API Examples

**Run Research Cycle:**
```bash
curl -X POST http://localhost:8006/research/run-cycle \
  -H "Content-Type: application/json" \
  -d '{
    "goals": [
      "psilocybin 5-HT2A receptor depression",
      "ketamine NMDA antagonist"
    ]
  }'
```

**Search Articles:**
```bash
curl -X POST http://localhost:8006/articles/search \
  -H "Content-Type: application/json" \
  -d '{"keyword": "psilocybin", "limit": 5}'
```

**Get All Discoveries:**
```bash
curl http://localhost:8006/discoveries/all
```

---

## 🔄 Daily Automation

### Automated Research Cycle

Schedule daily research using cron:

```bash
# Edit crontab
crontab -e

# Add this line (runs daily at 9 AM)
0 9 * * * /path/to/run_daily_research.sh
```

The script `run_daily_research.sh` is provided and configured for:
- 5 research goals
- 50 API calls/day limit
- Automatic logging
- Error handling

---

## 📁 Data Files

### Research Articles
- `RESEARCH_ARTICLES_DATABASE.json` - Complete article database
- `RESEARCH_ARTICLES_DATABASE.csv` - Spreadsheet format
- `RESEARCH_ARTICLES_README.md` - Documentation

### Analog Discoveries
- `MASTER_ANALOG_DISCOVERIES.json` - All discoveries with timestamps
- `PUBLIC_ANALOG_DISCOVERY_REGISTRY.md` - Public IP registry

### Session Logs
- `research_session_*.json` - Research cycle logs
- `rdkit_sync_*.json` - RDKit integration logs

---

## 🔧 Common Operations

### Using Docker Compose

```bash
# Start services
docker-compose up -d

# Stop services
docker-compose down

# View logs
docker-compose logs -f

# View specific service logs
docker-compose logs -f research-engine

# Rebuild a service
docker-compose up --build research-engine

# Restart a service
docker-compose restart research-engine
```

### Using Makefile

```bash
make up          # Start all services
make down        # Stop all services
make restart     # Restart all services
make logs        # View all logs
make health      # Check service health
make test        # Run integration tests
make clean       # Clean up containers and volumes
```

---

## 🧪 Testing

### Health Checks

```bash
# Check all services
curl http://localhost:8080/health

# Check individual services
curl http://localhost:8001/health  # Compound Analysis
curl http://localhost:8002/health  # Analog Generation
curl http://localhost:8003/health  # ML Models
curl http://localhost:8004/health  # Quantum Calculator
curl http://localhost:8005/health  # Auth Service
curl http://localhost:8006/health  # Research Engine
```

### Integration Tests

```bash
# Run full test suite
python3 test_microservices.py

# Run specific module tests
python3 test_modules_direct.py

# Test external APIs
python3 test_external_apis.py
```

---

## 📚 Documentation

### Complete Documentation Set

- **README.md** (this file) - Main documentation
- **MICROSERVICES_RESEARCH_INTEGRATION.md** - Integration details
- **MICROSERVICES_DEPLOYMENT.md** - Deployment guide
- **AUTONOMOUS_RESEARCH_SYSTEM_COMPLETE.md** - Research system docs
- **API_DOCUMENTATION.md** - API reference
- **SECURITY_SUMMARY.md** - Security considerations
- **QUICK_START.md** - Getting started guide

---

## 🔒 Security

### Authentication

- JWT-based authentication via Auth Service
- Role-based access control (RBAC)
- Secure password hashing
- Token expiration and refresh

### Data Protection

- PostgreSQL with connection pooling
- Redis for secure session management
- Environment variable configuration
- Secrets management

### Network Security

- Internal Docker network isolation
- API Gateway as single entry point
- CORS configuration
- Rate limiting

---

## 📈 Statistics

### Research Data

- **Research Articles:** 26 (1985-2024)
- **DOI Coverage:** 96%
- **Unique Journals:** 26

### Analog Discoveries

- **Total Analogs:** 119
- **Patent-Free:** 104 (87%)
- **High IP Opportunity:** 110 (92%)
- **Compounds Covered:** 12 parent compounds

### API Usage

- **Daily Limit:** 50 calls
- **Typical Usage:** 9-10 calls/cycle (18-20%)
- **Cost:** Minimal (free PubMed API)

---

## 🌐 Deployment

### Production Deployment

The platform is ready for deployment to:
- **Railway** (recommended)
- **Render**
- **Fly.io**
- **Heroku**
- **AWS ECS**
- **Google Cloud Run**
- **Azure Container Instances**

### Configuration Files

- `docker-compose.yml` - Local development
- `Dockerfile` - Each service has its own
- `.env.example` - Environment template
- `railway.json` - Railway configuration

---

## 🆘 Troubleshooting

### Services Won't Start

```bash
# Check Docker is running
docker ps

# Check logs for errors
docker-compose logs

# Rebuild all services
docker-compose up --build
```

### Database Connection Issues

```bash
# Check PostgreSQL is healthy
docker-compose ps db

# View database logs
docker-compose logs db

# Restart database
docker-compose restart db
```

### Port Conflicts

```bash
# Check what's using ports
netstat -tuln | grep -E "8080|8001|8002|8003|8004|8005|8006"

# Stop conflicting services or change ports in docker-compose.yml
```

---

## 🤝 Contributing

### Development Workflow

1. Create feature branch
2. Make changes
3. Run tests
4. Update documentation
5. Submit pull request

### Code Style

- Python: PEP 8
- Type hints encouraged
- Docstrings for all functions
- Comments for complex logic

---

## 📞 Support

For questions or issues:
- **GitHub Issues:** https://github.com/justincihi/pharmasight-platform/issues
- **Documentation:** See docs folder
- **Email:** Submit feedback at https://help.manus.im

---

## 📄 License

Proprietary - All Rights Reserved

---

## 🎉 Acknowledgments

- **RDKit** for molecular calculations
- **FastAPI** for microservices framework
- **PostgreSQL** for reliable data storage
- **Docker** for containerization
- **PubMed** for research article access

---

**Version:** 2.0  
**Branch:** microservices-with-research  
**Status:** Production Ready ✅  
**Last Updated:** 2025-11-11

