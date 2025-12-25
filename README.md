# PharmaSight™ Platform - Combined Edition

![Status](https://img.shields.io/badge/status-integrated-success)
![Version](https://img.shields.io/badge/version-4.0.0--combined-blue)
![License](https://img.shields.io/badge/license-Proprietary-red)

**Comprehensive AI-Powered Drug Discovery Platform with Integrated Admin Dashboard**

---

## 🎯 Overview

PharmaSight™ Combined Edition integrates the complete drug discovery platform with a modern administrative dashboard, providing researchers and administrators with a unified, powerful system for pharmaceutical research and platform management.

### What's New in Combined Edition

This branch merges the best features from multiple development streams:

✅ **Complete Drug Discovery Platform** (from Replit-December)
- 40+ specialized drug discovery modules
- BioTransformer metabolism prediction integration
- Advanced molecular editor (Phase 2)
- Comprehensive data export capabilities (Phase 3)
- Full RDKit integration for cheminformatics

✅ **Modern Admin Dashboard** (from manus-December)
- TypeScript/React-based management interface
- User management and authentication
- Database administration tools
- LLM API configuration interface
- Analytics and reporting

✅ **Production-Ready Infrastructure**
- Docker containerization
- Multi-subdomain architecture
- Shared PostgreSQL database
- Redis caching layer
- Nginx reverse proxy

---

## 🏗️ Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    PharmaSight Platform                      │
├─────────────────────────────────────────────────────────────┤
│                                                               │
│  ┌──────────────────────┐      ┌──────────────────────┐    │
│  │   Flask Application  │      │   Admin Dashboard    │    │
│  │  (Drug Discovery)    │      │   (Management UI)    │    │
│  │                      │      │                      │    │
│  │  • 40+ Modules       │      │  • User Management   │    │
│  │  • BioTransformer    │◄────►│  • Analytics         │    │
│  │  • Molecular Editor  │      │  • LLM Config        │    │
│  │  • Data Export       │      │  • DB Admin          │    │
│  │                      │      │                      │    │
│  │  Port: 5000          │      │  Port: 3000          │    │
│  │  app.pharmasight.com │      │  admin.pharmasight.com│   │
│  └──────────┬───────────┘      └──────────┬───────────┘    │
│             │                               │                │
│             └───────────┬───────────────────┘                │
│                         │                                    │
│              ┌──────────▼──────────┐                        │
│              │  PostgreSQL Database │                        │
│              │   (Shared Storage)   │                        │
│              └──────────────────────┘                        │
│                                                               │
└─────────────────────────────────────────────────────────────┘
```

---

## 🚀 Quick Start

### Option 1: Docker Compose (Recommended)

```bash
# Clone repository
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
git checkout Combined

# Configure environment
cp .env.example .env
# Edit .env with your settings

# Start all services
docker-compose up --build

# Access applications
# Main app: http://localhost:5000
# Admin dashboard: http://localhost:3000
```

### Option 2: Manual Setup

**Flask Application:**
```bash
# Install dependencies
pip install -r requirements.txt
pip install -r requirements-rdkit.txt

# Run application
python app.py
```

**Admin Dashboard:**
```bash
# Navigate to admin dashboard
cd admin-dashboard

# Install dependencies
pnpm install

# Start development server
pnpm dev
```

---

## 📦 Features

### Drug Discovery Platform (Flask Application)

#### Core Capabilities
- **Compound Analysis** - Comprehensive molecular property calculation
- **Analog Generation** - AI-powered chemical analog discovery
- **Drug-Drug Interaction** - Interaction prediction and risk assessment
- **ADMET Prediction** - Absorption, distribution, metabolism, excretion, toxicity
- **Molecular Docking** - Protein-ligand binding prediction
- **Patent Intelligence** - Patent landscape analysis

#### Advanced Features
- **BioTransformer Integration** - Metabolism prediction using open-source tool
- **Molecular Editor** - Interactive structure editing with RDKit
- **Data Export** - Export to CSV, Excel, PDF formats
- **Quantum Calculations** - DFT calculations for molecular properties
- **Pharmacophore Modeling** - 3D pharmacophore generation
- **Retrosynthesis** - Synthetic route planning
- **Off-Target Prediction** - Unintended target identification
- **PK/PD Simulation** - Pharmacokinetic/pharmacodynamic modeling

#### Research Tools
- **Autonomous Research Engine** - 24/7 automated literature mining
- **Research Findings Database** - Comprehensive research data storage
- **Receptor Profiling** - Multi-receptor binding analysis
- **SAR Explorer** - Structure-activity relationship analysis
- **Virtual Screening** - High-throughput compound screening

### Admin Dashboard (TypeScript/React Application)

#### Management Features
- **User Management** - Create, edit, delete user accounts
- **Role-Based Access Control** - Granular permission management
- **Database Administration** - View and manage compound data
- **Analytics Dashboard** - Platform usage statistics
- **LLM API Configuration** - Configure OpenAI, Perplexity, Gemini keys
- **System Monitoring** - Health checks and performance metrics
- **Audit Logging** - Track all administrative actions

#### Integration Features
- **Flask API Integration** - Seamless communication with main app
- **Shared Database** - Unified data access across both applications
- **JWT Authentication** - Secure cross-application authentication
- **Real-time Updates** - Live data synchronization

---

## 🛠️ Technology Stack

### Flask Application
- **Backend**: Python 3.11+, Flask 3.1.1
- **Chemistry**: RDKit 2025.9.1
- **ML/AI**: scikit-learn, scipy, numpy
- **Docking**: AutoDock Vina
- **Database**: SQLAlchemy ORM
- **Server**: Gunicorn

### Admin Dashboard
- **Frontend**: React 18, Vite, TailwindCSS
- **Backend**: Node.js 18+, Express
- **Database**: Drizzle ORM
- **Language**: TypeScript
- **Package Manager**: pnpm

### Infrastructure
- **Containerization**: Docker, Docker Compose
- **Database**: PostgreSQL 15
- **Cache**: Redis 7
- **Reverse Proxy**: Nginx
- **SSL**: Let's Encrypt (certbot)

---

## 📚 Documentation

- **[Deployment Guide](COMBINED_DEPLOYMENT_GUIDE.md)** - Complete deployment instructions
- **[Admin Dashboard README](admin-dashboard/README.md)** - Dashboard-specific documentation
- **[API Documentation](API_DOCUMENTATION.md)** - API endpoints and usage
- **[Branch Analysis](../branch_analysis.md)** - Detailed comparison of merged branches

---

## 🌐 Deployment Options

### Local Development
- Docker Compose for full stack
- Individual services for development
- Hot reload enabled

### Production Deployment

#### VPS/Dedicated Server
- Docker Compose deployment
- Nginx reverse proxy
- SSL with Let's Encrypt
- Subdomain routing

#### Render (PaaS)
- Separate web services for Flask and Admin
- Shared PostgreSQL database
- Custom domain support
- Automatic SSL

#### AWS/GCP/Azure
- Container orchestration (ECS, GKE, AKS)
- Managed databases (RDS, Cloud SQL)
- Load balancing
- Auto-scaling

---

## 🔒 Security

- JWT-based authentication
- Role-based access control (RBAC)
- Environment variable secrets
- HTTPS/SSL encryption
- CORS configuration
- Rate limiting
- SQL injection protection
- XSS prevention

---

## 📊 Database Schema

Both applications share a PostgreSQL database with the following main tables:

- **users** - User accounts and authentication
- **compounds** - Chemical compound data
- **research_findings** - Research results and analytics
- **analog_discoveries** - Generated chemical analogs
- **ddi_interactions** - Drug-drug interaction data
- **patents** - Patent information
- **audit_logs** - System activity tracking

---

## 🧪 Testing

```bash
# Test Flask application
python test_core_features.py

# Test admin dashboard
cd admin-dashboard
pnpm test

# Integration tests
docker-compose up -d
python test_integration.py
```

---

## 📈 Performance

- **Response Time**: <2 seconds for most operations
- **Concurrent Users**: Supports 100+ simultaneous users
- **Database**: Optimized queries with indexing
- **Caching**: Redis for expensive calculations
- **Scalability**: Horizontal scaling with load balancer

---

## 🔄 Updates and Maintenance

### Updating the Platform

```bash
# Pull latest changes
git pull origin Combined

# Rebuild containers
docker-compose up -d --build

# Run migrations if needed
docker-compose exec flask-app python migrate.py
docker-compose exec admin-dashboard pnpm drizzle-kit push:pg
```

### Database Backups

```bash
# Backup
docker exec pharmasight-db pg_dump -U pharmasight pharmasight > backup.sql

# Restore
docker exec -i pharmasight-db psql -U pharmasight pharmasight < backup.sql
```

---

## 🤝 Contributing

This is a proprietary platform. For authorized contributors:

1. Create feature branch from Combined
2. Make changes and test thoroughly
3. Submit pull request with detailed description
4. Await code review and approval

---

## 📝 License

Copyright © 2025 PharmaSight™ Platform. All rights reserved.

This is proprietary software. Unauthorized copying, modification, or distribution is prohibited.

---

## 🎓 Key Integrations

### Open Source Tools
- **BioTransformer 3.0** - Metabolism prediction (GPL v2.1)
- **RDKit** - Cheminformatics toolkit (BSD)
- **AutoDock Vina** - Molecular docking (Apache 2.0)

### External APIs (Optional)
- **PubChem** - Chemical database
- **ChEMBL** - Bioactivity database
- **OpenAI** - LLM capabilities
- **Perplexity** - Research augmentation
- **Google Gemini** - AI analysis

---

## 🎯 Roadmap

### Immediate (Q1 2025)
- [ ] Enhanced LLM integration in admin dashboard
- [ ] Real-time collaboration features
- [ ] Mobile-responsive admin interface
- [ ] Advanced analytics and reporting

### Short-term (Q2 2025)
- [ ] Clinical trial integration
- [ ] Enhanced molecular visualization
- [ ] API rate limiting and quotas
- [ ] Multi-language support

### Long-term (Q3-Q4 2025)
- [ ] Machine learning model training interface
- [ ] Laboratory information system integration
- [ ] Advanced workflow automation
- [ ] Mobile applications (iOS/Android)

---

## 💡 Use Cases

1. **Drug Discovery Research**
   - Analyze compounds and generate analogs
   - Predict ADMET properties
   - Identify patent-free opportunities

2. **Platform Administration**
   - Manage user accounts and permissions
   - Monitor system performance
   - Configure AI integrations

3. **Data Management**
   - Export research data in multiple formats
   - Backup and restore databases
   - Audit trail review

4. **Collaboration**
   - Share research findings
   - Team-based compound analysis
   - Project management

---

## 📞 Support

For questions, issues, or feature requests:
- Review documentation first
- Check existing GitHub issues
- Contact development team
- Submit detailed bug reports

---

## 🏆 Achievements

✅ 40+ specialized drug discovery modules integrated  
✅ Modern admin dashboard with full CRUD operations  
✅ BioTransformer metabolism prediction integrated  
✅ Multi-subdomain architecture implemented  
✅ Docker containerization complete  
✅ Comprehensive documentation provided  
✅ Production-ready deployment configurations  

---

**PharmaSight™ Combined Edition** - Where cutting-edge drug discovery meets modern platform management.

*Advancing pharmaceutical research through artificial intelligence and innovative software architecture.*
