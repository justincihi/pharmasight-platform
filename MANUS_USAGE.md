# Manus Agent Usage Guide

## 🎯 Quick Start

### Branch to Use
```bash
claude/fix-todo-comment-8Pkt3
```

### One-Command Test
```bash
./manus-test.sh
```

This automated script will:
- ✅ Test all 12 services
- ✅ Verify ketamine pipeline
- ✅ Test BioTransformer metabolite prediction
- ✅ Run integration tests
- ✅ Generate comprehensive test report

---

## 📋 Configuration File

The Manus agent can read the configuration from:
```
.manus/pharmasight-test-config.json
```

This JSON file defines:
- **11 test tasks** with priorities
- **Environment variables** needed
- **Service ports** and health checks
- **Expected results** and success criteria
- **Troubleshooting** tips

---

## 🚀 Manual Task Execution

If you prefer to run tests manually, use these commands:

### 1. Test Ketamine Pipeline
```bash
cd ketamine-pipeline
pip3 install -r requirements.txt
python3 main.py
```

**Expected Output:**
- `results/ketamine_aryl_analogs_ip_labels.json` - IP classified analogs
- `data/ketamine_3d.sdf` - 3D conformer for docking

### 2. Start All Services
```bash
docker-compose up -d
```

**Verify:**
```bash
docker-compose ps
# All services should show "Up" and "healthy"
```

### 3. Test BioTransformer
```bash
# Health check
curl http://localhost:8007/health

# Predict metabolites
curl -X POST http://localhost:8007/predict \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",
    "metabolism_type": "human",
    "steps": 2
  }'
```

### 4. Setup Admin Dashboard
```bash
cd admin-dashboard

# Create .env file
cat > .env << EOF
DATABASE_URL=mysql://root:pharmasight_mysql_2024@localhost:3306/pharmasight
PLATFORM_API_KEY=test-api-key-12345
NODE_ENV=development
EOF

# Install dependencies
npm install

# Run database migrations
npx drizzle-kit push

# Start dev server
npm run dev
```

**Access:** http://localhost:3000

### 5. Test All Microservices
```bash
# Test each health endpoint
curl http://localhost:8001/health  # Compound Service
curl http://localhost:8002/health  # Analog Service
curl http://localhost:8003/health  # ML Service
curl http://localhost:8004/health  # Quantum Calculator
curl http://localhost:8005/health  # Auth Service
curl http://localhost:8006/health  # PopHIVE Connector
curl http://localhost:8007/health  # BioTransformer
curl http://localhost:8080/health  # API Gateway
```

---

## 📊 Service Ports Reference

| Service | Port | Health Check URL |
|---------|------|------------------|
| Admin Dashboard | 3000 | http://localhost:3000/api/health |
| Compound Service | 8001 | http://localhost:8001/health |
| Analog Service | 8002 | http://localhost:8002/health |
| ML Service | 8003 | http://localhost:8003/health |
| Quantum Calculator | 8004 | http://localhost:8004/health |
| Auth Service | 8005 | http://localhost:8005/health |
| PopHIVE Connector | 8006 | http://localhost:8006/health |
| BioTransformer | 8007 | http://localhost:8007/health |
| API Gateway | 8080 | http://localhost:8080/health |
| MySQL | 3306 | - |
| PostgreSQL | 5432 | - |
| Redis | 6379 | - |

---

## 🎯 Priority Tasks for Manus Agent

### HIGH Priority

1. **Admin Dashboard Setup & Test** (15 min)
   - Install npm dependencies
   - Run database migrations
   - Start dev server
   - Test chatbot interface
   - **Value:** Validates most complex component

2. **BioTransformer Integration Test** (5 min)
   - Start service
   - Test metabolite prediction
   - Verify mock mode works
   - **Value:** Confirms new feature works

3. **End-to-End Workflow** (15 min)
   - Generate analog with ketamine pipeline
   - Predict metabolites with BioTransformer
   - Verify data flows through system
   - **Value:** Proves complete integration

### MEDIUM Priority

4. **Microservices Health Check** (5 min)
   - Start all Python services
   - Verify health endpoints
   - Check Docker resource usage
   - **Value:** Ensures base infrastructure works

5. **Database Setup** (5 min)
   - Start MySQL, PostgreSQL, Redis
   - Verify connectivity
   - Run migrations
   - **Value:** Foundation for all services

### LOW Priority

6. **Performance Testing** (20 min)
   - Batch analog generation
   - API response time testing
   - Resource usage monitoring
   - **Value:** Nice to have, not critical

---

## 🔧 Environment Variables

### Required (Included in Script)
```bash
MYSQL_ROOT_PASSWORD=pharmasight_mysql_2024
MYSQL_DATABASE=pharmasight
POSTGRES_USER=pharmasight_user
POSTGRES_PASSWORD=pharmasight_pass_2024
POSTGRES_DB=pharmasight_db
PLATFORM_API_KEY=test-api-key-12345
```

### Optional (For Full Functionality)
```bash
# LLM APIs for chatbot
GEMINI_API_KEY=your_key_here
ANTHROPIC_API_KEY=your_key_here
SONAR_API_KEY=your_key_here

# PubMed for autonomous research
PUBMED_API_KEY=your_key_here
```

---

## 📝 Expected Test Results

### Success Criteria
- ✅ All 12 services start successfully
- ✅ All health checks pass
- ✅ Ketamine pipeline generates 3+ analogs
- ✅ BioTransformer predicts metabolites (mock mode OK)
- ✅ Admin dashboard accessible at http://localhost:3000
- ✅ API response times < 2000ms
- ✅ No critical errors in logs

### Test Report Output
After running `./manus-test.sh`, you'll get:
- **MANUS_TEST_REPORT.md** - Comprehensive test report
- **manus-test-results.log** - Detailed execution log

---

## ⚠️ Common Issues & Solutions

### Issue: MySQL Connection Refused
**Solution:** Wait 30 seconds after `docker-compose up` for MySQL to initialize
```bash
docker-compose up -d mysql
sleep 30
docker-compose exec mysql mysqladmin ping
```

### Issue: Admin Dashboard npm Install Fails
**Solution:** Use Node.js 18+
```bash
nvm install 18
nvm use 18
cd admin-dashboard && npm install
```

### Issue: BioTransformer Returns Mock Data
**Solution:** This is expected! BioTransformer JAR needs manual download
- Mock mode is perfect for testing integration
- Service works correctly, just returns sample metabolites
- Check response for `"mock_mode": true`

### Issue: Port Already in Use
**Solution:** Stop conflicting services
```bash
docker-compose down
sudo lsof -ti:3000 | xargs kill -9  # Kill process on port 3000
sudo lsof -ti:8080 | xargs kill -9  # Kill process on port 8080
```

### Issue: RDKit Import Errors
**Solution:** Install RDKit
```bash
pip3 install rdkit
```

---

## 📚 Documentation References

- **Platform Overview:** `INTEGRATION_COMPLETE_SUMMARY.md`
- **Dragonfly Integration:** `DRAGONFLY_INTEGRATION_PLAN.md`
- **Ketamine Pipeline:** `ketamine-pipeline/README.md`
- **BioTransformer:** `services/biotransformer/README.md`
- **Platform API:** `PLATFORM_API.md`

---

## 🎁 What to Focus On

### Most Valuable Tests for Manus
1. **Admin Dashboard** - Complex TypeScript/React UI
2. **Database Setup** - MySQL schema and migrations
3. **Integration Testing** - Full workflow validation
4. **Performance Metrics** - Response times and resource usage

### Skip These (Better Done Locally)
- ❌ Heavy ML model training
- ❌ Large file downloads (BioTransformer JAR, Dragonfly models)
- ❌ GPU-intensive tasks (protein folding, MD simulations)

---

## 🚀 Quick Test Commands

### Fast Health Check
```bash
# Check all services in one command
for port in 3000 8001 8002 8003 8004 8005 8006 8007 8080; do
  echo -n "Port $port: "
  curl -s -f http://localhost:$port/health > /dev/null && echo "✅" || echo "❌"
done
```

### View All Logs
```bash
docker-compose logs -f
```

### Restart Everything
```bash
docker-compose down
docker-compose up -d
```

### Clean Slate
```bash
docker-compose down -v  # WARNING: Deletes all data!
docker system prune -f
```

---

## 📞 Support

If tests fail, check:
1. **Test Report:** `MANUS_TEST_REPORT.md`
2. **Full Log:** `manus-test-results.log`
3. **Service Logs:** `docker-compose logs [service-name]`
4. **Troubleshooting:** `.manus/pharmasight-test-config.json` → `troubleshooting` section

---

## ✅ Success Checklist

After Manus testing, you should have:
- [ ] All services running (`docker-compose ps` shows all "Up")
- [ ] Admin dashboard accessible at http://localhost:3000
- [ ] Ketamine pipeline generates analogs
- [ ] BioTransformer predicts metabolites
- [ ] Test report shows 0 critical errors
- [ ] Database contains 150+ analogs
- [ ] All health endpoints return 200 OK

**When complete:** Review `MANUS_TEST_REPORT.md` and proceed with production deployment!

---

**Version:** 1.0
**Last Updated:** January 24, 2026
**Branch:** claude/fix-todo-comment-8Pkt3
