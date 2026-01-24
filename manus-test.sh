#!/bin/bash

################################################################################
# PharmaSight Platform - Manus Agent Integration Test Script
# Branch: claude/fix-todo-comment-8Pkt3
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test results tracking
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# Log file
LOG_FILE="manus-test-results.log"
echo "PharmaSight Platform Test Run - $(date)" > "$LOG_FILE"

################################################################################
# Helper Functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOG_FILE"
    ((TESTS_PASSED++))
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
    ((TESTS_FAILED++))
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

log_skip() {
    echo -e "${YELLOW}[SKIP]${NC} $1" | tee -a "$LOG_FILE"
    ((TESTS_SKIPPED++))
}

wait_for_service() {
    local service_name=$1
    local url=$2
    local max_attempts=30
    local attempt=0

    log_info "Waiting for $service_name to be ready..."

    while [ $attempt -lt $max_attempts ]; do
        if curl -s -f "$url" > /dev/null 2>&1; then
            log_success "$service_name is ready"
            return 0
        fi
        sleep 2
        ((attempt++))
    done

    log_error "$service_name failed to start after $max_attempts attempts"
    return 1
}

################################################################################
# Test 1: Environment Setup
################################################################################

test_environment_setup() {
    log_info "========================================="
    log_info "Test 1: Environment Setup"
    log_info "========================================="

    # Check Docker
    if command -v docker &> /dev/null; then
        log_success "Docker is installed: $(docker --version)"
    else
        log_error "Docker is not installed"
        return 1
    fi

    # Check Docker Compose
    if command -v docker-compose &> /dev/null; then
        log_success "Docker Compose is installed: $(docker-compose --version)"
    else
        log_error "Docker Compose is not installed"
        return 1
    fi

    # Check Node.js
    if command -v node &> /dev/null; then
        log_success "Node.js is installed: $(node --version)"
    else
        log_warning "Node.js is not installed (required for admin dashboard)"
    fi

    # Check Python
    if command -v python3 &> /dev/null; then
        log_success "Python3 is installed: $(python3 --version)"
    else
        log_error "Python3 is not installed"
        return 1
    fi

    # Create .env file if it doesn't exist
    if [ ! -f .env ]; then
        log_info "Creating .env file with default values..."
        cat > .env << EOF
# Database Configuration
MYSQL_ROOT_PASSWORD=pharmasight_mysql_2024
MYSQL_DATABASE=pharmasight
POSTGRES_USER=pharmasight_user
POSTGRES_PASSWORD=pharmasight_pass_2024
POSTGRES_DB=pharmasight_db

# Platform API
PLATFORM_API_KEY=test-api-key-12345

# LLM API Keys (optional - add your own)
# GEMINI_API_KEY=your_key_here
# ANTHROPIC_API_KEY=your_key_here
# SONAR_API_KEY=your_key_here
# PUBMED_API_KEY=your_key_here

# Service Configuration
SECRET_KEY=change-this-secret-key-in-production-2024
EOF
        log_success "Created .env file with default configuration"
    else
        log_success ".env file already exists"
    fi

    # Create admin-dashboard .env
    if [ ! -f admin-dashboard/.env ]; then
        log_info "Creating admin-dashboard/.env file..."
        cat > admin-dashboard/.env << EOF
DATABASE_URL=mysql://root:pharmasight_mysql_2024@localhost:3306/pharmasight
PLATFORM_API_KEY=test-api-key-12345
NODE_ENV=development
EOF
        log_success "Created admin-dashboard/.env file"
    else
        log_success "admin-dashboard/.env file already exists"
    fi
}

################################################################################
# Test 2: Ketamine Pipeline
################################################################################

test_ketamine_pipeline() {
    log_info "========================================="
    log_info "Test 2: Ketamine Analog Pipeline"
    log_info "========================================="

    cd ketamine-pipeline || return 1

    # Install dependencies
    log_info "Installing Python dependencies..."
    if pip3 install -r requirements.txt >> "$LOG_FILE" 2>&1; then
        log_success "Dependencies installed"
    else
        log_error "Failed to install dependencies"
        cd ..
        return 1
    fi

    # Run pipeline
    log_info "Running ketamine pipeline..."
    if python3 main.py >> "$LOG_FILE" 2>&1; then
        log_success "Ketamine pipeline executed successfully"
    else
        log_error "Ketamine pipeline failed"
        cd ..
        return 1
    fi

    # Verify outputs
    if [ -f "results/ketamine_aryl_analogs_ip_labels.json" ]; then
        log_success "IP classification JSON generated"

        # Count analogs
        analog_count=$(python3 -c "import json; data=json.load(open('results/ketamine_aryl_analogs_ip_labels.json')); print(len(data.get('analogs', [])))")
        log_info "Generated $analog_count analogs"

        if [ "$analog_count" -ge 3 ]; then
            log_success "Sufficient analogs generated (>= 3)"
        else
            log_warning "Only $analog_count analogs generated"
        fi
    else
        log_error "IP classification JSON not found"
    fi

    if [ -f "data/ketamine_3d.sdf" ]; then
        log_success "3D SDF file generated"
    else
        log_error "3D SDF file not found"
    fi

    cd ..
}

################################################################################
# Test 3: Start Databases
################################################################################

test_databases() {
    log_info "========================================="
    log_info "Test 3: Database Services"
    log_info "========================================="

    log_info "Starting database services..."
    if docker-compose up -d mysql postgres redis >> "$LOG_FILE" 2>&1; then
        log_success "Database services started"
    else
        log_error "Failed to start database services"
        return 1
    fi

    # Wait for databases
    sleep 15

    # Test MySQL
    if docker-compose exec -T mysql mysqladmin ping -h localhost -uroot -ppharmasight_mysql_2024 >> "$LOG_FILE" 2>&1; then
        log_success "MySQL is healthy"
    else
        log_error "MySQL health check failed"
    fi

    # Test PostgreSQL
    if docker-compose exec -T postgres pg_isready -U pharmasight_user >> "$LOG_FILE" 2>&1; then
        log_success "PostgreSQL is healthy"
    else
        log_error "PostgreSQL health check failed"
    fi

    # Test Redis
    if docker-compose exec -T redis redis-cli ping >> "$LOG_FILE" 2>&1; then
        log_success "Redis is healthy"
    else
        log_error "Redis health check failed"
    fi
}

################################################################################
# Test 4: BioTransformer Service
################################################################################

test_biotransformer() {
    log_info "========================================="
    log_info "Test 4: BioTransformer Service"
    log_info "========================================="

    log_info "Starting BioTransformer service..."
    if docker-compose up -d biotransformer >> "$LOG_FILE" 2>&1; then
        log_success "BioTransformer service started"
    else
        log_error "Failed to start BioTransformer service"
        return 1
    fi

    # Wait for service
    if wait_for_service "BioTransformer" "http://localhost:8007/health"; then

        # Test health endpoint
        health_response=$(curl -s http://localhost:8007/health)
        log_info "Health response: $health_response"

        # Test metabolism types
        log_info "Testing metabolism types endpoint..."
        if curl -s http://localhost:8007/metabolism-types > /dev/null; then
            log_success "Metabolism types endpoint working"
        else
            log_error "Metabolism types endpoint failed"
        fi

        # Test prediction endpoint
        log_info "Testing metabolite prediction..."
        prediction_response=$(curl -s -X POST http://localhost:8007/predict \
            -H "Content-Type: application/json" \
            -d '{"smiles":"CCN(C1CCCCC1=O)c2cccc(F)c2Cl","metabolism_type":"human","steps":1}')

        if echo "$prediction_response" | grep -q "success"; then
            log_success "Metabolite prediction endpoint working"

            # Check if mock mode
            if echo "$prediction_response" | grep -q "mock_mode"; then
                log_warning "BioTransformer running in MOCK mode (JAR not available)"
            else
                log_success "BioTransformer running with real JAR file"
            fi
        else
            log_error "Metabolite prediction failed"
        fi
    else
        log_error "BioTransformer service did not become healthy"
    fi
}

################################################################################
# Test 5: Python Microservices
################################################################################

test_microservices() {
    log_info "========================================="
    log_info "Test 5: Python Microservices"
    log_info "========================================="

    log_info "Starting all Python microservices..."
    if docker-compose up -d compound-service analog-service ml-service \
                         quantum-calculator auth-service api-gateway \
                         pophive-connector >> "$LOG_FILE" 2>&1; then
        log_success "Microservices started"
    else
        log_error "Failed to start microservices"
        return 1
    fi

    sleep 20

    # Test each service
    declare -A services=(
        ["compound-service"]="8001"
        ["analog-service"]="8002"
        ["ml-service"]="8003"
        ["quantum-calculator"]="8004"
        ["auth-service"]="8005"
        ["pophive-connector"]="8006"
        ["api-gateway"]="8080"
    )

    for service in "${!services[@]}"; do
        port=${services[$service]}

        if curl -s -f "http://localhost:$port/health" > /dev/null 2>&1; then
            log_success "$service is healthy (port $port)"
        else
            log_error "$service health check failed (port $port)"
        fi
    done
}

################################################################################
# Test 6: Admin Dashboard
################################################################################

test_admin_dashboard() {
    log_info "========================================="
    log_info "Test 6: Admin Dashboard"
    log_info "========================================="

    cd admin-dashboard || return 1

    # Check if Node.js is available
    if ! command -v node &> /dev/null; then
        log_skip "Node.js not available - skipping admin dashboard tests"
        cd ..
        return 0
    fi

    # Install dependencies
    log_info "Installing npm dependencies (this may take a while)..."
    if npm install >> "$LOG_FILE" 2>&1; then
        log_success "npm dependencies installed"
    else
        log_error "npm install failed"
        cd ..
        return 1
    fi

    # Run database migrations
    log_info "Running database migrations..."
    if npx drizzle-kit push --config=drizzle.config.ts >> "$LOG_FILE" 2>&1; then
        log_success "Database schema created"
    else
        log_warning "Database migrations may have failed (check if already applied)"
    fi

    # Start dev server in background
    log_info "Starting admin dashboard..."
    npm run dev >> "$LOG_FILE" 2>&1 &
    DASHBOARD_PID=$!

    # Wait for dashboard
    if wait_for_service "Admin Dashboard" "http://localhost:3000/api/health"; then
        log_success "Admin dashboard is running"

        # Test API health
        dashboard_health=$(curl -s http://localhost:3000/api/health 2>&1)
        log_info "Dashboard health: $dashboard_health"

        # Kill dashboard
        kill $DASHBOARD_PID 2>/dev/null || true
    else
        log_error "Admin dashboard failed to start"
        kill $DASHBOARD_PID 2>/dev/null || true
    fi

    cd ..
}

################################################################################
# Test 7: Integration Test
################################################################################

test_integration() {
    log_info "========================================="
    log_info "Test 7: End-to-End Integration"
    log_info "========================================="

    # Test workflow: Ketamine pipeline -> BioTransformer -> Database

    # 1. Get analog from ketamine pipeline
    if [ -f "ketamine-pipeline/results/ketamine_aryl_analogs_ip_labels.json" ]; then
        analog_smiles=$(python3 -c "import json; data=json.load(open('ketamine-pipeline/results/ketamine_aryl_analogs_ip_labels.json')); print(data['analogs'][0]['smiles'])")
        log_info "Testing with analog: $analog_smiles"

        # 2. Predict metabolites
        metabolite_result=$(curl -s -X POST http://localhost:8007/predict \
            -H "Content-Type: application/json" \
            -d "{\"smiles\":\"$analog_smiles\",\"metabolism_type\":\"human\",\"steps\":1}")

        if echo "$metabolite_result" | grep -q "success"; then
            num_metabolites=$(echo "$metabolite_result" | python3 -c "import sys, json; print(json.load(sys.stdin).get('num_metabolites', 0))")
            log_success "Predicted $num_metabolites metabolites for analog"
        else
            log_error "Metabolite prediction failed in integration test"
        fi
    else
        log_skip "No ketamine pipeline results available for integration test"
    fi

    # Test compound service API
    log_info "Testing compound analysis service..."
    if curl -s -f http://localhost:8001/health > /dev/null 2>&1; then
        log_success "Compound analysis service accessible"
    else
        log_error "Compound analysis service not accessible"
    fi

    # Test API gateway
    log_info "Testing API gateway..."
    if curl -s -f http://localhost:8080/health > /dev/null 2>&1; then
        log_success "API gateway accessible"
    else
        log_error "API gateway not accessible"
    fi
}

################################################################################
# Test 8: Performance Check
################################################################################

test_performance() {
    log_info "========================================="
    log_info "Test 8: Performance Metrics"
    log_info "========================================="

    # Check Docker resource usage
    log_info "Docker resource usage:"
    docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}" | tee -a "$LOG_FILE"

    # Test API response times
    log_info "Testing API response times..."

    services_to_test=(
        "http://localhost:8001/health:Compound Service"
        "http://localhost:8007/health:BioTransformer"
        "http://localhost:8080/health:API Gateway"
    )

    for service_test in "${services_to_test[@]}"; do
        IFS=':' read -r url name <<< "$service_test"

        response_time=$(curl -o /dev/null -s -w '%{time_total}\n' "$url" 2>/dev/null)
        response_ms=$(echo "$response_time * 1000" | bc)

        if (( $(echo "$response_time < 2.0" | bc -l) )); then
            log_success "$name response time: ${response_ms}ms"
        else
            log_warning "$name response time: ${response_ms}ms (>2000ms)"
        fi
    done
}

################################################################################
# Generate Test Report
################################################################################

generate_report() {
    log_info "========================================="
    log_info "Generating Test Report"
    log_info "========================================="

    cat > MANUS_TEST_REPORT.md << EOF
# PharmaSight Platform - Manus Test Report

**Date:** $(date)
**Branch:** claude/fix-todo-comment-8Pkt3

## Test Summary

- **Tests Passed:** $TESTS_PASSED ✅
- **Tests Failed:** $TESTS_FAILED ❌
- **Tests Skipped:** $TESTS_SKIPPED ⏭️
- **Total Tests:** $((TESTS_PASSED + TESTS_FAILED + TESTS_SKIPPED))

## Service Status

| Service | Port | Status |
|---------|------|--------|
| MySQL | 3306 | $(docker-compose ps mysql | grep -q "Up" && echo "✅ Running" || echo "❌ Stopped") |
| PostgreSQL | 5432 | $(docker-compose ps postgres | grep -q "Up" && echo "✅ Running" || echo "❌ Stopped") |
| Redis | 6379 | $(docker-compose ps redis | grep -q "Up" && echo "✅ Running" || echo "❌ Stopped") |
| Admin Dashboard | 3000 | $(curl -s -f http://localhost:3000/api/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "⚠️ Not running") |
| Compound Service | 8001 | $(curl -s -f http://localhost:8001/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| Analog Service | 8002 | $(curl -s -f http://localhost:8002/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| ML Service | 8003 | $(curl -s -f http://localhost:8003/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| Quantum Calculator | 8004 | $(curl -s -f http://localhost:8004/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| Auth Service | 8005 | $(curl -s -f http://localhost:8005/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| PopHIVE Connector | 8006 | $(curl -s -f http://localhost:8006/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| BioTransformer | 8007 | $(curl -s -f http://localhost:8007/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |
| API Gateway | 8080 | $(curl -s -f http://localhost:8080/health > /dev/null 2>&1 && echo "✅ Healthy" || echo "❌ Unhealthy") |

## Pipeline Tests

### Ketamine Analog Pipeline
- **Status:** $([ -f ketamine-pipeline/results/ketamine_aryl_analogs_ip_labels.json ] && echo "✅ Working" || echo "❌ Failed")
- **Analogs Generated:** $([ -f ketamine-pipeline/results/ketamine_aryl_analogs_ip_labels.json ] && python3 -c "import json; print(len(json.load(open('ketamine-pipeline/results/ketamine_aryl_analogs_ip_labels.json')).get('analogs', [])))" || echo "0")
- **3D SDF Generated:** $([ -f ketamine-pipeline/data/ketamine_3d.sdf ] && echo "✅ Yes" || echo "❌ No")

### BioTransformer Service
- **Status:** $(curl -s http://localhost:8007/health | grep -q "healthy" && echo "✅ Running" || echo "❌ Not running")
- **Mode:** $(curl -s http://localhost:8007/health | grep -q "mock_mode.*true" && echo "Mock (JAR not available)" || echo "Production")

## Issues Found

$(grep "\[ERROR\]" "$LOG_FILE" | sed 's/^/- /')

## Warnings

$(grep "\[WARNING\]" "$LOG_FILE" | sed 's/^/- /')

## Recommendations

EOF

    if [ $TESTS_FAILED -eq 0 ]; then
        cat >> MANUS_TEST_REPORT.md << EOF
✅ **All tests passed!** The platform is ready for use.

### Next Steps:
1. Configure LLM API keys in .env for full chatbot functionality
2. Download BioTransformer JAR for production metabolite prediction
3. Download receptor PDB files for AutoDock Vina setup
4. Implement Dragonfly_gen and Dragonfly optimizer services

EOF
    else
        cat >> MANUS_TEST_REPORT.md << EOF
⚠️ **Some tests failed.** Review the issues above and consult the troubleshooting guide.

### Troubleshooting:
- Check full log: \`cat $LOG_FILE\`
- Verify Docker services: \`docker-compose ps\`
- Check service logs: \`docker-compose logs [service-name]\`
- Review documentation: INTEGRATION_COMPLETE_SUMMARY.md

EOF
    fi

    cat >> MANUS_TEST_REPORT.md << EOF
## Full Test Log

See \`$LOG_FILE\` for complete test output.

---

**Generated by:** Manus Agent Test Script
**Report Version:** 1.0
EOF

    log_success "Test report generated: MANUS_TEST_REPORT.md"
}

################################################################################
# Cleanup Function
################################################################################

cleanup() {
    log_info "========================================="
    log_info "Cleanup"
    log_info "========================================="

    log_info "Stopping all services..."
    docker-compose down >> "$LOG_FILE" 2>&1 || true

    # Kill any remaining processes
    pkill -f "npm run dev" 2>/dev/null || true

    log_success "Cleanup complete"
}

################################################################################
# Main Execution
################################################################################

main() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║   PharmaSight Platform - Integration Test Suite          ║"
    echo "║   Branch: claude/fix-todo-comment-8Pkt3                  ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    # Run tests
    test_environment_setup || true
    test_ketamine_pipeline || true
    test_databases || true
    test_biotransformer || true
    test_microservices || true
    test_admin_dashboard || true
    test_integration || true
    test_performance || true

    # Generate report
    generate_report

    # Summary
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║   Test Summary                                            ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Tests Passed:  $TESTS_PASSED ✅"
    echo "Tests Failed:  $TESTS_FAILED ❌"
    echo "Tests Skipped: $TESTS_SKIPPED ⏭️"
    echo ""
    echo "Report: MANUS_TEST_REPORT.md"
    echo "Log:    $LOG_FILE"
    echo ""

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}✅ All tests passed! Platform is ready.${NC}"
        return 0
    else
        echo -e "${RED}❌ Some tests failed. Review MANUS_TEST_REPORT.md${NC}"
        return 1
    fi
}

# Handle Ctrl+C
trap cleanup EXIT

# Run main
main "$@"
