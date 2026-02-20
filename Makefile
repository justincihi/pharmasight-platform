.PHONY: help build up down restart logs health test clean

# Default target
help:
	@echo "PharmaSight Microservices - Available Commands"
	@echo "=============================================="
	@echo "  make build      - Build all Docker images"
	@echo "  make up         - Start all services"
	@echo "  make down       - Stop all services"
	@echo "  make restart    - Restart all services"
	@echo "  make logs       - Show logs from all services"
	@echo "  make health     - Check health of all services"
	@echo "  make test       - Run integration tests"
	@echo "  make clean      - Stop and remove all containers and volumes"
	@echo ""
	@echo "Service-specific commands:"
	@echo "  make logs-gateway       - Show API gateway logs"
	@echo "  make logs-compound      - Show compound service logs"
	@echo "  make logs-analog        - Show analog service logs"
	@echo "  make logs-ml            - Show ML service logs"
	@echo "  make logs-quantum       - Show quantum calculator logs"
	@echo "  make logs-auth          - Show auth service logs"
	@echo "  make logs-db            - Show database logs"
	@echo "  make logs-redis         - Show Redis logs"

# Build all Docker images
build:
	docker-compose build

# Start all services
up:
	docker-compose up -d
	@echo "Waiting for services to be ready..."
	@sleep 10
	@make health

# Stop all services
down:
	docker-compose down

# Restart all services
restart: down up

# Show logs from all services
logs:
	docker-compose logs -f

# Service-specific logs
logs-gateway:
	docker-compose logs -f api-gateway

logs-compound:
	docker-compose logs -f compound-service

logs-analog:
	docker-compose logs -f analog-service

logs-ml:
	docker-compose logs -f ml-service

logs-quantum:
	docker-compose logs -f quantum-calculator

logs-auth:
	docker-compose logs -f auth-service

logs-db:
	docker-compose logs -f db

logs-redis:
	docker-compose logs -f redis

# Check health of all services
health:
	@echo "Checking service health..."
	@curl -s http://localhost:8080/health | python3 -m json.tool || echo "Gateway not responding"

# Run integration tests
test:
	@echo "Running integration tests..."
	python3 test_microservices.py

# Clean up everything
clean:
	docker-compose down -v --remove-orphans
	@echo "All containers and volumes removed"

# Development: rebuild and restart a specific service
rebuild-gateway:
	docker-compose up -d --build api-gateway

rebuild-compound:
	docker-compose up -d --build compound-service

rebuild-analog:
	docker-compose up -d --build analog-service

rebuild-ml:
	docker-compose up -d --build ml-service

rebuild-quantum:
	docker-compose up -d --build quantum-calculator

rebuild-auth:
	docker-compose up -d --build auth-service

# Database operations
db-shell:
	docker-compose exec db psql -U pharmasight_user -d pharmasight_db

redis-cli:
	docker-compose exec redis redis-cli

# Show running containers
ps:
	docker-compose ps

# Show resource usage
stats:
	docker stats --no-stream
