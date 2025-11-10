# PharmaSight™ - Microservices Architecture

This document outlines the new microservices-based architecture for the PharmaSight™ platform.

## 1. Overview

The platform has been refactored from a monolithic Flask application into a distributed system of specialized microservices, orchestrated with Docker Compose. This modern architecture provides enhanced scalability, maintainability, and performance. An API Gateway acts as the single entry point for all client requests, routing them to the appropriate backend service.

## 2. Services

The architecture consists of the following services:

| Service                | Port | Description                                                                                             |
| ---------------------- | ---- | ------------------------------------------------------------------------------------------------------- |
| **API Gateway**        | 8080 | The single entry point for all client requests. Routes traffic to the appropriate downstream service.     |
| **Authentication**     | 8005 | Manages user authentication (JWT), authorization, and role-based access control (RBAC).                   |
| **Compound Analysis**  | 8001 | Handles molecular analysis, property calculation (RDKit), and chemical structure visualization.           |
| **Analog Generation**  | 8002 | Generates and analyzes chemical analogs, including similarity scoring and patent analysis.                |
| **ML Models**          | 8003 | Serves pre-trained machine learning models for predicting ADMET properties, toxicity, and more.           |
| **Quantum Calculator** | 8004 | Performs high-precision quantum chemistry calculations (DFT) using PySCF for advanced property prediction. |
| **PostgreSQL Database**| 5432 | Primary relational database for storing compound information, user data, and research findings.         |
| **Redis Cache**        | 6379 | In-memory data store for caching expensive calculations and managing session data.                   |

## 3. Getting Started

### Prerequisites

- Docker (version 20.10+)
- Docker Compose (version 1.29+)
- Python 3.11+ (for running tests)

### Quick Start

1. **Clone the Repository:**

```bash
git clone https://github.com/justincihi/pharmasight-platform.git
cd pharmasight-platform
```

2. **Configure Environment (Optional):**

```bash
cp .env.example .env
# Edit .env with your preferred settings
```

3. **Build and Start Services:**
From the root of the project directory, run:

```bash
docker-compose up --build -d
```

Or use the Makefile:

```bash
make up
```

4. **Verify Services Are Running:**
Check the health of all services:

```bash
curl http://localhost:8080/health
```

Or use:

```bash
make health
```

5. **Run Integration Tests:**

```bash
python3 test_microservices.py
```

Or use:

```bash
make test
```

### Accessing the Platform

Once all services are running:

- **API Gateway:** [<http://localhost:8080>](http://localhost:8080)
- **Interactive API Docs:** [<http://localhost:8080/docs>](<http://localhost:8080/docs>)
- **Health Check:** [<http://localhost:8080/health>](<http://localhost:8080/health>)

### Common Operations

**View Logs:**

```bash
docker-compose logs -f                    # All services
docker-compose logs -f compound-service   # Specific service
make logs                                 # Using Makefile
make logs-compound                        # Specific service
```

**Stop Services:**

```bash
docker-compose down
make down
```

**Restart Services:**

```bash
docker-compose restart
make restart
```

**Clean Everything:**

```bash
docker-compose down -v
make clean
```

## 4. API Endpoints

All endpoints are accessed through the API Gateway. The gateway uses the path to route to the correct service.

**Example:** A `POST` request to `<http://localhost:8080/compound-analysis/analyze>` will be routed to the `/analyze` endpoint on the `compound-service`.

- **Authentication:**
  - `POST /auth-service/token`: Obtain a JWT access token.
  - `GET /auth-service/users/me`: Get information about the currently authenticated user.

- **Compound Analysis:**
  - `POST /compound-analysis/analyze`: Analyze a compound by name or SMILES string.

- **Analog Generation:**
  - `POST /analog-generation/generate`: Generate analogs for a parent compound.

- **ML Predictions:**
  - `POST /ml-models/predict/admet`: Predict ADMET properties for a list of molecules.
  - `POST /ml-models/predict/properties`: Predict other properties for a single molecule.

- **Quantum Calculations:**
  - `POST /quantum-calculator/calculate`: Perform a DFT calculation on a molecule.

## 5. Development

When developing a specific service, you can run just that service and its dependencies. For example, to work on the `compound-service`:

```bash
docker-compose up -d db redis compound-service
```

This will start the database, cache, and the compound service, allowing you to test it in isolation. Remember to rebuild the image if you make changes to the code:

```bash
docker-compose up -d --build compound-service
# Or using Makefile
make rebuild-compound
```

### Development Workflow

1. **Make code changes** in `services/[service-name]/`
2. **Rebuild the service:** `make rebuild-[service-name]`
3. **Check logs:** `make logs-[service-name]`
4. **Test changes:** Run integration tests or use curl

### Additional Resources

- **Comprehensive API Documentation:** See [API_DOCUMENTATION.md](API_DOCUMENTATION.md)
- **Deployment Guide:** See [MICROSERVICES_DEPLOYMENT.md](MICROSERVICES_DEPLOYMENT.md)
- **Environment Configuration:** See [.env.example](.env.example)

### Troubleshooting

See the [MICROSERVICES_DEPLOYMENT.md](MICROSERVICES_DEPLOYMENT.md) file for detailed troubleshooting steps.

Common issues:
- **Port conflicts:** Change ports in docker-compose.yml
- **Services won't start:** Check logs with `make logs`
- **Database connection errors:** Verify DATABASE_URL environment variable
- **Health check failures:** Ensure curl is installed in containers
