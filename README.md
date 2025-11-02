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

- Docker
- Docker Compose

### Running the Platform

1. **Build and Start Services:**
   From the root of the project directory, run the following command to build the Docker images and start all the services in detached mode:

   ```bash
   docker-compose up --build -d
   ```

2. **Accessing the Platform:**
   Once all services are running, the main API Gateway will be accessible at `http://localhost:8080`.

3. **Viewing Logs:**
   To view the logs for all running services, use:

   ```bash
   docker-compose logs -f
   ```

   To view logs for a specific service (e.g., `compound-service`):

   ```bash
   docker-compose logs -f compound-service
   ```

4. **Stopping the Platform:**
   To stop all running services, use the following command:

   ```bash
   docker-compose down
   ```

## 4. API Endpoints

All endpoints are accessed through the API Gateway. The gateway uses the path to route to the correct service.

**Example:** A `POST` request to `http://localhost:8080/compound-analysis/analyze` will be routed to the `/analyze` endpoint on the `compound-service`.

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
```
