from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import httpx
import logging
import redis
import os
from datetime import datetime, timedelta

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="PharmaSight API Gateway", version="1.0.0")

# CORS - allow the frontend to call the gateway (restrict in production)
app.add_middleware(
    CORSMiddleware,
    allow_origins=os.getenv("CORS_ALLOW_ORIGINS", "*").split(",") if os.getenv("CORS_ALLOW_ORIGINS") else ["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

SERVICE_URLS = {
    "compound-service": "http://compound-service:8000",
    "analog-service": "http://analog-service:8000",
    "ml-service": "http://ml-service:8000",
    "quantum-calculator": "http://quantum-calculator:8000",
    "auth-service": "http://auth-service:8000",
}

REQUEST_TIMEOUT = 30.0  # seconds

# Redis connection for rate limiting (optional)
try:
    redis_client = redis.Redis(
        host=os.getenv("REDIS_HOST", "redis"),
        port=int(os.getenv("REDIS_PORT", 6379)),
        db=0,
        decode_responses=True,
    )
    # quick check
    redis_client.ping()
except Exception as e:
    logger.warning(f"Redis unavailable, rate limiting disabled: {e}")
    redis_client = None

# Rate limiting configuration
RATE_LIMIT_REQUESTS = int(os.getenv("RATE_LIMIT_REQUESTS", 100))  # requests per window
RATE_LIMIT_WINDOW = int(os.getenv("RATE_LIMIT_WINDOW", 60))  # seconds

def check_rate_limit(client_ip: str) -> bool:
    """Check if client IP is within rate limits."""
    # If Redis is unavailable, allow requests (no rate limiting)
    if not redis_client:
        return True
    key = f"rate_limit:{client_ip}"
    current_time = datetime.utcnow().timestamp()

    # Use Redis pipeline for atomic operations
    pipe = redis_client.pipeline()
    pipe.zremrangebyscore(key, 0, current_time - RATE_LIMIT_WINDOW)
    pipe.zadd(key, {str(current_time): current_time})
    pipe.zcard(key)
    pipe.expire(key, RATE_LIMIT_WINDOW)
    results = pipe.execute()

    request_count = results[2]
    return request_count <= RATE_LIMIT_REQUESTS

@app.get("/health")
async def health_check():
    """Health check endpoint for the API Gateway."""
    service_health = {}
    async with httpx.AsyncClient(timeout=5.0) as client:
        for service_name, service_url in SERVICE_URLS.items():
            try:
                response = await client.get(f"{service_url}/health", timeout=5.0)
                service_health[service_name] = {
                    "status": "healthy" if response.status_code == 200 else "unhealthy",
                    "status_code": response.status_code
                }
            except Exception as e:
                service_health[service_name] = {
                    "status": "unreachable",
                    "error": str(e)
                }
    
    all_healthy = all(s["status"] == "healthy" for s in service_health.values())
    return {
        "gateway": "healthy",
        "services": service_health,
        "overall_status": "healthy" if all_healthy else "degraded"
    }

@app.api_route("/{service}/{path:path}", methods=["GET", "POST", "PUT", "DELETE"])
async def route_request(service: str, path: str, request: Request):
    if service not in SERVICE_URLS:
        logger.error(f"Service not found: {service}")
        raise HTTPException(status_code=404, detail=f"Service '{service}' not found")

    # Rate limiting check
    client_ip = request.client.host if request.client else "unknown"
    if not check_rate_limit(client_ip):
        logger.warning(f"Rate limit exceeded for IP: {client_ip}")
        raise HTTPException(status_code=429, detail="Too many requests")

    service_url = f"{SERVICE_URLS[service]}/{path}"
    logger.info(f"Routing {request.method} request to {service_url}")
    
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT) as client:
        method = request.method
        headers = dict(request.headers)
        # The host header is for the gateway, not the downstream service
        headers.pop("host", None)
        
        # Read raw body (may be empty). We'll try JSON first, otherwise forward raw bytes.
        body_bytes = await request.body()
        json_data = None
        if body_bytes:
            try:
                json_data = await request.json()
            except Exception:
                json_data = None

        try:
            if json_data is not None:
                response = await client.request(
                    method,
                    service_url,
                    headers=headers,
                    json=json_data,
                    params=dict(request.query_params),
                    timeout=REQUEST_TIMEOUT,
                )
            else:
                response = await client.request(
                    method,
                    service_url,
                    headers=headers,
                    content=body_bytes if body_bytes else None,
                    params=dict(request.query_params),
                    timeout=REQUEST_TIMEOUT,
                )

            # Preserve downstream status code and try to return JSON if available
            content_type = response.headers.get("content-type", "application/json")
            try:
                parsed = response.json()
                return JSONResponse(content=parsed, status_code=response.status_code)
            except Exception:
                logger.warning(f"Non-JSON response from {service_url}")
                return Response(content=response.content, status_code=response.status_code, media_type=content_type)
        except httpx.TimeoutException:
            logger.error(f"Timeout when calling {service_url}")
            raise HTTPException(status_code=504, detail=f"Service timeout: {service}")
        except httpx.RequestError as exc:
            logger.error(f"Request error when calling {service_url}: {exc}")
            raise HTTPException(status_code=503, detail=f"Service unavailable: {service}")

