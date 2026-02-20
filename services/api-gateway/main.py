from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
import httpx
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="PharmaSight API Gateway", version="1.0.0")

SERVICE_URLS = {
    "compound-analysis": "http://compound-service:8000",
    "analog-generation": "http://analog-service:8000",
    "ml-models": "http://ml-service:8000",
    "quantum-calculator": "http://quantum-calculator:8000",
    "auth-service": "http://auth-service:8000",
}

REQUEST_TIMEOUT = 30.0  # seconds

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

    service_url = f"{SERVICE_URLS[service]}/{path}"
    logger.info(f"Routing {request.method} request to {service_url}")
    
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT) as client:
        method = request.method
        headers = dict(request.headers)
        # The host header is for the gateway, not the downstream service
        headers.pop("host", None)
        
        json_data = None
        try:
            json_data = await request.json()
        except Exception:
            pass

        try:
            if json_data:
                response = await client.request(
                    method, service_url, headers=headers, 
                    json=json_data, params=request.query_params,
                    timeout=REQUEST_TIMEOUT
                )
            else:
                response = await client.request(
                    method, service_url, headers=headers, 
                    params=request.query_params,
                    timeout=REQUEST_TIMEOUT
                )
            
            # Handle different response types
            try:
                return response.json()
            except ValueError as e:
                # Response is not JSON
                logger.warning(f"Non-JSON response from {service_url}: {str(e)}")
                return JSONResponse(
                    content={"data": response.text, "warning": "Response was not JSON"},
                    status_code=response.status_code
                )
        except httpx.TimeoutException:
            logger.error(f"Timeout when calling {service_url}")
            raise HTTPException(status_code=504, detail=f"Service timeout: {service}")
        except httpx.RequestError as exc:
            logger.error(f"Request error when calling {service_url}: {exc}")
            raise HTTPException(status_code=503, detail=f"Service unavailable: {service}")

