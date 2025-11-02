from fastapi import FastAPI, HTTPException, Request
import httpx

app = FastAPI()

SERVICE_URLS = {
    "compound-analysis": "http://compound-service:8000",
    "analog-generation": "http://analog-service:8000",
    "ml-models": "http://ml-service:8000",
    "quantum-calculator": "http://quantum-calculator:8000",
    "auth-service": "http://auth-service:8000",
}

@app.api_route("/{service}/{path:path}", methods=["GET", "POST", "PUT", "DELETE"])
async def route_request(service: str, path: str, request: Request):
    if service not in SERVICE_URLS:
        raise HTTPException(status_code=404, detail="Service not found")

    service_url = f"{SERVICE_URLS[service]}/{path}"
    
    async with httpx.AsyncClient() as client:
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
                response = await client.request(method, service_url, headers=headers, json=json_data, params=request.query_params)
            else:
                response = await client.request(method, service_url, headers=headers, params=request.query_params)
            
            return response.json()
        except httpx.RequestError as exc:
            raise HTTPException(status_code=503, detail=f"Service unavailable: {exc}")

