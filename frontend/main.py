from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware
import os

app = FastAPI(title="PharmaSight Frontend", version="1.0.0")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")

# Templates
templates = Jinja2Templates(directory="templates")

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Serve the main frontend page"""
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/health")
async def health():
    """Health check endpoint"""
    return {"status": "healthy", "service": "frontend"}

@app.get("/PUBLIC_ANALOG_DISCOVERY_REGISTRY.md")
async def registry():
    """Serve the public IP registry"""
    registry_path = "../PUBLIC_ANALOG_DISCOVERY_REGISTRY.md"
    if os.path.exists(registry_path):
        return FileResponse(registry_path, media_type="text/markdown")
    return {"error": "Registry not found"}

@app.get("/docs-page", response_class=HTMLResponse)
async def docs_page(request: Request):
    """API documentation page"""
    return """
    <html>
        <head>
            <title>PharmaSight API Documentation</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; background: #0a0e27; color: #fff; }
                h1 { color: #00a8e8; }
                .endpoint { background: #1a1f3a; padding: 20px; margin: 20px 0; border-radius: 8px; }
                code { background: #0a0e27; padding: 4px 8px; border-radius: 4px; color: #00d9ff; }
            </style>
        </head>
        <body>
            <h1>PharmaSight™ API Documentation</h1>
            <div class="endpoint">
                <h2>Research Engine API</h2>
                <p><code>GET /research-engine/articles/all</code> - Get all research articles</p>
                <p><code>POST /research-engine/articles/search</code> - Search articles</p>
                <p><code>GET /research-engine/discoveries/all</code> - Get all analog discoveries</p>
                <p><code>POST /research-engine/research/run-cycle</code> - Run research cycle</p>
            </div>
            <div class="endpoint">
                <h2>Compound Analysis API</h2>
                <p><code>POST /compound-service/analyze</code> - Analyze compound</p>
            </div>
            <div class="endpoint">
                <h2>Analog Generation API</h2>
                <p><code>POST /analog-service/generate</code> - Generate analogs</p>
            </div>
            <p><a href="/" style="color: #00a8e8;">← Back to Platform</a></p>
        </body>
    </html>
    """

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=3000)

