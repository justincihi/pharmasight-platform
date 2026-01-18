"""
PharmaSight™ Web Frontend Service
FastAPI server for serving static files and handling file uploads
"""

from fastapi import FastAPI, File, UploadFile, HTTPException, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import os
import shutil
from pathlib import Path
from typing import List
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="PharmaSight Web Frontend",
    description="Futuristic web interface for PharmaSight drug discovery platform",
    version="1.0.0"
)

# CORS configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Setup paths
BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"
TEMPLATES_DIR = BASE_DIR / "templates"
UPLOAD_DIR = STATIC_DIR / "assets" / "images"

# Ensure upload directory exists
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Mount static files
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Setup templates
templates = Jinja2Templates(directory=TEMPLATES_DIR)

# Allowed file extensions for uploads
ALLOWED_EXTENSIONS = {'.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp'}
MAX_FILE_SIZE = 10 * 1024 * 1024  # 10MB


# ===================================
# Routes
# ===================================

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Serve the main landing page"""
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "web-frontend",
        "version": "1.0.0"
    }


@app.get("/dashboard", response_class=HTMLResponse)
async def dashboard(request: Request):
    """Redirect to API gateway or serve dashboard"""
    # In production, this would redirect to the admin dashboard
    return JSONResponse(
        content={
            "message": "Dashboard access",
            "redirect": "http://localhost:8080/health"
        }
    )


# ===================================
# File Upload Endpoints
# ===================================

@app.post("/api/upload/image")
async def upload_image(file: UploadFile = File(...)):
    """
    Upload an image to the assets directory

    This endpoint allows admins to upload graphics for the About page
    and other sections of the site.
    """
    try:
        # Validate file extension
        file_ext = Path(file.filename).suffix.lower()
        if file_ext not in ALLOWED_EXTENSIONS:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid file type. Allowed: {', '.join(ALLOWED_EXTENSIONS)}"
            )

        # Read file content
        content = await file.read()

        # Validate file size
        if len(content) > MAX_FILE_SIZE:
            raise HTTPException(
                status_code=400,
                detail=f"File too large. Maximum size: {MAX_FILE_SIZE / 1024 / 1024}MB"
            )

        # Generate safe filename
        safe_filename = "".join(c for c in file.filename if c.isalnum() or c in ('_', '-', '.'))
        file_path = UPLOAD_DIR / safe_filename

        # Save file
        with open(file_path, "wb") as f:
            f.write(content)

        logger.info(f"Image uploaded successfully: {safe_filename}")

        return {
            "success": True,
            "filename": safe_filename,
            "url": f"/static/assets/images/{safe_filename}",
            "size": len(content)
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Upload error: {str(e)}")
        raise HTTPException(status_code=500, detail="Upload failed")


@app.post("/api/upload/bulk")
async def upload_bulk_images(files: List[UploadFile] = File(...)):
    """Upload multiple images at once"""
    if len(files) > 10:
        raise HTTPException(
            status_code=400,
            detail="Maximum 10 files per upload"
        )

    results = []
    errors = []

    for file in files:
        try:
            # Validate and save each file
            file_ext = Path(file.filename).suffix.lower()
            if file_ext not in ALLOWED_EXTENSIONS:
                errors.append({
                    "filename": file.filename,
                    "error": "Invalid file type"
                })
                continue

            content = await file.read()

            if len(content) > MAX_FILE_SIZE:
                errors.append({
                    "filename": file.filename,
                    "error": "File too large"
                })
                continue

            safe_filename = "".join(c for c in file.filename if c.isalnum() or c in ('_', '-', '.'))
            file_path = UPLOAD_DIR / safe_filename

            with open(file_path, "wb") as f:
                f.write(content)

            results.append({
                "filename": safe_filename,
                "url": f"/static/assets/images/{safe_filename}",
                "size": len(content)
            })

        except Exception as e:
            errors.append({
                "filename": file.filename,
                "error": str(e)
            })

    return {
        "success": len(results) > 0,
        "uploaded": results,
        "errors": errors,
        "total": len(results)
    }


@app.get("/api/images")
async def list_images():
    """List all uploaded images"""
    try:
        images = []
        for file_path in UPLOAD_DIR.iterdir():
            if file_path.is_file() and file_path.suffix.lower() in ALLOWED_EXTENSIONS:
                stat = file_path.stat()
                images.append({
                    "filename": file_path.name,
                    "url": f"/static/assets/images/{file_path.name}",
                    "size": stat.st_size,
                    "modified": stat.st_mtime
                })

        # Sort by modification time (newest first)
        images.sort(key=lambda x: x['modified'], reverse=True)

        return {
            "images": images,
            "total": len(images)
        }

    except Exception as e:
        logger.error(f"Error listing images: {str(e)}")
        raise HTTPException(status_code=500, detail="Failed to list images")


@app.delete("/api/images/{filename}")
async def delete_image(filename: str):
    """Delete an uploaded image"""
    try:
        # Sanitize filename
        safe_filename = "".join(c for c in filename if c.isalnum() or c in ('_', '-', '.'))
        file_path = UPLOAD_DIR / safe_filename

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="Image not found")

        file_path.unlink()
        logger.info(f"Image deleted: {safe_filename}")

        return {
            "success": True,
            "message": f"Image {safe_filename} deleted successfully"
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Delete error: {str(e)}")
        raise HTTPException(status_code=500, detail="Failed to delete image")


# ===================================
# Integrations Documentation API
# ===================================

@app.get("/api/integrations")
async def get_integrations():
    """Return current and planned integrations"""
    return {
        "current": {
            "cheminformatics": [
                {"name": "RDKit", "description": "Open-source cheminformatics toolkit", "version": "2025.09.1"},
                {"name": "OpenBabel", "description": "Chemical file format converter", "version": "3.1"},
                {"name": "ChemAxon", "description": "Advanced chemistry tools", "version": "enterprise"}
            ],
            "databases": [
                {"name": "PubChem", "description": "Public chemical database", "compounds": "110M+"},
                {"name": "ChEMBL", "description": "Bioactivity database", "compounds": "2.3M+"},
                {"name": "ZINC", "description": "Drug-like compounds", "compounds": "230M+"}
            ],
            "molecular_dynamics": [
                {"name": "AutoDock Vina", "description": "Molecular docking", "version": "1.2"},
                {"name": "GROMACS", "description": "MD simulations", "version": "2024"},
                {"name": "BioTransformer", "description": "Metabolite prediction", "version": "3.0"}
            ],
            "quantum_computing": [
                {"name": "PySCF", "description": "Quantum chemistry", "version": "2.4"},
                {"name": "Qiskit", "description": "Quantum computing", "version": "1.0"},
                {"name": "Psi4", "description": "Quantum chemistry", "version": "1.9"}
            ],
            "ai_ml": [
                {"name": "OpenAI GPT-4", "description": "Generative AI", "status": "integrated"},
                {"name": "Google Gemini", "description": "Multimodal AI", "status": "integrated"},
                {"name": "Anthropic Claude", "description": "Constitutional AI", "status": "integrated"}
            ],
            "visualization": [
                {"name": "3Dmol.js", "description": "3D molecular viewer", "version": "2.0"},
                {"name": "Plotly", "description": "Interactive charts", "version": "5.0"},
                {"name": "BioRender", "description": "Scientific illustrations", "status": "planned"}
            ]
        },
        "planned": [
            {
                "name": "CRISPR Design Tools",
                "category": "Genetics",
                "description": "Genetic engineering capabilities for target validation",
                "eta": "Q2 2026"
            },
            {
                "name": "IBM RXN for Chemistry",
                "category": "Retrosynthesis",
                "description": "AI-powered retrosynthetic analysis",
                "eta": "Q1 2026"
            },
            {
                "name": "Flow Chemistry Platform",
                "category": "Synthesis",
                "description": "Automated synthesis workflow optimization",
                "eta": "Q3 2026"
            },
            {
                "name": "NVIDIA bioNEMO",
                "category": "Generative AI",
                "description": "Generative AI for biomolecules (pending licensing)",
                "eta": "Q2 2026"
            },
            {
                "name": "ASKCOS",
                "category": "Retrosynthesis",
                "description": "MIT's computer-aided synthesis planning",
                "eta": "Q2 2026"
            },
            {
                "name": "VCell/COPASI",
                "category": "Biosimulation",
                "description": "Open-source biological pathway simulation",
                "eta": "Q3 2026"
            },
            {
                "name": "PopHive",
                "category": "Population Health",
                "description": "Population-level drug response modeling",
                "eta": "Q4 2026"
            }
        ]
    }


# ===================================
# Main Entry Point
# ===================================

if __name__ == "__main__":
    import uvicorn

    port = int(os.getenv("PORT", 8090))

    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=port,
        reload=True,
        log_level="info"
    )
