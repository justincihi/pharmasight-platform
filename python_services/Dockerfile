FROM python:3.11-slim

# Install system dependencies for RDKit and AutoDock Vina
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    autodock-vina \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY analysis_microservice.py .
COPY ../server/python_modules ./python_modules

# Create non-root user
RUN useradd -m -u 1000 pharmasight && chown -R pharmasight:pharmasight /app
USER pharmasight

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import requests; requests.get('http://localhost:5000/health')"

# Run with gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "4", "--timeout", "120", "analysis_microservice:app"]
