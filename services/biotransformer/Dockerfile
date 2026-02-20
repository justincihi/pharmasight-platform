FROM openjdk:11-jre-slim

# Install Python and dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Download BioTransformer 3.0
WORKDIR /opt
RUN wget https://bitbucket.org/djoumbou/biotransformer/downloads/BioTransformer3.0.jar \
    -O BioTransformer3.0.jar || \
    echo "Note: BioTransformer download may require manual setup"

# Create app directory
WORKDIR /app

# Copy requirements and install Python packages
COPY requirements.txt /app/
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy application code
COPY wrapper.py /app/
COPY config.py /app/

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python3 -c "import requests; requests.get('http://localhost:8000/health')" || exit 1

# Run the application
CMD ["python3", "wrapper.py"]
