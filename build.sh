#!/bin/bash
# Build script for Render deployment

set -e  # Exit on error

echo "Starting PharmaSight Platform build..."

# Install Python dependencies
echo "Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt
pip install -r requirements-rdkit.txt

echo "Build completed successfully!"
