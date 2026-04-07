#!/bin/bash
set -e

echo "Setting up Python virtual environment for PharmaSight..."

# Create venv directory if it doesn't exist
mkdir -p server/python_modules/venv

# Create virtual environment with Python 3.11
python3.11 -m venv server/python_modules/venv

# Activate venv
source server/python_modules/venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required packages
echo "Installing RDKit and dependencies..."
pip install rdkit scipy numpy scikit-learn pandas google-generativeai

echo "✅ Python environment setup complete!"
echo "Installed packages:"
pip list | grep -E "rdkit|scipy|numpy|scikit-learn|pandas|google-generativeai"
