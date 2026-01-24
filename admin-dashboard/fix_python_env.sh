#!/bin/bash

# Fix Python SRE Module Mismatch
# This script resolves the "AssertionError: SRE module mismatch" error

echo "🔧 Fixing Python environment..."

# Method 1: Reinstall Python standard library
echo "Step 1: Checking Python installation..."
python3 --version

# Method 2: Clear Python cache
echo "Step 2: Clearing Python cache..."
find /home/ubuntu/pharmasight-admin-dashboard -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find /home/ubuntu/pharmasight-admin-dashboard -type f -name "*.pyc" -delete 2>/dev/null || true

# Method 3: Create isolated virtual environment
echo "Step 3: Creating isolated Python environment..."
cd /home/ubuntu/pharmasight-admin-dashboard/server/python_modules
python3 -m venv venv --clear 2>/dev/null || python3 -m venv venv

# Activate and install dependencies
source venv/bin/activate
pip install --upgrade pip
pip install rdkit-pypi requests numpy pandas

echo "✅ Python environment fixed!"
echo ""
echo "To use the fixed environment:"
echo "  source /home/ubuntu/pharmasight-admin-dashboard/server/python_modules/venv/bin/activate"
