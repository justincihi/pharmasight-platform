#!/bin/bash
# Automatically create Python venv and install dependencies if missing

VENV_DIR="$(dirname "$0")/venv"

if [ ! -d "$VENV_DIR" ]; then
  echo "[Python Setup] Creating virtual environment..."
  python3 -m venv "$VENV_DIR"
  
  echo "[Python Setup] Installing dependencies..."
  "$VENV_DIR/bin/pip" install --quiet rdkit scipy numpy google-generativeai
  
  echo "[Python Setup] ✅ Python environment ready"
else
  echo "[Python Setup] ✅ Virtual environment already exists"
fi
