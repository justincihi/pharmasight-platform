#!/usr/bin/env python3
"""
PharmaSight Platform - Entry Point for Deployment
This file imports and runs the main Flask application from src/pharmasight_complete.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

from pharmasight_complete import app

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=False)
