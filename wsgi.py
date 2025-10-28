#!/usr/bin/env python3.11
"""
WSGI entry point for PharmaSight Platform
Production deployment with gunicorn
"""

import sys
import os

# Add src directory to Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import the Flask app
try:
    from pharmasight_complete import app
    print("✅ Successfully imported pharmasight_complete app")
except ImportError as e:
    print(f"❌ Failed to import pharmasight_complete: {e}")
    print("Falling back to basic app.py")
    sys.path.insert(0, os.path.dirname(__file__))
    from app import app

# Make app available for gunicorn
application = app

if __name__ == "__main__":
    # For direct execution (development)
    app.run(host='0.0.0.0', port=5000, debug=False)

