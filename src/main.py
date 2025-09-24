# PharmaSight™ Ultimate - Production Entry Point
# This file ensures the ultimate platform is deployed correctly

import sys
import os

# Add parent directory to path to access main.py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the ultimate platform
from main import app

if __name__ == '__main__':
    print("🚀 Starting PharmaSight™ Ultimate Platform...")
    print("✅ Version: 4.0.0-ULTIMATE")
    print("✅ Features: All advanced research capabilities enabled")
    app.run(host='0.0.0.0', port=8080, debug=False)
