# PharmaSight™ Ultimate - Deployment Entry Point
# This file properly imports the ultimate platform for deployment

import sys
import os

# Add parent directory to Python path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

# Import the Flask app from the main module in parent directory
import importlib.util
spec = importlib.util.spec_from_file_location("main_app", os.path.join(parent_dir, "main.py"))
main_app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(main_app)

# Get the app instance
app = main_app.app

if __name__ == '__main__':
    print("🚀 Starting PharmaSight™ Ultimate Platform...")
    print("✅ Version: 4.0.0-ULTIMATE")
    print("✅ Features: All advanced research capabilities enabled")
    app.run(host='0.0.0.0', port=8080, debug=False)
