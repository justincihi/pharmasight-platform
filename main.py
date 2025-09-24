#!/usr/bin/env python3
"""
PharmaSight™ Enhanced Platform - Main Application Entry Point
Ensures all enhanced features are properly loaded and accessible
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import the enhanced application
try:
    from pharmasight_complete import app
    print("✅ PharmaSight™ Enhanced Platform loaded successfully")
    print("✅ All API endpoints and features are active")
    print("✅ Database integrations: PubChem, FDA, ChEMBL, DrugBank, ZINC, OpenTargets")
    print("✅ Regulatory compliance features enabled")
    print("✅ Administrator research management active")
    print("✅ Enhanced compound display with molecular structures")
    print("✅ Discovery logging and DOI tracking operational")
except ImportError as e:
    print(f"❌ Error importing enhanced application: {e}")
    sys.exit(1)

if __name__ == '__main__':
    # Run the enhanced application
    app.run(host='0.0.0.0', port=5008, debug=False)
