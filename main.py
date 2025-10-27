#!/usr/bin/env python3
"""
PharmaSight™ Platform - Main Entry Point
Runs the complete enterprise application on port 5000
"""

import sys
sys.path.insert(0, 'src')

from pharmasight_complete import app

if __name__ == '__main__':
    print("=" * 60)
    print("🧬 PharmaSight™ Platform Starting...")
    print("=" * 60)
    print("📊 Enterprise Drug Discovery Platform")
    print("🔬 RDKit Integration: ENABLED")
    print("🧪 All Features: OPERATIONAL")
    print("🌐 Server: http://0.0.0.0:5000")
    print("=" * 60)
    app.run(host='0.0.0.0', port=5000, debug=False)
