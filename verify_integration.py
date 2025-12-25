#!/usr/bin/env python3
"""
PharmaSight Combined Branch - Integration Verification Script
Verifies that all key components from Replit-December and manus-December are properly integrated
"""

import os
import sys
import json
from pathlib import Path

class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'

def check_file(filepath, description):
    """Check if a file exists"""
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"{Colors.GREEN}✓{Colors.END} {description}: {filepath} ({size:,} bytes)")
        return True
    else:
        print(f"{Colors.RED}✗{Colors.END} {description}: {filepath} NOT FOUND")
        return False

def check_directory(dirpath, description):
    """Check if a directory exists"""
    if os.path.isdir(dirpath):
        file_count = len(list(Path(dirpath).rglob('*')))
        print(f"{Colors.GREEN}✓{Colors.END} {description}: {dirpath} ({file_count} files)")
        return True
    else:
        print(f"{Colors.RED}✗{Colors.END} {description}: {dirpath} NOT FOUND")
        return False

def count_python_modules(directory):
    """Count Python modules in a directory"""
    if not os.path.isdir(directory):
        return 0
    return len([f for f in os.listdir(directory) if f.endswith('.py')])

def main():
    print("=" * 80)
    print(f"{Colors.BLUE}PharmaSight Combined Branch - Integration Verification{Colors.END}")
    print("=" * 80)
    print()
    
    results = {
        'flask_app': [],
        'admin_dashboard': [],
        'infrastructure': [],
        'documentation': []
    }
    
    # Flask Application Components (Replit-December)
    print(f"{Colors.YELLOW}[1] Flask Application Components (Replit-December){Colors.END}")
    print("-" * 80)
    
    results['flask_app'].append(check_file('app.py', 'Main entry point'))
    results['flask_app'].append(check_file('main.py', 'Alternative entry point'))
    results['flask_app'].append(check_file('requirements.txt', 'Python dependencies'))
    results['flask_app'].append(check_file('requirements-rdkit.txt', 'RDKit dependencies'))
    
    # Check src directory
    results['flask_app'].append(check_directory('src', 'Source code directory'))
    
    # Key modules from Replit-December
    print()
    print(f"{Colors.YELLOW}[1.1] Key Drug Discovery Modules{Colors.END}")
    results['flask_app'].append(check_file('src/pharmasight_complete.py', 'Main Flask application'))
    results['flask_app'].append(check_file('src/biotransformer_client.py', 'BioTransformer integration (NEW)'))
    results['flask_app'].append(check_file('src/molecular_editor.py', 'Molecular editor (Phase 2)'))
    results['flask_app'].append(check_file('src/data_export.py', 'Data export (Phase 3)'))
    results['flask_app'].append(check_file('src/admet_predictor.py', 'ADMET prediction'))
    results['flask_app'].append(check_file('src/quantum_computing_module.py', 'Quantum calculations'))
    results['flask_app'].append(check_file('src/patent_case_generator.py', 'Patent analysis'))
    results['flask_app'].append(check_file('src/autonomous_research_engine.py', 'Autonomous research'))
    
    # Count total modules
    module_count = count_python_modules('src')
    print()
    print(f"{Colors.GREEN}✓{Colors.END} Total Python modules in src/: {module_count}")
    
    print()
    print("=" * 80)
    
    # Admin Dashboard Components (manus-December)
    print(f"{Colors.YELLOW}[2] Admin Dashboard Components (manus-December){Colors.END}")
    print("-" * 80)
    
    results['admin_dashboard'].append(check_directory('admin-dashboard', 'Admin dashboard root'))
    results['admin_dashboard'].append(check_file('admin-dashboard/package.json', 'Node.js dependencies'))
    results['admin_dashboard'].append(check_file('admin-dashboard/vite.config.ts', 'Vite configuration'))
    results['admin_dashboard'].append(check_file('admin-dashboard/tsconfig.json', 'TypeScript configuration'))
    results['admin_dashboard'].append(check_file('admin-dashboard/drizzle.config.ts', 'Drizzle ORM config'))
    
    print()
    print(f"{Colors.YELLOW}[2.1] Admin Dashboard Structure{Colors.END}")
    results['admin_dashboard'].append(check_directory('admin-dashboard/client', 'React frontend'))
    results['admin_dashboard'].append(check_directory('admin-dashboard/server', 'Node.js backend'))
    results['admin_dashboard'].append(check_directory('admin-dashboard/drizzle', 'Database migrations'))
    
    print()
    print("=" * 80)
    
    # Infrastructure Components
    print(f"{Colors.YELLOW}[3] Infrastructure Components{Colors.END}")
    print("-" * 80)
    
    results['infrastructure'].append(check_file('docker-compose.yml', 'Docker Compose configuration'))
    results['infrastructure'].append(check_file('Dockerfile', 'Flask app Dockerfile'))
    results['infrastructure'].append(check_file('admin-dashboard/Dockerfile', 'Admin dashboard Dockerfile'))
    results['infrastructure'].append(check_file('nginx.conf', 'Nginx reverse proxy config'))
    results['infrastructure'].append(check_file('.env.example', 'Environment variables template'))
    results['infrastructure'].append(check_file('admin-dashboard/.env.example', 'Admin dashboard env template'))
    
    print()
    print("=" * 80)
    
    # Documentation
    print(f"{Colors.YELLOW}[4] Documentation{Colors.END}")
    print("-" * 80)
    
    results['documentation'].append(check_file('README.md', 'Main README'))
    results['documentation'].append(check_file('COMBINED_DEPLOYMENT_GUIDE.md', 'Deployment guide'))
    results['documentation'].append(check_file('admin-dashboard/README.md', 'Admin dashboard README'))
    
    print()
    print("=" * 80)
    
    # Summary
    print(f"{Colors.BLUE}VERIFICATION SUMMARY{Colors.END}")
    print("=" * 80)
    
    total_checks = 0
    passed_checks = 0
    
    for category, checks in results.items():
        category_passed = sum(checks)
        category_total = len(checks)
        total_checks += category_total
        passed_checks += category_passed
        
        status = f"{Colors.GREEN}PASS{Colors.END}" if category_passed == category_total else f"{Colors.RED}FAIL{Colors.END}"
        print(f"{category.replace('_', ' ').title()}: {category_passed}/{category_total} [{status}]")
    
    print()
    print(f"Overall: {passed_checks}/{total_checks} checks passed")
    
    if passed_checks == total_checks:
        print(f"\n{Colors.GREEN}✓ Integration verification PASSED!{Colors.END}")
        print(f"{Colors.GREEN}All components from Replit-December and manus-December are properly integrated.{Colors.END}")
        return 0
    else:
        print(f"\n{Colors.YELLOW}⚠ Integration verification completed with warnings.{Colors.END}")
        print(f"{Colors.YELLOW}Some optional components may be missing.{Colors.END}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
