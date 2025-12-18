#!/usr/bin/env python3
"""
GitHub Push Script using Replit's GitHub Integration
Pushes changes to the specified branch on GitHub
"""

import os
import subprocess
import requests
import json
import base64
from datetime import datetime

def get_github_access_token():
    """Get GitHub access token from Replit's connector API"""
    hostname = os.environ.get('REPLIT_CONNECTORS_HOSTNAME')
    repl_identity = os.environ.get('REPL_IDENTITY')
    web_repl_renewal = os.environ.get('WEB_REPL_RENEWAL')
    
    if repl_identity:
        x_replit_token = f'repl {repl_identity}'
    elif web_repl_renewal:
        x_replit_token = f'depl {web_repl_renewal}'
    else:
        raise Exception('X_REPLIT_TOKEN not found for repl/depl')
    
    if not hostname:
        raise Exception('REPLIT_CONNECTORS_HOSTNAME not found')
    
    response = requests.get(
        f'https://{hostname}/api/v2/connection?include_secrets=true&connector_names=github',
        headers={
            'Accept': 'application/json',
            'X_REPLIT_TOKEN': x_replit_token
        }
    )
    
    data = response.json()
    connection = data.get('items', [{}])[0] if data.get('items') else {}
    settings = connection.get('settings', {})
    
    access_token = settings.get('access_token') or settings.get('oauth', {}).get('credentials', {}).get('access_token')
    
    if not access_token:
        raise Exception('GitHub not connected or access token not found')
    
    return access_token

def configure_git_credentials(token):
    """Configure git to use the access token"""
    subprocess.run(['git', 'config', '--global', 'credential.helper', 'store'], check=True)
    
    result = subprocess.run(['git', 'remote', 'get-url', 'origin'], capture_output=True, text=True)
    remote_url = result.stdout.strip()
    
    if 'github.com' in remote_url:
        if remote_url.startswith('https://'):
            parts = remote_url.replace('https://', '').split('/')
            if len(parts) >= 2:
                repo_path = '/'.join(parts[1:])
                new_url = f'https://oauth2:{token}@github.com/{repo_path}'
                subprocess.run(['git', 'remote', 'set-url', 'origin', new_url], check=True)
                return True
    return False

def push_to_github(branch='Replit-December'):
    """Push current branch to GitHub"""
    print(f"🚀 Pushing to GitHub branch: {branch}")
    print("=" * 50)
    
    try:
        print("1. Getting GitHub access token...")
        token = get_github_access_token()
        print("   ✅ Token retrieved successfully")
        
        print("2. Configuring git credentials...")
        configure_git_credentials(token)
        print("   ✅ Credentials configured")
        
        print(f"3. Pushing to origin/{branch}...")
        result = subprocess.run(
            ['git', 'push', 'origin', branch],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print("   ✅ Push successful!")
            print("\nOutput:")
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr)
            return True
        else:
            print(f"   ❌ Push failed with code {result.returncode}")
            print("Error:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Error: {e}")
        return False
    finally:
        result = subprocess.run(['git', 'remote', 'get-url', 'origin'], capture_output=True, text=True)
        original_url = result.stdout.strip()
        if '@github.com' in original_url:
            clean_url = original_url.split('@')[1] if '@' in original_url else original_url
            clean_url = f'https://{clean_url}'
            subprocess.run(['git', 'remote', 'set-url', 'origin', clean_url], capture_output=True)

if __name__ == '__main__':
    import sys
    branch = sys.argv[1] if len(sys.argv) > 1 else 'Replit-December'
    success = push_to_github(branch)
    sys.exit(0 if success else 1)
