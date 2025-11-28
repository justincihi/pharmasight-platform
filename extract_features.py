#!/usr/bin/env python3
import subprocess
import json

branches = {
    "feature/pk-core-v2": "Most recent (Nov 18)",
    "replit-enhanced": "Replit version (Nov 12)", 
    "analog-discoveries-ip-protected": "IP protected (Nov 11)",
    "11/10": "Microservices (Nov 10)",
    "10-28": "Oct 28 version",
    "branch-10": "v2.1.0",
    "branch-13": "500+ compounds",
    "branch-14": "Advanced research",
    "branch-30": "Version 30"
}

features = {}

for branch, desc in branches.items():
    print(f"\n=== Analyzing {branch} ({desc}) ===")
    features[branch] = {
        "description": desc,
        "unique_files": [],
        "python_modules": [],
        "frontend_files": [],
        "images": [],
        "docs": []
    }
    
    # Get unique Python files
    cmd = f"git ls-tree -r origin/{branch} --name-only | grep '\\.py$' | grep -v venv | grep -v __pycache__"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    py_files = result.stdout.strip().split('\n') if result.stdout else []
    features[branch]["python_modules"] = [f for f in py_files if f and 'src/' in f][:10]
    
    # Get images
    cmd = f"git ls-tree -r origin/{branch} --name-only | grep -E '\\.(png|jpg|jpeg|svg)$' | grep -v venv"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    features[branch]["images"] = result.stdout.strip().split('\n') if result.stdout else []
    
    # Get docs
    cmd = f"git ls-tree -r origin/{branch} --name-only | grep -E '\\.md$' | head -20"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    features[branch]["docs"] = result.stdout.strip().split('\n') if result.stdout else []
    
    print(f"  Python modules: {len(features[branch]['python_modules'])}")
    print(f"  Images: {len(features[branch]['images'])}")
    print(f"  Docs: {len(features[branch]['docs'])}")

# Save results
with open('branch_features.json', 'w') as f:
    json.dump(features, f, indent=2)

print("\n✅ Feature extraction complete! Saved to branch_features.json")
