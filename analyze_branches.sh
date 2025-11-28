#!/bin/bash
echo "# PharmaSight Branch Analysis Report"
echo "Generated: $(date)"
echo ""

branches=(
  "10-28"
  "11/10"
  "analog-discoveries-ip-protected"
  "branch-10"
  "branch-13"
  "branch-14"
  "branch-30"
  "branch-5"
  "branch-8"
  "claude/create-readme-011CUSE9T8QRX6nD8beW9cAT"
  "claude/merge-all-features-011CUSE9T8QRX6nD8beW9cAT"
  "copilot/fix-end-of-file-error"
  "copilot/fixscipy-compat"
  "feature/pk-core-v2"
  "main"
  "main-enhanced"
  "replit-enhanced"
)

for branch in "${branches[@]}"; do
  echo "## Branch: $branch"
  echo ""
  
  # Get last commit info
  echo "**Last Commit:**"
  git log origin/$branch -1 --format="- Date: %ai%n- Author: %an%n- Message: %s" 2>/dev/null || echo "Error accessing branch"
  echo ""
  
  # Count files
  echo "**File Count:**"
  git ls-tree -r origin/$branch --name-only 2>/dev/null | wc -l
  echo ""
  
  # List key directories
  echo "**Key Directories:**"
  git ls-tree -r origin/$branch --name-only 2>/dev/null | cut -d/ -f1 | sort -u | head -10
  echo ""
  echo "---"
  echo ""
done
