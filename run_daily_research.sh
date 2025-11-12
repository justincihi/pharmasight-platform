#!/bin/bash
#
# PharmaSight™ Daily Autonomous Research Script
# Run this script once per day to perform automated research activities
#
# Usage: ./run_daily_research.sh
# Or schedule with cron: 0 9 * * * /home/ubuntu/pharmasight-latest/run_daily_research.sh
#

echo "=========================================="
echo "PharmaSight™ Daily Research - $(date)"
echo "=========================================="

cd /home/ubuntu/pharmasight-latest

# Run autonomous research engine
echo "Running autonomous research engine..."
python3.11 src/autonomous_research_engine.py

# Update public README
echo "Updating public article database README..."
python3.11 -c "from src.research_article_database import ResearchArticleDatabase; db = ResearchArticleDatabase(); db.generate_public_readme()"

# Optional: Commit to GitHub (uncomment if you want automatic commits)
# echo "Committing updates to GitHub..."
# git add RESEARCH_ARTICLES_DATABASE.* RESEARCH_ARTICLES_README.md research_session_*.json
# git commit -m "Automated research update - $(date +%Y-%m-%d)"
# git push origin analog-discoveries-ip-protected

echo "=========================================="
echo "Daily research complete!"
echo "=========================================="

