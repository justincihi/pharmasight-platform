# Manual Autonomous Research Trigger Guide

## Quick Method (Recommended)

### Option 1: Via Dashboard UI
1. Navigate to **Scheduler Dashboard** in the admin panel
2. Click **"Run Now"** button (if available)
3. Monitor progress in real-time

### Option 2: Via Command Line

```bash
# Navigate to project directory
cd /home/ubuntu/pharmasight-admin-dashboard

# Run the trigger script
npx tsx trigger_research.ts
```

### Option 3: Via tRPC API Call

```bash
curl -X POST https://your-dashboard-url.com/api/trpc/scheduler.runNow \
  -H "Content-Type: application/json" \
  -H "Cookie: your-session-cookie" \
  -d '{}'
```

## Current Issue: Python SRE Module Mismatch

The autonomous research engine is failing due to a Python regex module conflict:

```
AssertionError: SRE module mismatch
```

### Root Cause
The Python environment has conflicting regex (SRE) module versions, preventing the research engine from running.

### Temporary Solution
Until the Python environment is fixed, you can:

1. **Import analogs manually** via the Platform API:
   ```bash
   python python_integration/pharmasight_dashboard_client.py
   ```

2. **Use the dashboard chatbot** to generate analog suggestions

3. **Manually add analogs** through the dashboard UI

### Permanent Fix Needed

The autonomous research engine (`server/runAutonomousResearch.ts`) needs to:
1. Fix Python environment dependencies
2. Update regex module to compatible version
3. Test with a simple discovery run

## How Manual Runs Work

When you trigger a manual run, the system:

1. **Executes Python research engine** (`src/analog_generation_fix.py` in pharmasight-platform repo)
2. **Discovers new analogs** based on configured parent compounds
3. **Imports discoveries** into the dashboard database
4. **Updates master_analogs.json** file
5. **Sends notifications** for high-confidence discoveries (>85%)

## Expected Output

A successful run should:
- Discover 5-50 new analogs (depending on research scope)
- Add them to the database with `discovered_by = 'autonomous_system'`
- Update the master analog list
- Send notification to dashboard owner

## Troubleshooting

### No new analogs after run
- Check Python error logs in server console
- Verify research engine configuration
- Ensure parent compounds are defined

### Import fails
- Check Platform API key is set (`PLATFORM_API_KEY`)
- Verify database connection
- Check master_analogs.json file permissions

### Scheduler not running automatically
- Verify `SCHEDULER_ENABLED` is not set to "false"
- Check cron schedule: `SCHEDULER_CRON="0 9 * * *"` (9 AM daily)
- Restart the server to reinitialize scheduler

## Future Improvements

1. Add "Run Now" button to Scheduler Dashboard UI
2. Show real-time progress during research runs
3. Display research run history and results
4. Add configurable research goals (target compounds, therapeutic areas)
5. Implement parallel research runs for multiple goals
