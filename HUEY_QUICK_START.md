# Huey Production Setup - Quick Start Guide

## Overview

This guide shows how to test the complete Huey production setup with SQLite backend and background worker.

## Quick Start: Development Mode (Default)

The easiest way to get started is development mode (synchronous execution):

```bash
export SECRET_KEY=test-secret-key
python -m flask --app wsgi:app run --debug
# Navigate to http://localhost:5000/export
# Export tasks execute immediately (no worker needed)
```

## Production Setup: Async Execution with Worker

For production-like testing with background task processing:

### Step 1: Terminal 1 - Start Flask App

```bash
export SECRET_KEY=test-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
python -m flask --app wsgi:app run
```

Output should show:
```
Huey task queue initialized with Memory backend (immediate/synchronous)
# This is normal - each terminal inherits different environment
```

### Step 2: Terminal 2 - Start Huey Worker

```bash
export SECRET_KEY=test-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export HUEY_THREADS=4

# On macOS/Linux:
./run_huey_worker.sh

# On Windows:
python -m huey_consumer poi_broker.extensions.huey --workers=4 --worker-type=thread --verbose
```

Expected output:
```
[2026-06-08 14:30:45] Huey v3.0.1 | Worker: poi_broker | PID: 12345
[2026-06-08 14:30:45] Worker started with 4 threads
[2026-06-08 14:30:45] Consuming tasks from...
```

### Step 3: Terminal 3 - Test the Export Feature

```bash
# Open browser to http://localhost:5000/
# Login to your account
# Navigate to /export
# Build a query and click "Start Export"
```

### Step 4: Monitor Progress

Watch Terminal 2 (Huey worker):
```
[2026-06-08 14:31:02] EXECUTING: poi_broker.tasks.create_export_file
  task_id: 550e8400-e29b-41d4-a716-446655440000
  user_id: 1
[2026-06-08 14:31:03] Task returned: None
[2026-06-08 14:31:03] Status updated to: SUCCESS
```

## Testing the Setup

### Run Huey-Specific Tests

```bash
# Development mode tests
python -m pytest tests/test_huey_tasks.py::test_huey_imports -v
python -m pytest tests/test_huey_tasks.py::test_app_with_memory_huey -v

# Production mode tests
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./test_huey.db
python -m pytest tests/test_huey_tasks.py::test_app_with_sqlite_huey -v
```

### Verify Task Queue

**Check if Huey database exists:**
```bash
ls -la ./instance/huey.db
```

**Inspect task history (with SQLite):**
```bash
sqlite3 ./instance/huey.db '.tables'
sqlite3 ./instance/huey.db 'SELECT id, status, created_at FROM task LIMIT 5;'
```

## Environment Variables Reference

| Variable | Example | Purpose |
|----------|---------|---------|
| `HUEY_BACKEND` | `sqlite` | Use `memory` for dev, `sqlite` for production |
| `HUEY_SQLITE_PATH` | `./instance/huey.db` | Location of SQLite database |
| `HUEY_IMMEDIATE` | `false` | Set to `false` to run tasks async |
| `HUEY_THREADS` | `4` | Number of worker threads |
| `HUEY_LOGFILE` | `huey.log` | Worker logfile path |

## Troubleshooting

### Issue: "huey_consumer: command not found"

**Solution:** Use Python module instead:
```bash
python -m huey_consumer poi_broker.extensions.huey --workers=4 --worker-type=thread --verbose
```

### Issue: Tasks run synchronously in production

**Cause:** `HUEY_BACKEND=sqlite` not set in both app and worker terminals

**Solution:**
```bash
# Verify setting
echo $HUEY_BACKEND  # Should output: sqlite

# Set it if missing
export HUEY_BACKEND=sqlite
```

### Issue: "database is locked" errors

**Solution:** Ensure only one Huey worker is running
```bash
# Find and kill any stuck processes
ps aux | grep huey_consumer
kill -9 <PID>

# Clean up lock files
rm ./instance/huey.db-shm ./instance/huey.db-wal
```

## Verifying the Complete Setup

Use this checklist to verify everything works:

- [ ] Flask app boots without errors: `python -m flask --app wsgi:app routes | grep export`
- [ ] Huey imports work: `python -c "from poi_broker.extensions import huey; print(f'Huey: {type(huey).__name__}')" `
- [ ] ExportTask model exists: `python -c "from poi_broker.models import ExportTask; print('OK')" `
- [ ] Development mode works: Start app with `HUEY_BACKEND=memory`, submit export, watch it complete immediately
- [ ] Production mode works: Start worker with `HUEY_BACKEND=sqlite`, start app same way, submit export, watch worker process it
- [ ] Task database created: `ls -la ./instance/huey.db` shows file

## Files Created/Modified

- **New Files:**
  - `poi_broker/huey_config.py` - Flexible Huey backend configuration
  - `run_huey_worker.sh` - Linux/macOS startup script
  - `run_huey_worker.bat` - Windows startup script
  - `tests/test_huey_tasks.py` - Comprehensive test suite
  - `HUEY_SETUP.md` - Detailed setup documentation

- **Modified Files:**
  - `poi_broker/extensions.py` - Uses new huey_config module
  - `poi_broker/__init__.py` - Updated Huey initialization

## Next Steps

After confirming the setup works:

1. **Production Deployment:** Update startup scripts for your server (systemd, Docker, etc.)
2. **Monitoring:** Add logs and metrics collection
3. **Scaling:** Run multiple Huey workers for higher throughput
4. **Persistence:** Consider backing up the SQLite task database

See `HUEY_SETUP.md` for detailed production configuration.
