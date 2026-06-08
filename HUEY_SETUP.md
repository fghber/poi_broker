# Huey Task Queue Configuration

## Overview

This application uses **Huey** as its task queue for asynchronous job processing, specifically for the bulk data export feature. The setup supports both development and production modes:

- **Development** (default): Synchronous execution using `MemoryHuey`
- **Production**: Asynchronous execution using `SqliteHuey` with background worker

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Flask Application                        │
│  - Handles HTTP requests                                         │
│  - Enqueues tasks to Huey                                        │
│  - Routes: POST /export, GET /export/download/<task_id>          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ Task Queue
                              │ (SQLite DB)
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Huey Consumer Worker                        │
│  - Polls task queue                                              │
│  - Executes background tasks                                     │
│  - Updates task status in database                               │
│  - Writes exported CSV files                                     │
└─────────────────────────────────────────────────────────────────┘
```

## Environment Configuration

### Development Mode (Default)

```bash
export SECRET_KEY=your-secret-key
# HUEY_BACKEND=memory is the default, no need to set explicitly
# Tasks execute synchronously (immediately when called)
```

**Pros:**
- No external dependencies (Redis, etc.)
- Easier debugging
- Suitable for local development

**Cons:**
- Blocks the web request until task completes
- Not suitable for long-running tasks
- No task persistence

### Production Mode

```bash
export SECRET_KEY=your-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
export ALERTS_DB_PATH=/var/lib/poi_broker/ztf_alerts_stream.db
export USERS_DB_PATH=/var/lib/poi_broker/users.db
export HUEY_THREADS=4
export HUEY_LOGFILE=/var/log/poi_broker/huey.log
```

**Pros:**
- Non-blocking: web app returns immediately
- Persistent queue: tasks survive restarts
- Scalable: multiple workers can process tasks
- Suitable for long-running operations

**Cons:**
- Requires background worker process
- Slightly more complex setup
- Database contention with Flask app

## Environment Variables

| Variable | Default | Values | Purpose |
|----------|---------|--------|---------|
| `HUEY_BACKEND` | `memory` | `memory`, `sqlite` | Task queue backend |
| `HUEY_SQLITE_PATH` | `instance/huey.db` | file path | SQLite database location |
| `HUEY_IMMEDIATE` | `true` | `true`, `false` | Synchronous vs async execution |
| `HUEY_THREADS` | `4` | integer | Number of worker threads |
| `HUEY_LOGFILE` | `huey.log` | file path | Worker process logfile |

## Installation & Setup

### 1. Ensure Huey is Installed

```bash
pip install -r requirements.txt
# huey[sqlite]==3.* should be installed
```

### 2. Development: Run with Synchronous Execution

```bash
export SECRET_KEY=test-secret-key
python -m flask --app wsgi:app run --debug
```

Tasks execute immediately when called (no separate worker needed).

### 3. Production: Run with Background Worker

**Terminal 1 - Start the Flask App:**

```bash
export SECRET_KEY=your-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
export ALERTS_DB_PATH=/var/lib/poi_broker/ztf_alerts_stream.db
export USERS_DB_PATH=/var/lib/poi_broker/users.db

# Using Gunicorn (recommended)
gunicorn -w 4 -b 0.0.0.0:8000 wsgi:app

# Or Flask development server
python -m flask --app wsgi:app run
```

**Terminal 2 - Start the Huey Worker:**

```bash
# Set the same environment variables as the Flask app
export SECRET_KEY=your-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
export ALERTS_DB_PATH=/var/lib/poi_broker/ztf_alerts_stream.db
export USERS_DB_PATH=/var/lib/poi_broker/users.db

# On Linux/Mac
./run_huey_worker.sh

# On Windows
run_huey_worker.bat

# Or run directly
huey_consumer poi_broker.extensions.huey --workers=4 --worker-type=thread --verbose
```

## Testing

### Run All Tests

```bash
# Development mode (default)
python -m pytest tests/test_huey_tasks.py -v

# Production mode with SQLite backend
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./test_huey.db
python -m pytest tests/test_huey_tasks.py -v
```

### Manual Testing Workflow

1. **Start the worker** (in separate terminal):
   ```bash
   export HUEY_BACKEND=sqlite
   export SECRET_KEY=test-secret-key
   ./run_huey_worker.sh
   ```

2. **Start the Flask app** (in another terminal):
   ```bash
   export HUEY_BACKEND=sqlite
   export SECRET_KEY=test-secret-key
   python -m flask --app wsgi:app run --debug
   ```

3. **Open the application**:
   - Navigate to http://localhost:5000/
   - Login to your account
   - Go to `/export` page

4. **Submit an export query**:
   - Build a query using the visual query builder
   - Click "Start Export"
   - Watch the status update in real-time

5. **Monitor the worker**:
   ```
   [2026-06-08 10:30:45] EXECUTING: poi_broker.tasks.create_export_file
   [2026-06-08 10:30:47] Task: create_export_file returned: None
   ```

6. **Verify task completion**:
   - Status should change from PENDING → RUNNING → SUCCESS
   - Download link should appear on the /export page

### Inspect Task Queue

```bash
# View queued tasks (development mode, in-memory)
python -c "from poi_broker.extensions import huey; print(huey.all_tasks())"

# View task database (production mode with SQLite)
sqlite3 ./instance/huey.db
sqlite> SELECT * FROM task;
sqlite> SELECT * FROM taskresult;
```

## Performance Tuning

### For High-Volume Exports

Adjust worker threads and database settings:

```bash
export HUEY_THREADS=8              # More worker threads
export HUEY_SQLITE_PATH=./huey.db  # Fast local storage
export HUEY_LOGFILE=/dev/null      # Disable logging for speed
```

### SQLite WAL Mode

The Huey SQLite backend is configured with:
- `journal_mode='wal'`: Write-Ahead Logging for better concurrency
- `timeout=5`: 5-second timeout for lock waits
- `cache_mb=8`: 8MB cache for performance
- `fsync=False`: Relax fsync for speed (at cost of durability)

Adjust these in `poi_broker/huey_config.py` if needed.

## Troubleshooting

### Tasks Not Executing

**Problem:** Export task stays in PENDING forever

**Solutions:**
1. Verify worker is running: `ps aux | grep huey_consumer`
2. Check worker logfile: `tail -f huey.log`
3. Verify `HUEY_BACKEND=sqlite` is set (not just in app, but also in worker)
4. Restart both app and worker with same environment variables

### SQLite Database Locked

**Problem:** "database is locked" errors

**Solutions:**
1. Ensure only one Huey worker is running
2. Remove stale lock files: `rm instance/huey.db-shm instance/huey.db-wal`
3. Stop all processes, wait 5 seconds, restart
4. Check for stuck Python processes: `ps aux | grep python`

### Tasks Running Synchronously in Production

**Problem:** Exports still block in production

**Cause:** `HUEY_BACKEND` not set correctly or `HUEY_IMMEDIATE=true`

**Solution:**
```bash
# Verify settings
echo $HUEY_BACKEND  # Should be 'sqlite'
echo $HUEY_IMMEDIATE  # Should be 'false' or unset

# Restart with correct settings
export HUEY_BACKEND=sqlite
./run_huey_worker.sh
```

### Worker Process Crashes

**Problem:** Huey consumer exits with traceback

**Solutions:**
1. Check SECRET_KEY is set: `echo $SECRET_KEY`
2. Check database paths exist: `ls -la $ALERTS_DB_PATH $USERS_DB_PATH`
3. View full error in logfile: `cat huey.log`
4. Run with verbose logging: `huey_consumer ... --verbose --logfile=-` (output to console)

## Advanced Configuration

### Use a Custom Backend

To use a different backend (e.g., Redis), modify `poi_broker/huey_config.py`:

```python
if backend_type == 'redis':
    from huey import RedisHuey
    return RedisHuey(
        'poi_broker',
        host=os.environ.get('REDIS_HOST', 'localhost'),
        port=int(os.environ.get('REDIS_PORT', 6379)),
    )
```

### Multiple Worker Processes

For better concurrency, run multiple workers:

```bash
# Terminal 1
huey_consumer poi_broker.extensions.huey --workers=4 --worker-type=thread

# Terminal 2
huey_consumer poi_broker.extensions.huey --workers=4 --worker-type=thread
```

### Process Manager (systemd, supervisord)

Example systemd service:

```ini
# /etc/systemd/system/poi_broker_huey.service
[Unit]
Description=POI Broker Huey Worker
After=network.target

[Service]
Type=simple
User=poi_broker
WorkingDirectory=/opt/poi_broker
Environment="HUEY_BACKEND=sqlite"
Environment="HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db"
Environment="SECRET_KEY=..."
ExecStart=/opt/poi_broker/.venv/bin/huey_consumer poi_broker.extensions.huey --workers=4 --logfile=/var/log/poi_broker/huey.log
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

Enable and start:
```bash
sudo systemctl enable poi_broker_huey
sudo systemctl start poi_broker_huey
sudo systemctl status poi_broker_huey
```

## Files

- `poi_broker/huey_config.py` - Huey backend configuration
- `poi_broker/extensions.py` - Extension initialization
- `poi_broker/tasks.py` - Background task definitions
- `run_huey_worker.sh` - Startup script (Linux/Mac)
- `run_huey_worker.bat` - Startup script (Windows)
- `tests/test_huey_tasks.py` - Test suite

## References

- [Huey Documentation](https://huey.readthedocs.io/)
- [SQLAlchemy with Flask](https://flask-sqlalchemy.palletsprojects.com/)
- [Gunicorn WSGI Server](https://gunicorn.org/)
