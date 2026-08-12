How to test "like production"
Locally (closest to production, no server needed): run the app and worker as separate processes with the SQLite backend, exactly as production does:

```bash
# Terminal 1 — app (sqlite backend, async)
export SECRET_KEY=test-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export ALERTS_DB_PATH=../_broker_db/ztf_alerts_stream.db
export USERS_DB_PATH=../_broker_db/users.db
python -m flask --app wsgi:app run

# Terminal 2 — worker (same env)
export SECRET_KEY=test-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export ALERTS_DB_PATH=../_broker_db/ztf_alerts_stream.db
export USERS_DB_PATH=../_broker_db/users.db
python -m huey.bin.huey_consumer poi_broker.worker.huey --workers=2 --worker-type=thread --verbose
```

Then log in, go to /export, build a query, click Start Export. Watch Terminal 2 process it and the page flip to "Export Completed" with a download link. This exercises the exact same code path as production (SQLite queue + separate worker), just without Nginx/systemd.