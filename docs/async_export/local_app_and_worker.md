# Test async export locally (app + Huey worker)

Run the Flask app and Huey worker as separate processes with the SQLite
backend — the same code path as production, without Nginx or systemd.

Both terminals need the same environment:

```bash
export SECRET_KEY=test-secret-key
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export ALERTS_DB_PATH=../_broker_db/ztf_alerts_stream.db
export USERS_DB_PATH=../_broker_db/users.db
```

**Terminal 1 — app**

```bash
python -m flask --app wsgi:app run
```

**Terminal 2 — worker** (target `poi_broker.worker.huey`, not `extensions.huey`)

```bash
python -m huey.bin.huey_consumer poi_broker.worker.huey --workers=2 --worker-type=thread --verbose
```

Log in, go to `/export`, build a query, click **Start Export**. The page should return immediately. The worker should log:

```
EXECUTING: poi_broker.tasks.create_export_file
```

The page then shows **Export Completed** with a download link.

Confirm the shared queue DB exists and has tables:

```bash
ls -la ./instance/huey.db
sqlite3 ./instance/huey.db .tables
```

Production (existing Ubuntu site): `docs/async_export/ubuntu_add_worker.md`.
Decision record (Huey + SQLite vs Celery + Redis): `docs/async_export/adr_huey_sqlite.md`.
