# POI Broker — Adding the Async Export Feature to an Existing Production Site

This guide is for a **Linux production site that is already fully configured**
(Ubuntu 24.04+, Python + venv + dependencies + Nginx + Gunicorn already running
the POI Broker web app). **Only the Huey/export feature is new.**

It covers:
1. What the new feature adds to the running system.
2. The minimal steps to deploy it (new env vars, new systemd unit, no Nginx change).
3. How to test the feature "like production" before/after deploying.

If you are setting up the site from scratch instead, see the full guide in
`HUEY_SETUP.md` / `HUEY_QUICK_START.md`.

---

## What the feature adds

```
Browser
   │  HTTPS
   ▼
Nginx (already running)  ── proxy ──▶ Gunicorn (already running, 4 workers)
                                          │  Flask app (wsgi:app)
                                          │
                                          │  NEW: enqueues export jobs
                                          ▼
                                  Huey SQLite queue
                                  /var/lib/poi_broker/huey.db   (NEW)
                                          ▲
                                          │  polls & executes
                                  Huey consumer (NEW systemd: poi-broker-huey)
                                          │
                                          ▼
                                  writes CSV to /var/lib/poi_broker/exports/
```

- **Gunicorn (unchanged):** the 4 existing workers now also import the shared
  `poi_broker.extensions.huey` instance, which resolves to the SQLite queue
  (`HUEY_SQLITE_PATH`). All 4 workers coordinate on **one** task queue.
- **Huey consumer (NEW):** a separate long-running process (systemd) that polls
  the queue and executes `create_export_file`. It must be supervised so a crash
  cannot strand exports.
- **ExportTask (NEW):** rows in `users.db` track user-facing status
  (PENDING → RUNNING → SUCCESS/FAILED). The stale-task guard resets stuck tasks
  after 30 minutes.

**Nginx does not need any change** — the export feature is served by the same
Gunicorn app and the same `/export` routes. The existing `X-Forwarded-*` headers
are already in place.

---

## 1. Update the code & dependencies

Pull the latest code into the existing deployment and install the new
dependency (`huey`):

```bash
cd /opt/poi_broker
sudo -u poi-broker git pull
sudo -u poi-broker .venv/bin/pip install -r requirements.txt   # adds huey
```

> The alerts DB and users DB are unchanged and stay where they already are
> (`ALERTS_DB_PATH` / `USERS_DB_PATH`). The new `export_task` table is created
> automatically in `users.db` by SQLAlchemy on first app start.

## 2. Add the Huey environment variables

The web app needs two new env vars. Add them to the **existing** env file that
Gunicorn already reads (e.g. `/etc/poi-broker/web.env`):

```bash
HUEY_BACKEND=sqlite
HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
```

Create a **new** env file for the worker, `/etc/poi-broker/huey.env`. It must
set the **same** `HUEY_SQLITE_PATH` and DB paths as the web app so both
processes share one queue:

```bash
SECRET_KEY="<same as web.env>"
ALERTS_DB_PATH=/var/lib/poi_broker/ztf_alerts_stream.db
USERS_DB_PATH=/var/lib/poi_broker/users.db
HUEY_BACKEND=sqlite
HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
HUEY_THREADS=2
HUEY_LOGFILE=/var/log/poi_broker/huey.log
```

> **Important:** `HUEY_SQLITE_PATH` must be an **absolute path** and identical in
> both files. If the worker cannot resolve it (no Flask app context), it raises a
> clear error rather than silently creating a second empty queue.

Set ownership and permissions:

```bash
sudo chown root:poi-broker /etc/poi-broker/huey.env
sudo chmod 640 /etc/poi-broker/huey.env
```

## 3. Add the worker systemd unit (NEW)

Only the **worker** unit is new. The web app keeps its existing systemd unit —
just restart it after the env change.

```bash
sudo cp deploy/poi-broker-huey.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now poi-broker-huey
```

Restart the web app so it picks up the new env vars:

```bash
sudo systemctl restart poi-broker-web
```

Verify both:

```bash
sudo systemctl status poi-broker-web
sudo systemctl status poi-broker-huey
tail -f /var/log/poi_broker/huey.log
```

> The worker unit points the consumer at **`poi_broker.worker.huey`** (not
> `poi_broker.extensions.huey`). This is essential: `worker.py` imports
> `poi_broker.tasks`, which registers `create_export_file` on the shared Huey
> instance. A worker started against `extensions.huey` runs but never executes
> any task.

### Worker sizing (how many threads?)

The Huey consumer runs **one process** with `HUEY_THREADS` worker threads. This
is **independent** of the Gunicorn web worker count — they don't need to match.

- **`gunicorn -w 4`** = 4 web processes that *enqueue* export jobs.
- **`huey_consumer --workers=N`** = N threads in one process that *execute* jobs.

For a modest deployment, **1 consumer with 2 threads** is the recommended
default (already set in the run scripts and deploy units). This gives:

- 2 concurrent export slots (plus huey's independent periodic scheduler for
  cleanup/prune tasks).
- Minimal resource footprint.
- No SQLite lock contention from multiple consumer processes.

Only increase `HUEY_THREADS` (e.g. to 3–4) if you observe exports queueing up
under concurrent users, or if a single export is very long. You generally do
**not** need multiple `huey_consumer` processes unless throughput demands it.

## 4. Nginx — no change required

The existing Nginx config already proxies `/export` to Gunicorn and already sets
the `X-Forwarded-*` headers the app needs. **Do not change Nginx.**

---

## 5. Verify the full flow

1. Open the site and log in.
2. Navigate to **/export**, build a query, click **Start Export**.
3. The page should return immediately ("Export started").
4. Watch the worker process it:

```bash
tail -f /var/log/poi_broker/huey.log
# [..] EXECUTING: poi_broker.tasks.create_export_file
# [..] ExportTask <id> completed successfully
```

5. The export page should show **Export Completed** with a download link.

### Sanity checks

```bash
# Worker has the task registered (should print a task object, not None)
sudo -u poi-broker /opt/poi_broker/.venv/bin/python -c \
  "from poi_broker.worker import huey; print(huey._registry.string_to_task('poi_broker.tasks.create_export_file'))"

# Queue DB exists and is shared
ls -la /var/lib/poi_broker/huey.db
```

---

## Troubleshooting

### Export stays PENDING forever
1. Is the worker running? `sudo systemctl status poi-broker-huey`
2. Check the worker log: `tail -f /var/log/poi_broker/huey.log`
3. Is `HUEY_BACKEND=sqlite` set in **both** env files?
4. Is `HUEY_SQLITE_PATH` the **same absolute path** in both files?
5. Does the worker target `poi_broker.worker.huey`? (see unit file)

### "database is locked"
- Ensure only **one** Huey consumer is running: `ps aux | grep huey_consumer`
- Remove stale lock files: `rm /var/lib/poi_broker/huey.db-shm /var/lib/poi_broker/huey.db-wal`
- Restart the worker: `sudo systemctl restart poi-broker-huey`

### Exports run synchronously (block the request)
- `HUEY_BACKEND` is not `sqlite`, or `HUEY_IMMEDIATE=true` is set. Verify both
  env files.

### Stale exports blocking a user
- The stale-task guard resets `PENDING`/`RUNNING` exports older than 30 minutes
  to `FAILED` (at app startup and every 5 min via the periodic task). This frees
  the user's active-export slot so they can retry.

---

## Files in `deploy/`

| File | Purpose |
|------|---------|
| `poi-broker-huey.service` | **NEW** systemd unit for the Huey consumer |
| `supervisord-poi-broker-huey.conf` | supervisord alternative for the worker |
| `poi-broker-web.service` | reference only — your existing web unit is unchanged |
| `nginx-poi-broker.conf` | reference only — your existing Nginx config is unchanged |
