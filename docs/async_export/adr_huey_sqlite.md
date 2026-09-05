# ADR: Huey + SQLite for async CSV export

- **Status:** Accepted
- **Date:** 2026-06
- **Feature:** Visual-query-builder CSV export (`/export`)

## Context

Exporting alert rows matching a query-builder filter can take long enough to
block a Gunicorn worker. The site already runs on Ubuntu with Nginx + Gunicorn
and stores data in SQLite (`ALERTS_DB_PATH`, `USERS_DB_PATH`). The team is
small; adding infrastructure that needs its own server process is costly.

The queue had to:

- Run export jobs off the HTTP request path.
- Persist jobs across process restarts.
- Fit the existing SQLite operations model.
- Avoid a new networked service (Redis, RabbitMQ, etc.).

## Decision

Use **Huey** with the **SQLite** backend (`SqliteHuey`) in production, and
**MemoryHuey** (immediate/synchronous) for local development and pytest.

Production adds **one** extra process: the Huey consumer, supervised by
**systemd** (`Restart=always`). See `docs/async_export/ubuntu_add_worker.md`.

## Alternatives considered

### Celery + Redis (rejected)

Celery is the usual Flask choice. It would require Redis (or RabbitMQ) as a
broker: another package, another daemon, another failure domain, and ongoing
maintenance. That extra server process was the deal-breaker. Celery’s footprint
is also larger than needed for a single job type (CSV export) at modest
concurrency.

### Huey + Redis (rejected)

Huey can use Redis. That still adds a Redis process. The only gain over
Huey + SQLite is higher throughput, which this feature does not need.

### In-request / thread pool only (rejected)

A thread or `multiprocessing` pool inside Gunicorn does not survive worker
restarts, does not give a durable job record, and still ties export CPU to the
web process.

## Consequences

- **Smaller footprint:** `huey[sqlite]` in `requirements.txt`; queue file is
  another SQLite DB (`HUEY_SQLITE_PATH`). No Redis.
- **One new process, not a new service:** the consumer is a Python process
  started with `python -m huey.bin.huey_consumer poi_broker.worker.huey`.
  systemd is enough; do not add supervisord, wrapper `.sh`/`.bat` scripts, or
  a second consumer on the same queue.
- **SQLite limits:** one consumer (default 2 threads). WAL is enabled in
  `poi_broker/huey_config.py`. App and worker must share the **same absolute**
  `HUEY_SQLITE_PATH` or they silently use two empty queues.
- **Dev vs prod:** default `HUEY_BACKEND=memory` runs tasks in-process (no
  worker). In that mode `POST /export` blocks until the CSV is written; keep
  local queries small or lower `EXPORT_MAX_ROWS`. There is no consumer, so
  periodic stale-cleanup never runs: unstick a Ctrl+C with
  `python -m flask --app wsgi:app fail-stale-exports --force` (same env as
  the app). Production sets `HUEY_BACKEND=sqlite`. `HUEY_IMMEDIATE` applies
  only to the memory backend.
- **User-facing state** lives in `ExportTask` (`users.db`), not in Huey’s
  result table. The consumer’s periodic stale-task guard
  (`cleanup_stale_export_tasks`, every 5 minutes) fails all stale
  `PENDING`/`RUNNING` rows after 30 minutes so a dead worker cannot block a
  user’s export slot forever. Independently, `POST /export` may fail **the
  current user’s** stale active row (same age window) before returning 409, so
  a down consumer cannot strand a user forever. Web processes do not run
  global cleanup at startup. The worker only transitions from expected
  statuses (compare-and-swap) and heartbeats `updated_at` while `RUNNING`, so
  a live long export is not false-failed and cannot be resurrected after a
  stale `FAILED`. Keep `EXPORT_STALE_MAX_AGE_SECONDS >= 2 *
  EXPORT_HEARTBEAT_SECONDS` (defaults 1800 vs 60) so a healthy export is not
  false-failed between heartbeats.
- **One Flask app per consumer process:** Huey handlers reuse a lazy
  process-wide app (`_get_worker_app` in `tasks.py`) instead of calling
  `create_app()` on every task.

## Invariants (for implementers and coding agents)

These are easy to get wrong; the consumer will look "fine" while doing nothing.

1. Point the consumer at **`poi_broker.worker.huey`**, never
   `poi_broker.extensions.huey`. `worker.py` imports `poi_broker.tasks` so
   `@huey.task` handlers register. `extensions.huey` only creates the instance.
2. Start the consumer with **`python -m huey.bin.huey_consumer`**. The name
   `huey_consumer` is a console script, not an importable module.
3. Do not introduce Celery, Redis, or a second queue backend without a new ADR.
4. Do not run export work in the request thread when `HUEY_BACKEND=sqlite`.
5. Do not call `reset_stale_export_tasks()` from `create_app()` / `init_huey()`.
   The Huey consumer’s periodic `cleanup_stale_export_tasks` is the **global**
   owner. `POST /export` may call it with `user_id=` (age-aware only, never
   `ignore_age`) to free that user’s stale slot. Leave periodic scheduling
   enabled (do not pass `--no-periodic`).
6. Automated tests: `tests/test_huey_tasks.py` (memory backend via `conftest.py`).
   Manual two-process check: `docs/async_export/local_app_and_worker.md`.

## Related code

| Path | Role |
|------|------|
| `poi_broker/huey_config.py` | Memory vs SQLite instance |
| `poi_broker/extensions.py` | Shared `huey` object |
| `poi_broker/worker.py` | Consumer entrypoint (registers tasks) |
| `poi_broker/tasks.py` | `create_export_file`, stale-task / retention jobs |
| `poi_broker/cli.py` | `fail-stale-exports` one-off (local `--force` unstick) |
| `poi_broker/routes/export.py` | Enqueue + download; age-aware per-user stale fail |
| `docs/async_export/poi-broker-huey.service` | systemd unit template |
