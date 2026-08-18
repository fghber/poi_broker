I've completed a thorough inspection of the Huey async export infrastructure. Here's my analysis:

## Summary of Findings

### ✅ Correct Aspects

1. **Initialization order fix is correct** — `init_huey(app)` now runs after `db.init_app(app)` and the SQLite pragma configuration, so `reset_stale_export_tasks()` has a bound DB when it queries `ExportTask`.

2. **Compare-and-swap status transitions** — `_transition_export_task()` uses `status IN from_statuses` guards, preventing a stale-guard `FAILED` from being overwritten back to `RUNNING`/`SUCCESS` by a still-running worker.

3. **Heartbeat mechanism** — The worker refreshes `updated_at` every `EXPORT_HEARTBEAT_SECONDS` (default 60s) while `RUNNING`, so long exports aren't false-failed by the 30-minute stale guard.

4. **Partial file cleanup** — Failed exports unlink any partially-written CSV, preventing orphan files.

5. **Worker app singleton** — `_get_worker_app()` uses a process-wide lazy singleton with a lock, avoiding `create_app()` per task.

6. **Consumer entrypoint** — `worker.py` correctly imports `tasks` to register `@huey.task` handlers on the shared Huey instance.

### ⚠️ Issues Found

1. **Enqueue path ambiguity** — In `routes/export.py`, `create_export_file(...)` is called directly (not `.enqueue()`). In `immediate=True` (memory) mode this runs synchronously in the request thread. In `immediate=False` (sqlite) mode, Huey's `Task.__call__` enqueues the task. This works but is subtle — the route comment says "Enqueue the background task" but in dev mode it actually blocks. This is documented behavior but worth noting.

2. **Periodic tasks only run in the consumer** — `cleanup_stale_export_tasks` and `prune_huey_results` are `@huey.periodic_task` decorators. These only execute when the Huey consumer is running with periodic scheduling enabled. The startup `reset_stale_export_tasks()` in `init_huey()` covers the web app side, but `cleanup_expired_exports()` (file + DB row cleanup) only runs in the consumer.

3. **SQLite concurrency** — Both the web app and Huey consumer write to `users.db` (for `ExportTask` rows). WAL mode is enabled, but with multiple Gunicorn workers + the Huey consumer, there's potential for write contention on the `uix_export_task_one_active_per_user` partial unique index. The `IntegrityError` handling in `export_submit` catches the race for concurrent export creation, but the compare-and-swap transitions in the worker don't have retry logic for `sqlite3.OperationalError` (database locked).

4. **No retry on DB lock** — `_transition_export_task()` catches exceptions and rolls back, but doesn't retry on `sqlite3.OperationalError` (database is locked). Under high concurrency, a heartbeat or status transition could fail transiently.

5. **Test coverage gap** — `test_export_routes_exist` checks for `/export/status/<int:task_id>` but the route is registered at `/export/status/<int:task_id>` (not `/api/export/status/<int:task_id>` as the docstring says). The test passes because it checks for `'export'` substring, but the docstring is misleading.

### 🔒 Security Assessment

- **No path traversal** — `file_path` is never exposed to clients; downloads go through `/export/download/<id>` with ownership checks.
- **Ownership checks** — Both `export_download` and `export_status` verify `export_task.user_id == current_user.id`.
- **CSRF protection** — The `POST /export` route is behind `@login_required` and Flask-WTF CSRF (JSON 400 for `/api/*` paths, but `/export` is not under `/api/` prefix — it uses form-based CSRF).
- **Rate limiting** — `POST /export` has `READ_RATE_LIMIT_LAX` (30/min).

### 🐳 Parallel Execution Safety (Nginx + Gunicorn + Huey)

- **Multiple Gunicorn workers** can each call `create_export_file()` — Huey's SQLite storage handles concurrent enqueues with its own locking.
- **The `uix_export_task_one_active_per_user` partial unique index** prevents two workers from creating concurrent active tasks for the same user (the `IntegrityError` handler returns 409).
- **Huey consumer with `--workers=2 --worker-type=thread`** — two threads share one Huey instance and one `_worker_app`. The `_worker_app_lock` ensures only one `create_app()` call.
- **Potential issue**: If two Gunicorn workers both try to enqueue for the same user simultaneously, one will hit `IntegrityError` on the unique index and return 409. This is correct behavior.

### Recommendations

1. **Add retry logic** for `sqlite3.OperationalError` in `_transition_export_task()` — transient lock errors during heartbeat/status transitions should retry, not silently fail.

2. **Clarify the enqueue path** — Consider using `create_export_file.enqueue()` explicitly for clarity, or document that direct call is intentional (Huey handles both modes).

3. **Fix the `export_status` docstring** — It says `GET /api/export/status/<task_id>` but the actual route is `/export/status/<int:task_id>`.

4. **Consider `HUEY_THREADS` tuning** — The systemd unit uses `--workers=${HUEY_THREADS}` (default 2). For large exports, more threads could help, but SQLite write locking limits gains.

The initialization-order fix is solid and the overall architecture is sound for the stated use case. The main risk is SQLite write contention under high concurrency, which is inherent to the chosen stack (Huey + SQLite, per the ADR).