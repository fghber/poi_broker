

The fix is verified correct and the codebase is in good shape. Here's my final verdict on the four angles you asked about.

## Correctness ✅

The init-order fix is sound:
- `db.init_app(app)` + the SQLite pragma/app-context block now run **before** `init_huey(app)`, so `reset_stale_export_tasks()` queries `ExportTask` against a bound engine.
- The worker path is unaffected: `_get_worker_app()` builds its own app via `create_app()` (which now has the correct order internally), so the consumer's startup reset also works.
- Compare-and-swap transitions (`_transition_export_task`) are correctly scoped — a stale-guard `FAILED` cannot be resurrected, and a live `RUNNING` export heartbeats `updated_at` so it isn't false-failed. The `uix_export_task_one_active_per_user` partial unique index backs the "one active export per user" invariant at the DB level.

## Security ✅

No regressions, and the export path is clean:
- `export_status` deliberately omits `file_path` (server filesystem path); downloads go through a separate route that re-checks `export_task.user_id == current_user.id`.
- `export_submit` validates the query payload, enforces `MAX_EXPORT_ROWS` via a pre-count, and handles `IntegrityError` → 409.
- Query-builder input is validated (operator whitelist, field format, LIKE escaping, indexed-ID operator restrictions) so user rules can't widen patterns or hit unindexed columns.
- Session cookies are `HttpOnly`/`Secure`/`SameSite`, CSRF is active, rate limits apply, and email links use `PUBLIC_BASE_URL` (falling back to `X-Forwarded-Host` only with `x_host=1`).

## Gaps ⚠️ (minor, pre-existing)

1. **`logging.basicConfig(FileHandler("app.log"))`** — relative path resolved against the worker's cwd. Under Gunicorn without `preload_app`, each worker calls `create_app()`; `basicConfig` is a no-op once handlers exist, so only the first worker writes `app.log`. Not a crash, but logging is silently incomplete in multi-worker setups. Consider an absolute, configurable log path.
2. **Web-app startup reset runs in every Gunicorn worker** — `reset_stale_export_tasks()` runs once per worker at boot. It's idempotent and wrapped in try/except, but 4 concurrent writers on `users.db` at boot can transiently hit `database is locked` (harmless, logged, returns 0). The worker's 5-minute periodic task is the real guard; the web-side call is just a fast-path.

## Parallel safety ✅

- **Nginx ↔ Gunicorn:** stateless app; `ProxyFix(x_for=1, x_proto=1, x_host=1)` trusts exactly one hop, and `PUBLIC_BASE_URL` prevents Host-poisoning of emailed links.
- **Gunicorn workers:** share `users.db`/`alerts.db` via SQLite WAL (`journal_mode=WAL`, `synchronous=NORMAL`). Reads don't block writes. The only cross-worker write contention is the one-time startup reset and request-time session commits — both short.
- **Huey consumer:** single process, `HUEY_THREADS` threads (default 2), one SQLite queue. The ADR's "one consumer, no second queue" rule is enforced by `get_huey_sqlite_path()` raising (not silently falling back to cwd) when no app context exists. Task handlers reuse a process-wide app (`_get_worker_app`, double-checked locked singleton — covered by `test_get_worker_app_is_process_singleton_under_concurrency`).
- **Web ↔ worker coordination:** user-facing state lives in `ExportTask`, not Huey results. Enqueue (web) and execution (worker) are decoupled; the stale-guard + compare-and-swap transitions keep the two processes from clobbering each other.

**Bottom line:** the fix is correct and safe. The two gaps above are cosmetic/observability items worth a follow-up, not blockers. Full suite: **501 passed**, `flask --app wsgi:app routes` boots cleanly.