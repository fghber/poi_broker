All implementation details confirmed. Here is the comprehensive review.

---

# Huey Async CSV Export — Correctness Review

## Overall Assessment

The pipeline is **functionally correct** for the common case and follows a clean architecture (route → service → task → housekeeping). The dev/prod backend split, worker registration pattern, and periodic housekeeping are well-designed. However, there are several issues ranging from **low** to **high** severity that should be addressed before relying on this in production.

---

## 1. Race Condition: Active-Task Check (🔴 High)

**Location:** `routes/export.py` — `export_submit`

The flow is:
1. `SELECT` for any active (PENDING/RUNNING) task for this user
2. If none, `INSERT` a new ExportTask
3. Enqueue the huey task

Steps 1 and 2 are **not atomic**. Two concurrent `POST /export` requests from the same user can both pass the check and both create tasks. There is also **no unique constraint** on `ExportTask` to enforce one-active-task-per-user at the DB level.

**Recommendation:**
- Add a partial unique index (SQLite supports this via `CREATE UNIQUE INDEX ... WHERE status IN ('PENDING','RUNNING')`) or use an application-level advisory lock.
- Alternatively, wrap the check-and-create in a transaction with `SELECT ... FOR UPDATE` — though note SQLite locking semantics differ from Postgres. The simplest robust fix for SQLite is the partial unique index plus catching `IntegrityError` on insert.

---

## 2. Memory: `.all()` Loads Entire Result Set (🔴 High)

**Location:** tasks.py — `create_export_file`

```python
results = filtered_query.all()  # loads ALL rows into memory
```

With `EXPORT_MAX_ROWS` defaulting to 1,000,000, a worst-case export loads ~1M rows (each with ~30+ columns) into a Python list of dicts before writing a single CSV line. This risks **OOM kills** on the worker process.

**Recommendation:** Use SQLAlchemy streaming with `yield_per()` or server-side cursors:

```python
for row in filtered_query.yield_per(1000):
    writer.writerow(build_export_row(row[0], row[1]))
```

This keeps memory bounded regardless of result size. Note: SQLite's Python driver doesn't support true server-side cursors, but `yield_per` still batches ORM object creation and is far better than `.all()`.

---

## 3. Security: `file_path` Leaked in Status API (🟡 Medium)

**Location:** `routes/export.py` — `export_status`

```python
'file_path': export_task.file_path,
```

The status JSON response includes the **server filesystem path** to the CSV file. This is an information leak — the user learns the server's directory structure and instance path. The download endpoint already serves the file content; the path itself should never be exposed.

**Recommendation:** Remove `file_path` from the JSON response. If a path is needed for debugging, gate it behind an admin check or omit it entirely for regular users.

---

## 4. Enqueue Failure Leaves Orphaned PENDING Task (🟡 Medium)

**Location:** `routes/export.py` — `export_submit`

If `create_export_file(...)` raises (e.g., queue full, serialization error) **after** the ExportTask row is committed, the task row remains `PENDING` forever — until the 5-minute stale reset eventually fails it. The user sees a stuck "pending" status for up to 30 minutes (`EXPORT_STALE_MAX_AGE_SECONDS`).

**Recommendation:** Wrap the enqueue in try/except and delete or mark the ExportTask as FAILED on enqueue failure:

```python
try:
    create_export_file(query_params, user_id, task.id)
except Exception:
    export_task.status = 'FAILED'
    export_task.error_message = 'Failed to enqueue export task'
    db.session.commit()
    return jsonify({'error': 'Failed to start export'}), 500
```

---

## 5. No Partial File Cleanup on Write Failure (🟡 Medium)

**Location:** tasks.py — `create_export_file`

If the CSV write fails mid-stream (disk full, exception during iteration), the exception handler sets status to `FAILED` but does **not** delete the partial file left on disk. These accumulate as orphans until `cleanup_expired_exports` runs (which only cleans `SUCCESS` tasks).

**Recommendation:** In the exception handler, attempt to delete any partial file:

```python
except Exception as exc:
    if file_path and os.path.exists(file_path):
        try:
            os.remove(file_path)
        except OSError:
            logger.warning(f'Could not remove partial file {file_path}')
    # ... set FAILED status
```

---

## 6. Stale Task Reset Based on `updated_at` (🟢 Low — Correct, but worth noting)

**Location:** tasks.py — `reset_stale_export_tasks`

The cutoff is computed from `updated_at`, which is correct because:
- `PENDING` tasks get `updated_at` at creation
- `RUNNING` tasks get `updated_at` when the task sets status to RUNNING
- `onupdate` fires on every status transition

This is sound. The only edge case: if the worker crashes **between** committing `PENDING` and the task starting (so `updated_at` never advances from creation time), the task is correctly failed after the stale window. ✅

---

## 7. Startup Double-Reset (🟢 Low — Harmless)

**Location:** __init__.py — `init_huey`

`reset_stale_export_tasks()` runs at startup inside `init_huey()`. Since both the web app and the worker call `create_app()` (the worker imports tasks which import extensions, but the consumer process itself doesn't call `create_app`), this is actually **only** running in the web process. If you ever run multiple web workers (gunicorn with multiple workers), each will run the reset — but since it's idempotent (only fails tasks past the cutoff), this is harmless. ✅

---

## 8. Dev Immediate Mode Behavior (🟢 Low — Correct)

**Location:** huey_config.py

With `MemoryHuey` + `immediate=True`, `create_export_file()` executes **synchronously** within the `export_submit` request. The task sets `RUNNING`, then immediately overwrites with `SUCCESS` or `FAILED`. The user's subsequent status poll will see the final state. This is correct — no stuck `RUNNING` state in dev. ✅

One caveat: in immediate mode, the `POST /export` request blocks until the entire CSV is generated. For large exports this means a long HTTP request. Consider documenting this or setting a lower `EXPORT_MAX_ROWS` in dev.

---

## 9. `prune_huey_results` Correctness (🟢 Low — Correct)

**Location:** tasks.py — `prune_huey_results`

```python
flushed = huey.storage.flush_results()
huey.storage.flush_schedule()
```

This correctly bounds the `huey.db` size. The comment is accurate — Huey never prunes `taskresult` by default. The try/except with logging is appropriate. Running every 6 hours is reasonable. ✅

**Note:** `flush_results()` flushes **all** results, including those for tasks that haven't been consumed yet. Since export tasks don't use return values (the status is tracked via `ExportTask` in the app DB), this is safe. If you ever add tasks that rely on `.get()` to retrieve results, this would break them.

---

## 10. Query Building Correctness (🟢 Low — Correct)

**Location:** `services/query_service.py`

- `build_query_from_rules` uses `outerjoin(Classification)` — correct, includes Ztf alerts without classification.
- `get_query_match_count` uses `with_entities(db.func.count())` + `order_by(None)` — correct, avoids loading rows and removes unnecessary ordering.
- The count is correct with the outerjoin: `COUNT(*)` counts rows in the joined result set, and since the outerjoin doesn't duplicate Ztf rows (Classification has one row per alert_id), the count equals the number of matching Ztf alerts. ✅

---

## 11. Worker Registration (🟢 Low — Correct)

**Location:** worker.py

```python
from .extensions import huey
from . import tasks
```

This is the correct pattern. The consumer command `huey_consumer poi_broker.worker.huey` imports `worker.huey` (the shared instance), and importing `tasks` triggers `@huey.task` decorator registration on that same instance. ✅

---

## Summary Table

| # | Issue | Severity | Status |
|---|-------|----------|--------|
| 1 | Active-task check-then-create race (no atomicity, no DB constraint) | 🔴 High | Fix recommended |
| 2 | `.all()` loads entire result set — OOM risk for large exports | 🔴 High | Fix recommended |
| 3 | `file_path` leaked in status API JSON | 🟡 Medium | Fix recommended |
| 4 | Enqueue failure leaves orphaned PENDING task | 🟡 Medium | Fix recommended |
| 5 | No partial file cleanup on write failure | 🟡 Medium | Fix recommended |
| 6 | Stale reset based on `updated_at` | 🟢 Low | ✅ Correct |
| 7 | Startup double-reset | 🟢 Low | ✅ Harmless |
| 8 | Dev immediate mode | 🟢 Low | ✅ Correct |
| 9 | `prune_huey_results` | 🟢 Low | ✅ Correct |
| 10 | Query building / count | 🟢 Low | ✅ Correct |
| 11 | Worker registration | 🟢 Low | ✅ Correct |

---

## Recommended Priority Order

1. **Fix #2 (streaming)** — highest risk of production outage (worker OOM)
2. **Fix #1 (race condition)** — add partial unique index on `ExportTask`
3. **Fix #3 (file_path leak)** — trivial one-line removal
4. **Fix #4 + #5 (failure cleanup)** — enqueue rollback + partial file deletion

Would you like me to implement any of these fixes?