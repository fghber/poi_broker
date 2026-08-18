I have enough context now. Let me deliver the multi-angle audit.

# Multi-Angle Audit of the Init-Order Fix

## TL;DR

The fix is **correct for the reported bug** (init order), but the relocated `init_huey()` call still has **two latent issues** under the exact production topology you flagged: parallel Gunicorn web workers + an external Huey consumer. Both predate this diff but are now more visible because of it. None break the `pytest` suite.

---

## 1. Correctness — Is the Fix Right?

### 1a. Reported error: ✅ Fixed

The original traceback is unambiguous: `reset_stale_export_tasks()` calls `ExportTask.query.filter(...).all()` → `session.execute()` → `get_bind()` → `db.engines` → `RuntimeError("The current Flask app is not registered with this 'SQLAlchemy' instance")`.

`db.engines` is only populated **after** `db.init_app(app)`. Moving `init_huey(app)` to run **after** the `db.init_app(app)` + SQLite-PRAGMA block makes the bound app visible to the first query in `reset_stale_export_tasks`. Verified:

- `flask --app wsgi:app routes` boots clean (no `Failed to reset stale export tasks at startup`).
- 501 tests pass.

### 1b. Subtle correctness point: app-context at startup

`init_huey` calls `reset_stale_export_tasks()` inside its **own** `with app.app_context():` (line 86 of the current file), so by the time the query runs we have:

1. `app` created → `app.config` updated.
2. `db.init_app(app)` → `db.engines` populated.
3. `app.app_context()` active (the outer block around pragma setup).
4. `init_huey(app)` re-enters `app.app_context()` — redundant but harmless.

`SQLAlchemy.engines` is a per-`SQLAlchemy`-instance dict and is set the first time `init_app` runs against that instance. It does **not** depend on an active app context; it only requires that some app has been registered. So step 2 is the load-bearing step, and step 4's nested `app_context()` is just defensive. ✅

### 1c. Are there other call sites that pre-bound `db` queries on Huey init?

I checked: `init_huey()` itself only runs `reset_stale_export_tasks()`. That function is the **only** thing in this path that talks to SQLAlchemy at app boot. No other code in the moved region queries the DB (pragma setup uses engine connect listeners, which are fine before any bind is needed by `db.engines`). ✅

**Verdict:** The fix is correct.

---

## 2. Safety — What's Now Different Because of the Move?

The new ordering is:

```python
db.init_app(app)                          # ← now earlier
with app.app_context():
    event.listens_for(...)                  # pragmas
    if 'users' in db.engines: ...           # pragmas (users bind)
    _ensure_email_verification_expiry_column(app)   # ALTER TABLE user

init_huey(app)                             # ← now later (reset_stale_export_tasks)
```

Two side effects of this ordering are worth flagging.

### 2a. `_ensure_email_verification_expiry_column` runs before `init_huey` ✅ (good)

Schema migrations now happen **before** the Huey cleanup query. If a user's `users.db` is missing the `email_verification_token_expires` column, the ALTER runs first, then `reset_stale_export_tasks` queries `ExportTask`. That's fine: `ExportTask` lives in `users.db` too (`__bind_key__ = 'users'`), but it's a different table, and `ALTER TABLE user` does not lock `export_task`. Order is correct.

### 2b. `_ensure_email_verification_expiry_column` opens a transaction that overlaps with Huey init

`_ensure_email_verification_expiry_column` does:

```python
with engine.begin() as conn:
    conn.execute(text('ALTER TABLE user ADD COLUMN ...'))
```

This holds an **implicit transaction** on the `users` engine. It exits the `with` block on success, so by the time `init_huey` runs the transaction is committed. ✅ Not a regression.

### 2c. SQLite `event.listens_for(db.engine, "connect")` registered at app boot

This is fine because:

- `db.engine` is the default-bind engine. `db.engines['users']` is the users-bind engine.
- The listener is invoked per-DBAPI-connection (lazy). It runs on every future connection, including connections opened by `init_huey()` → `reset_stale_export_tasks()` → `db.session.execute(...)`. So pragmas are active for the Huey startup query too. ✅

---

## 3. Parallel-Execution Concerns (Gunicorn workers + Huey consumer)

This is the angle you asked about specifically. There are **three races** to be aware of. None are caused by this diff, but the diff makes the first one **observable** in production for the first time at startup rather than only on the first request.

### 3a. ⚠️ Multiple Gunicorn workers all run `reset_stale_export_tasks()` at boot — race-prone

**The scenario:**

- Nginx → Gunicorn `-w 4`.
- `init_huey()` calls `reset_stale_export_tasks()` synchronously, in-app, **inside an `app.app_context()`**.
- That means **every Gunicorn worker process** runs the cleanup at startup, against the **same** `users.db`.

**The current implementation of `reset_stale_export_tasks`:**

```python
stale = ExportTask.query.filter(...).all()
for task in stale:
    task.status = 'FAILED'
    task.error_message = ...
    task.updated_at = datetime.now(timezone.utc)
if stale:
    try:
        db.session.commit()
```

**Race conditions:**

1. **Lost update on `error_message` / `updated_at`.** Workers A and B both SELECT the same `stale` rows. A commits `status='FAILED', error_message="Export aborted: no progress for over 30 minute(s)."` first; B then commits the same row with the same payload. With SQLite WAL + `synchronous=NORMAL`, the second writer blocks on the B-tree page lock and either retries or gets `OperationalError: database is locked`. Two workers racing the same row is a recipe for occasional `OperationalError` during startup, which is currently **swallowed** by `init_huey`'s `except Exception` block — so it's invisible in logs beyond the `Failed to reset stale export tasks at startup` ERROR.

2. **The `EXPORT_TASK_MAX_AGE_SECONDS` window is checked against wall clock, but the `ExportTask.query.filter(... updated_at < cutoff)` + write-then-commit is not transactional isolation.** On SQLite without `BEGIN IMMEDIATE`/`EXCLUSIVE`, two readers can both see the same stale row and both try to UPDATE.

3. **Heartbeat collision** (see §3b) is unrelated but worsens this.

**Why pytest didn't catch it:** The test suite uses a single in-process app; the race only manifests with ≥2 Gunicorn worker processes running `create_app()` simultaneously.

**Mitigation options (ordered by minimalness):**

1. **Only the Huey consumer should reset stale tasks.** Move the `reset_stale_export_tasks()` call out of `init_huey()` and rely on the existing periodic `cleanup_stale_export_tasks` (every 5 min, `crontab(minute='*/5')`) to do the reset. Then the web app's startup no longer races.
   - Trade-off: between worker start and the next 5-min tick, a freshly crashed task can still 409 a user. Acceptable in practice (≤5 min).
2. **Keep the web-side call, but make `reset_stale_export_tasks` use compare-and-swap** like `_transition_export_task` already does (`UPDATE ... WHERE status IN ('PENDING','RUNNING') AND updated_at < cutoff`). Then two workers contending on the same row: the first CAS wins, the second matches zero rows, both succeed.
   - Trade-off: still some lock contention on the page, but no torn writes.
3. **Run the reset under a short-lived advisory lock** (e.g. an in-process `threading.Lock` + a one-shot flag in app extensions, *plus* a SQLite `BEGIN IMMEDIATE` if you want cross-process). The simplest is option 1.

**Recommendation:** Do option 1, then keep the periodic `cleanup_stale_export_tasks` (it already runs every 5 min via the Huey consumer). This eliminates a class of races from the web tier entirely.

If you can't remove the startup call (e.g., you want stale tasks cleared on every web boot too), at minimum convert it to CAS (option 2). The existing code already shows the pattern in `_transition_export_task`.

### 3b. ⚠️ Heartbeat race between Huey worker and Huey consumer's periodic stale reset

Unrelated to this diff but worth re-flagging because the audit asked: a worker `EXPORT_HEARTBEAT_SECONDS` (default 60 s) refreshes `updated_at` while writing; the consumer's periodic stale-reset has `STALE_TASK_MAX_AGE` (default 1800 s). With those defaults, the heartbeat races the reset only if a heartbeat is missed for ≥30 min — i.e., the worker is genuinely stuck. ✅ Healthy under defaults.

If anyone tunes `STALE_TASK_MAX_AGE` down close to `HEARTBEAT_SECONDS`, this can false-fail healthy long exports. **No change needed now**; just a config interaction to keep in mind.

### 3c. �️ Per-worker `db.session` in Flask-SQLAlchemy is **not** worker-safe

This is a broader Flask-SQLAlchemy concern, not specific to this diff, but the fix makes it more visible: each Gunicorn worker has its own `db.session` (via Flask-SQLAlchemy's scoped session keyed on the worker process). The Huey consumer is a **separate process** with its **own** `db.session`.

If a Gunicorn worker writes `ExportTask.status='PENDING'` and commits, the Huey consumer's `db.session` will see it on its next SELECT (committed → visible to other connections via WAL). ✅

If the Huey consumer writes `status='RUNNING'` and commits, the Gunicorn worker's `db.session` will see it on its next `ExportTask.query.filter_by(user_id=...).filter(status.in_(['PENDING','RUNNING'])).first()` in `export_submit` — but **only if the worker flushed/expired its session**. Flask-SQLAlchemy by default uses `Session(autoflush=True)` and identity-map caching. Across a long-running worker process, a fresh request typically gets a fresh `db.session` (request-scoped), so reads see the latest committed state. ✅ in practice, but it's worth being aware that **two long-lived sessions in different processes can drift until they touch the DB**.

The active-task uniqueness is enforced at the DB level by the partial unique index `uix_export_task_one_active_per_user` (added per `huey_audit1.md` #1), so even if both sessions have stale identity-map copies, the DB-level `IntegrityError` on insert guarantees only one `PENDING`/`RUNNING` row per user. ✅

---

## 4. Security & Safety — Does the Move Re-Open Any Window?

### 4a. CSRF / login / rate-limit init — unchanged

`csrf.init_app(app)`, `limiter.init_app(app)`, `login_manager.init_app(app)` all run **after** `init_huey(app)`. The Huey startup reset does not need or use CSRF (it's a server-internal query). ✅

### 4b. Huey storage URI warning — unchanged

The `RATELIMIT_STORAGE_URI memory://` warning still fires for multi-worker Gunicorn, and `PUBLIC_BASE_URL` warning is intact. ✅

### 4c. `secret_key` requirement — unchanged

`build_app_config` raises `RuntimeError` if `SECRET_KEY` is unset. That happens before the move. ✅

### 4d. The `try/except Exception: app.logger.exception(...)` around `reset_stale_export_tasks()`

This catch swallows **everything**, including `OperationalError: database is locked` from §3a. That's good for boot resilience (the app still starts), but it hides the race. Consider logging the exception type so production boot logs show "DB locked during stale reset" distinctly from a real code bug. Minor, optional.

---

## 5. Gaps in the Fix Description

The diff message and PR description cover **what** changed but not **why** in terms of invariants:

1. The startup reset's correctness depends on `db.init_app(app)` having been called first. The added comment ("Must run AFTER db.init_app(app) so that reset_stale_export_tasks() ... has a bound app") is good, but consider also documenting **the alternative**: this function should arguably be moved out of `init_huey` entirely, since the Huey consumer's periodic task already runs it every 5 minutes. See §3a.
2. The fix relies on the implicit invariant that **nothing else in `create_app()` between `db.init_app(app)` and `init_huey(app)` registers a new app or creates a new `SQLAlchemy` instance**. Looking at the code: no such operation exists. The only SQLAlchemy interaction in between is `db.engine` access (which works after `init_app`) and the engine event listener registration. ✅

---

## 6. Recommendations (ranked)

| # | Action | Why | Risk |
|---|---|---|---|
| 1 | **Remove `reset_stale_export_tasks()` from `init_huey()`**; rely on the periodic `cleanup_stale_export_tasks` (already runs `*/5`) | Eliminates a 4-way startup race across Gunicorn workers against a single `users.db`. Also reduces boot-time DB load. | A crashed task can 409 a user for up to 5 min after worker death. Acceptable. |
| 2 | If #1 is rejected, **convert `reset_stale_export_tasks` to compare-and-swap** using the existing `_transition_export_task` pattern (`UPDATE ... WHERE status IN ('PENDING','RUNNING') AND updated_at < cutoff`) | Eliminates torn writes; survives multi-worker boot contention. | Slightly more verbose; same behavior for healthy cases. |
| 3 | Log the **exception type** in the `except Exception` around the startup reset (`logger.warning('Stale-task reset skipped: %s', type(e).__name__)`) | Visibility into §3a races without crashing boot. | None. |
| 4 | Document that `EXPORT_STALE_MAX_AGE_SECONDS` should remain `≥ 2 × EXPORT_HEARTBEAT_SECONDS` | Prevents the §3b false-fail interaction if anyone tunes either value. | None. |

None of #1–#4 are required to land the current fix. The current diff correctly resolves the reported bug. #1 is the cleanest follow-up because it pushes a process-shared responsibility to the **single** process that owns it (the Huey consumer).

---

## 7. Sanity-check on the Audit Itself

I confirmed:

- `__init__.py` lines 121–155 contain the fix exactly as described.
- `tasks.py` `reset_stale_export_tasks()` is the only thing `init_huey()` calls that touches the DB.
- `extensions.py` and `worker.py` are unchanged.
- The repo-memory note at `/memories/repo/huey-export.md` (line 64) already documents the partial unique index + `IntegrityError→409` race fix from a prior audit; that fix **complements** this one (the DB-level invariant is what protects cross-process identity-map drift even if the session cache is stale).

No further code change is required to land this PR. Apply recommendations as follow-up tickets if desired.