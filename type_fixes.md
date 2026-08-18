# Follow-up Plan: Pre-existing Pylance Type-Checking Noise

**Status:** Follow-up (non-blocking). These are **pre-existing** Pylance diagnostics, not
introduced by the Huey init-order / stale-cleanup work. They do not affect runtime —
the full suite passes (498) and `flask --app wsgi:app routes` boots cleanly. This plan
tracks cleaning them up so the Problems panel is quiet and the code is more type-safe.

All three categories below are cosmetic/type-hygiene. None change behavior. Each is
independent and can be done in its own small PR.

---

## 1. `from .app import register_blueprints` shadows the `app` variable

**File:** `poi_broker/__init__.py` (inside `create_app()`)

**Problem:** The local variable `app = Flask(__name__)` is shadowed by the module
import `from .app import register_blueprints`. Pylance then resolves `app.register_blueprint`,
`app.context_processor`, `app.after_request`, `app.errorhandler`, `app.config`, etc. against
the **module** `.app` instead of the `Flask` instance, producing a cascade of
`"register_blueprint" is not a known attribute of module ".app"` errors.

**Fix options (pick one):**

- **Option A (recommended, minimal):** Import the function under an alias so it no longer
  collides with the `app` variable:
  ```python
  from .app import register_blueprints as _register_blueprints
  ...
  _register_blueprints(app)
  ```
- **Option B:** Move the blueprint registration into a helper that takes `app` as a
  parameter and is imported at module top (outside `create_app`), avoiding the local
  shadow entirely.

**Scope of the fix:** Resolves the bulk of the `__init__.py` errors (lines ~181–237):
`register_blueprints`, `app.register_blueprint`, `login_manager.login_view`,
`context_processor`, `config`, `after_request`, `errorhandler`.

**Validation:** `get_errors` on `poi_broker/__init__.py` should show far fewer
`module ".app"` errors; `python -m pytest -q` still green; `flask --app wsgi:app routes`
still boots.

---

## 2. Legacy `Model.query` usage (SQLAlchemy 2.0 deprecation)

**Files:** `tests/test_huey_tasks.py` (lines ~286, 292), `tests/test_export_routes.py`
(line ~93), and any other `Model.query.get(...)` / `Model.query` call sites.

**Problem:** `Query.get()` is legacy in SQLAlchemy 2.0 and emits a `LegacyAPIWarning`
("The Query.get() method is considered legacy... now available as Session.get()").
Pylance also flags `"status" is not a known attribute of "None"` because `Query.get()`
returns `Optional`.

**Fix:** Replace `Model.query.get(pk)` with `db.session.get(Model, pk)`:
```python
# before
assert ExportTask.query.get(task.id).status == 'RUNNING'
# after
assert db.session.get(ExportTask, task.id).status == 'RUNNING'
```
For the `"status" is not a known attribute of "None"` diagnostics, either assert
non-None first or use `db.session.get(...)` and add a `assert task is not None`.

**Scope:** Sweep the repo for `\.query\.get\(` and `\.query\.` legacy patterns. Prefer
`db.session.get(Model, pk)` / `db.session.scalars(...)`.

**Validation:** `python -m pytest -q` green; the `LegacyAPIWarning` summary disappears
from the test run.

---

## 3. `db.engines` dict-access style

**File:** `poi_broker/__init__.py` — `_ensure_email_verification_expiry_column` (line ~29)

**Problem:** Pylance suggests using `dict.get()` instead of an `if "key" in dict` +
index pattern:
```python
engine = db.engines['users'] if 'users' in db.engines else db.engine
```

**Fix:**
```python
engine = db.engines.get('users', db.engine)
```

**Validation:** `get_errors` on `poi_broker/__init__.py` no longer reports the
`Use db.engines.get(...)` suggestion at line 29.

---

## 4. `_transition_export_task` `update()` values typing (optional)

**File:** `poi_broker/tasks.py` (line ~80)

**Problem:** Pylance flags the `dict[str, Unknown]` passed to `Query.update(values, ...)`
because the `values` parameter is typed `Dict[_DMLColumnArgument, Any]` (invariant key
type). This is a type-annotation mismatch, not a runtime issue.

**Fix (optional):** Annotate the `values` parameter with the expected column-key type,
e.g. `dict[str, Any]` won't satisfy the invariant; instead build the dict with
`Column` keys or cast. Given the code already works and is covered by tests, this is
low priority — leave as-is unless the Problems panel noise is bothersome.

---

## Suggested order

1. **#1** (biggest win — clears the `__init__.py` cascade).
2. **#2** (removes deprecation warnings from every test run).
3. **#3** (one-line, trivial).
4. **#4** (optional, only if desired).

Each is independent; none are required for the current Huey fix to be correct.
