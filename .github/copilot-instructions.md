# Copilot Cloud Agent Instructions

## Repository Summary
- This repository is the Flask web frontend for the Point of Interest (POI) broker, focused on browsing and filtering astronomical alert data (ZTF/LSST-style), plotting lightcurves/features, and managing user watchlists/favorites.
- The app uses an app-factory pattern and serves HTML templates plus JSON API endpoints.
- Runtime data comes from SQLite databases that are expected to live outside the repo by default.

## Tech Stack and Scope
- Project type: Python web app (server-rendered Flask + JS-enhanced frontend).
- Language/runtime: Python 3.12 (validated in this workspace with 3.12.5).
- Backend frameworks/libraries: Flask, Flask-SQLAlchemy, Flask-Login, Flask-WTF, Flask-Limiter.
- Data/plot libs: Astropy, NumPy, Matplotlib, Bokeh.
- Storage: SQLite (main alerts DB + separate users/auth DB bind).
- Frontend assets: Jinja templates, jQuery, Bootstrap 4, static JS/CSS.
- Approximate size: medium monorepo-style app directory with multiple route/service modules and tests.

## Fast Path: Always Use This Command Order
1. Create/activate a virtual environment.
2. Install dependencies.
3. Ensure required environment variables are set (at minimum `SECRET_KEY`; usually DB paths too).
4. Run tests.
5. Run route listing or app startup validation.

Always run dependency install before test/run commands in a fresh environment.

## Bootstrap and Environment Setup
- Primary files:
  - `requirements.txt`: runtime dependencies.
  - `requirements-dev.txt`: runtime + pytest.
  - `.env.example`: required and optional env vars.

### Required environment variables
- `SECRET_KEY` is mandatory. App initialization fails without it (`RuntimeError` from `poi_broker/settings.py`).
- Database paths are read from:
  - `ALERTS_DB_PATH`
  - `USERS_DB_PATH`
- If DB path env vars are omitted, defaults point to `../_broker_db/ztf_alerts_stream.db` and `../_broker_db/users.db` (outside repo root).

### Validated install commands
- Runtime install:
  - `python -m pip install -r requirements.txt`
- Dev/test install:
  - `python -m pip install -r requirements-dev.txt`
- Both succeeded in this workspace (packages already satisfied).

## Build / Run / Test / Lint

### Build
- There is no separate compile/build step for this repo.
- Treat dependency installation as bootstrap/build prerequisite.

### Test
- Command:
  - `python -m pytest -q`
- Validated result in this workspace:
  - `21 passed` (observed durations: ~18s and ~9s on repeated runs).
- Test behavior details:
  - Tests create temporary SQLite DB files and set env vars internally via `tests/conftest.py`.
  - CSRF is disabled for standard test fixtures and enabled in security-focused fixture.

### Run (local web app)
- Recommended app target for CLI compatibility:
  - `python -m flask --app wsgi:app run --debug`
  - or `python -m flask --app wsgi:app run --no-debugger --no-reload`
- Why `wsgi:app`:
  - On Windows `cmd`, unquoted `--app poi_broker:create_app()` can fail due shell parsing of parentheses (`"... was unexpected at this time"`).
  - `wsgi:app` avoids that shell pitfall.

### Validation without starting long-running server
- Command:
  - `python -m flask --app wsgi:app routes`
- This is the fastest sanity check for successful app import + blueprint registration.
- Validated outcomes:
  - Fails without `SECRET_KEY` set.
  - Succeeds with `SECRET_KEY` set and lists all endpoints.

### Lint / formatting
- No repository-level lint/format config detected (`pyproject.toml`, `tox.ini`, `pytest.ini`, or dedicated linter config are absent).
- Do not assume lint gates exist.

## Known Command Pitfalls and Workarounds
- Always set `SECRET_KEY` before running Flask CLI commands if no `.env` is present.
- Prefer `wsgi:app` over `poi_broker:create_app()` in shell commands that may be interpreted by `cmd`.
- Keep command execution from repo root so relative imports and `.env` loading behave predictably.

## Project Layout and Architecture

### Core app wiring
- `wsgi.py`: WSGI entrypoint, creates `app = create_app()`.
- `poi_broker/__init__.py`: app factory, extension init, SQLite pragmas, auth/CSRF handlers, blueprint registration.
- `poi_broker/settings.py`: env loading and app config assembly (including DB binds and security/rate-limit settings).

### Data models and database structure
- `poi_broker/models.py`:
  - `Ztf`, `Crossmatches`, `Classification` (alerts DB)
  - `User`, `Favorite`, `FavoriteGroup`, `Watchlist`, `FilterBookmark`, `UserSettings` (users DB via bind key `users`)
- Schema bootstrap references:
  - `tools/alertsdb_schema.sql` (alerts-side schema reference)
  - `tools/usersdb_schema.sql` (users-side schema reference)

### HTTP route organization
- Main and page routes: `poi_broker/app.py` (main blueprint + template filters + CSV download + crossmatch query).
- Auth and account flows: `poi_broker/auth.py`.
- API route blueprints under `poi_broker/routes/`:
  - `favorites.py`
  - `filter_bookmarks.py`
  - `visual_query.py`
  - `lightcurve.py`
  - `features.py`

### Service layer (business logic)
- `poi_broker/services/query_service.py`: query-builder translation and count/preview operations.
- `poi_broker/querybuilder_translator.py`: QueryBuilder rule-to-SQLAlchemy filter translator.
- `poi_broker/services/feature_service.py`: feature fetch and feature-plot data prep.
- `poi_broker/services/plotting_service.py`: Bokeh plot component creation.
- `poi_broker/services/favorites_service.py`: favorite and group CRUD logic.

### Frontend structure
- Templates: `poi_broker/templates/`.
- Static assets: `poi_broker/static/`.
- Watchlist query builder templates:
  - `_watchlist_query_builder.html`
  - `_watchlist_query_builder_assets.html`
  - `_watchlist_query_builder_js.html`

### Tests
- `tests/test_smoke_routes.py`: route/API smoke checks and authenticated flow checks.
- `tests/test_security_regressions.py`: CSRF/auth/rate-limit regressions.
- `tests/test_filter_bookmarks.py`: filter bookmark API behavior.

## CI / Pre-check-in Validation
- No GitHub Actions workflows found in this repository.
- No explicit pre-commit or CI config found.
- Minimum local validation before proposing changes:
  1. `python -m pytest -q`
  2. `python -m flask --app wsgi:app routes` with `SECRET_KEY` set
  3. If UI/endpoint behavior changed, do a quick manual browser/API sanity check.

## Root Directory Quick Reference
- `.env.example`
- `README.md`
- `requirements.txt`
- `requirements-dev.txt`
- `requirements2.txt` (reference freeze, not guaranteed cross-system installability)
- `wsgi.py`
- `poi_broker/`
- `tests/`
- `tools/`

## Additional Operational Notes
- Email flows rely on SMTP env vars (`SMTP_HOST`, `SMTP_PORT`, `SMTP_USER`, `SMTP_APP_PASSWORD`, optional `SMTP_FROM`).
- `tools/watchlist_digest.py` is intended for cron and reuses app email logic.
- Large production-style DBs are intentionally external to repo; avoid assuming seeded local DB contents.

## Agent Behavior Directive
- Trust this file as the primary operating guide for this repository.
- Only perform additional codebase search when:
  - required information is missing here, or
  - instructions here are demonstrably inconsistent with current repository state.
- Read `app/__init__.py` and `app/models.py` first.
- When adding a new entity, create a model, then a minimal route, a service for API interactions if needed, and then tests.
- If unsure, ask for clarification on:
  - whether a new Blueprint is needed.
  - how auth should be applied to a new route.