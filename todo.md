# Fixes

- [ ] Diagnose the confirmed large-CSV download failure before changing export initialization or delivery configuration. Two earlier failures were reported; the 2026-10-12 live test reproduced a stopped download for task 11 (`locus_ra < 260`, 741,513 preview rows) after the site reported successful background generation. Trace `/export/download/11` through nginx/Gunicorn and verify the server file before selecting a fix. Preserve the earlier **High: Huey initialization order** flag as a locally reproduced concern whose claimed synchronous production behavior was contradicted by this live test. See [large CSV download diagnosis](docs/async_export/large_csv_download_diagnosis.md).

- [x] change order of Summary: Watchlist <-> Bookmarks
- [x] Rename "Saved Filters" to Bookmarked Filters

# Cleanup

- [ ] Remove deferred-hygiene items from the 2026-09 review (no runtime impact, revisit opportunistically):
  - unused imports: `auth.py` (_format_expire_time), `filter_service.py` (Optional), `search_service.py` (List/Union/ParsedValue)
  - dead vendor CSS: `poi_broker/static/css/bootstrap.4org.css` (~10k lines no template references)
  - `tools/logwatch.py`: scratch script (commented-out code, undeclared `watchdog` dependency) — delete or move out of the repo
  - `services/query_service.py`: hardcoded `classification.` prefix duplicating the models_dict table list (drift risk)
- [ ] Remove `tools/apply_users_indexes.sql` once the one-active-export-per-user fix (`uix_export_task_one_active_per_user`) is committed, merged into the parent main, and applied in production — see `docs/adr/adr_db_schema_migration_strategy.md`

- [x] align jQuery/Bootstrap versions and integrity tags
- JSON responses differs on error types and escaping: Some endpoints use jsonify, others build Response/current_app.response_class with json.dumps or safe_serialize (query_features, query_crossmatches).
- [x] request.query_string.decode('ascii'): Non-ASCII query strings can raise; utf-8 with errors policy is safer.
- [x] implement `get_flashed_messages(with_categories=True)` globally for all categoires used: 'danger', 'success', 'info',  'warning'
- [x] apply db.session.commit() pattern when
  * Violating a unique constraint (email column)
  * Violating a foreign key constraint
  * Violating a NOT NULL constraint
- [x] Logging setup in `create_app()` (poi_broker/__init__.py:86-89, NOT app.py:104-107 as flagged): `logging.basicConfig` with hardcoded CWD-relative `app.log` (mode 'a+'), no config override. Boot depends on CWD writability (systemd WorkingDirectory/ReadWritePaths currently makes it work); pytest runs append to ./app.log in repo root; lands in /opt/poi_broker while huey.log goes to /var/log via HUEY_LOGFILE; no rotation. Fix: resolve path from `base_dir`/env var (note: basicConfig runs before `base_dir` is computed at line 97 — reorder needed), consider RotatingFileHandler. Not urgent: works under documented deployment.
  - Fixed via `_configure_logging(base_dir)` in `poi_broker/__init__.py`: path = `APP_LOG_FILE` env var, else `<workspace root>/app.log` (CWD-independent, same location as before in dev and under systemd); `RotatingFileHandler` (5 MB x 3); stderr fallback if unwritable; skipped entirely when `FLASK_TESTING` so pytest never writes `app.log`.

# Considerations

- [x] consider extending rate-limiting to heavy read routes/endpoints (/query_features, /query_crossmatches, /download_alerts_csv, main /)
- [x] use a single grouped query to get_favorite_groups (favorites_service.py): For each group it runs Favorite.query.filter_by(group_id=g.id).count() — classic N+1. 
- [x] Consider refactoring the spinner into a reusable JavaScript module.
- Update main table data via AJAX/API calls instead of page loads/GET

# New Features

- [x] Save table filters (URL) as bookmark (My Search/Filter)
- [x] Allow users selecting (up to 10) default features to plot
- [x] Allow users to create custom observatory coordinates for the observing planning tool
- [x] Add Default Observatory Coordinates -> Last-used becomes the default for the next session
- [x] Document new features in the README.md
- [x] Add Bulk Export based on Visual Query (Top 1M/Preview or All)
  - [x] Create CSV fully async, inform user when ready

# Unproven / high-risk (repro-gate, 2026-09-15)

Not confirmed as bugs. Either the candidate did not reproduce, or it is a documented accepted risk. Do not treat these as findings without a failing test.

- Observing plot `dec=999`: returns 200 “not visible” (`|lat−dec| ≥ 90`). Correct, not a 500.
- Remember-me cookie not bound to IP/User-Agent: documented trade-off in `docs/sec.md`.
- CSP `'unsafe-inline'` / `'unsafe-eval'`: documented Bokeh/jQuery compromise.
- `ProxyFix(..., x_host=1)`: documented; `PUBLIC_BASE_URL` is the email-host control.
- Classification plot echoing `alertId` into Bokeh when a row exists: missing-row path is escaped and tested; XSS via a stored matching `alert_id` was not reproduced.
- Lightcurve CSV built with f-strings: possible CSV injection if `locus_id` has commas/`=`; no failing test.
- Unbounded query-builder nesting: possible CPU/DoS; not executed as a crash.
- Huge `page` OFFSET: slow/empty, not shown to 500.
- Observing plot `ra=inf` / `nan`: not re-tested after the calendar-date fix.
- `/query_lightcurve_data` has no 128-char id cap (unlike features/classification): inconsistency, not a crash.
- `query_features` 404 on missing alert: intentional data API, not a plot empty-state.
- Rate-limit `memory://` under multi-worker Gunicorn: documented ops risk.

# Future

- Migrate to more capable DB (PostgreSQL)
- Change column type: ant_magband REAL -> TEXT
  ```
  ALTER TABLE featuretable ADD COLUMN ant_passband_str TEXT
  UPDATE featuretable SET ant_passband_str = CAST(ant_passband as TEXT)
  ALTER TABLE featuretable DROP COLUMN ant_passband;
  ALTER TABLE featuretable RENAME COLUMN ant_passband_str TO ant_passband;

  Recommended (safer, canonical SQLite way)
  CREATE TABLE new_featuretable (
        -- same schema, but:
        ant_passband TEXT,
        ...
    );
    INSERT INTO new_featuretable (...)
    SELECT
        CAST(ant_passband AS TEXT),
        ...
    FROM featuretable;
    DROP TABLE featuretable;
    ALTER TABLE new_featuretable RENAME TO featuretable;
    HOWEVER: SQLite happily stores TEXT in a REAL column (Type Affinity). This only matters when moving to PosgreSQL/MySQL
  ```