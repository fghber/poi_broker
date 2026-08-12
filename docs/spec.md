# POI Broker — Project Specification (Spec-Anchored Development Reference)

**Version:** 1.0.0
**Last updated:** 2026-08-03
**Status:** Draft — derived from codebase inspection (app v3.1.0)
**Scope:** HTTP endpoints only (internal service-layer contracts intentionally excluded)

---

## 1. Project Overview

### 1.1 Purpose
The **Point of Interest (POI) Community Broker** is a transient alert software (Rubin Observatory / LSST Alert and Community Broker), currently tested against the ZTF alert stream. This repository contains the **web frontend** for the broker: browsing and filtering astronomical alert data, plotting lightcurves and features, and managing user watchlists/favorites.

### 1.2 Scope and Boundaries
- **In scope:** Web frontend — server-rendered Flask app + JSON API endpoints, user auth/account management, alert browsing/filtering, visual query builder, lightcurve/feature/classification plotting, favorites/watchlists, filter bookmarks, user observatories, observing tool, CSV export.
- **Out of scope:** Alert ingestion, annotation, classification, and forwarding pipelines (backend broker processing). The main alerts SQLite DB is intentionally external to the repo.

### 1.3 Key Stakeholders
- Astronomers (variable-star observers) — primary users.
- Vanderbilt University development team.
- Rubin Observatory / LSST community broker program.

### 1.4 Tech Stack
- Python 3.12, Flask (app-factory pattern), Flask-SQLAlchemy, Flask-Login, Flask-WTF, Flask-Limiter.
- Astropy, NumPy, Matplotlib, Bokeh.
- SQLite (main alerts DB + separate users/auth DB bind).
- Jinja2 templates, jQuery, Bootstrap 4.

---

## 2. Source-of-Truth Spec Inventory

| Source | Path | Authoritative For | Notes |
|---|---|---|---|
| This spec | `docs/spec.md` | System behavior, API contracts, requirements | Primary reference for spec-anchored development |
| App config | `poi_broker/settings.py` | Env vars, rate limits, security config | `APP_VERSION='3.1.0'` |
| Env template | `.env.example` | Required/optional env vars | `SECRET_KEY` mandatory |
| Data models | `poi_broker/models.py` | Entities, columns, relationships, binds | Alerts + users DB |
| Alerts schema | `tools/alertsdb_schema.sql` | Alerts-side schema reference | |
| Users schema | `tools/usersdb_schema.sql` | Users-side schema reference | |
| Tests | `tests/*.py` | Enforced behavior (de-facto spec) | See §8 |
| Coverage policy | `tests/coverage.md` | Coverage targets | 50–75%/module, 75%+ total |
| Repo map | `Repomap.md` | Module/route overview | |
| README | `README.md` | Install/usage/ops guidance | Operational, not behavioral spec |

> **Note:** No formal external spec document exists. All requirements below are **inferred from code** unless a test explicitly enforces them (marked accordingly).

---

## 3. Requirements

### 3.1 Functional Requirements

| ID | Requirement | Enforced By |
|---|---|---|
| REQ-001 | User can sign up with email/name/password; email deliverability validated, min 8-char password | `tests/test_auth.py` |
| REQ-002 | User verifies email via emailed token (SHA-256 hashed) | `tests/test_auth.py` |
| REQ-003 | User can log in with email/password (optional remember-me) | `tests/test_auth.py` |
| REQ-004 | User can log out (secure session cookie cleared) | `tests/test_auth.py` |
| REQ-005 | User can request password reset (1-hour expiry token) and set new password | `tests/test_auth.py` |
| REQ-006 | User can change password (requires current password) | `tests/test_auth.py` |
| REQ-007 | User can browse/filter alerts on main page (date/MJD, alert_id, object_id, passband, locus_id, RA/Dec, magnitude, prob_class; multiple sort keys; 100 rows/page) | `tests/test_app.py`, `tests/test_alert_id_filter.py`, `tests/test_dec_filter.py`, `tests/test_ant_passband_filter.py` |
| REQ-008 | User can download filtered alerts as CSV | `tests/test_app.py` |
| REQ-009 | User can query crossmatches for a locus | `tests/test_app.py` |
| REQ-010 | User can build visual queries (QueryBuilder) and preview SQL / get match count | `tests/test_visual_query*.py` |
| REQ-011 | User can save/delete named watchlists (private, per-user) | `tests/test_visual_query*.py` |
| REQ-012 | User can view lightcurve plot and export locus CSV | `tests/test_lightcurve*.py` |
| REQ-013 | User can view feature values for an alert and feature plots for a locus | `tests/test_features*.py` |
| REQ-014 | User can toggle favorites per locus and organize into groups | `tests/test_favorites.py` |
| REQ-015 | User can save/delete named filter bookmarks (whitelisted keys) | `tests/test_filter_bookmarks.py` |
| REQ-016 | User can CRUD custom observatories (name/lat/lon, auto timezone) | `tests/test_user_observatories*.py` |
| REQ-017 | User can view observing plot (AltAz visibility) for builtin/custom observatory | `tests/test_observing_tool*.py` |
| REQ-018 | User can view classification radar chart for an alert | `tests/test_classification.py` |
| REQ-019 | User can configure default feature-plot columns and last-selected observatory in settings | `tests/test_user_settings*.py` |

### 3.2 Non-Functional Requirements

| ID | Requirement | Value / Source |
|---|---|---|
| NFR-001 | Auth rate limiting | login 10/min, signup 5/hr, forgot-pw 5/hr, reset-pw 10/hr, change-pw 10/hr (`AUTH_RATE_LIMIT_*`) |
| NFR-002 | Read rate limiting | `/`, `/query_crossmatches`, `/query_features` 30/min (LAX); `/download_alerts_csv` 15/min (MEDIUM) |
| NFR-003 | CSRF protection | Flask-WTF `CSRFProtect` on all POST forms; JSON 400 for `/api/*` |
| NFR-004 | Security headers | HSTS (HTTPS), `nosniff`, `X-Frame-Options: DENY`, Referrer-Policy, Permissions-Policy, strict CSP |
| NFR-005 | Cookie security | HttpOnly, Secure (prod), SameSite=Lax; remember cookie 14 days |
| NFR-006 | Password hashing | Werkzeug `generate_password_hash`/`check_password_hash` |
| NFR-007 | SQLite performance | WAL, `synchronous=NORMAL`, 64MB cache, 256MB mmap, `temp_store=MEMORY` |
| NFR-008 | Caching | `lru_cache` on MJD formatting (8192) and builtin observatories (1); per-request `UserSettings` reuse |
| NFR-009 | Coverage targets | 50–75% per module, 75%+ total (`tests/coverage.md`) |
| NFR-010 | Proxy support | `ProxyFix` enabled in production (trust `X-Forwarded-*`) |

### 3.3 Performance Notes

The main page (`start()`) loads an authenticated user's `UserSettings` row **once** per request and reuses it across the feature-plot-column and observatory-selection lookups, avoiding a redundant database query on every main-page load.

- `poi_broker/user_settings.py` exposes `get_user_settings(user_id)` to load the row once.
- `get_saved_feature_plot_columns()` and `get_saved_last_selected_observatory()` accept an optional pre-loaded `settings` row; when omitted they fall back to loading it themselves (backward compatible).
- This is a **per-request** optimization only — there is no cross-request cache, so no invalidation is needed when a user updates their settings in the profile UI. Each request reads the current row from the database.

---

## 4. Data Model

Two SQLite databases via SQLAlchemy binds: **alerts** (default bind) and **users** (`__bind_key__='users'`).

### 4.1 Alerts DB
| Model | Table | Key Columns / PK | Notes |
|---|---|---|---|
| `Ztf` | `featuretable` | Composite PK `(date_alert_mjd, alert_id, locus_id)` | ~200 feature columns; `locus_ra`, `locus_dec`, `ant_mag_corrected`, `ant_passband`, per-band `_magn_r/_g`, `_flux_r/_g`, `anomaly_score/mask/type`, `g_r_max/mean`, `best_period*`, `power_rate_*`, `mhps_*`, `drw_*`; rel `classification` (1:1 via alert_id, viewonly) |
| `Crossmatches` | `crossmatches` | `id` PK | `locus_id`, `catalog`, `object`, `ra_cat`, `dec_cat`, `separation` |
| `Classification` | `classification` | `alert_id` PK | `p_cvnova`, `p_e`, `p_lpv`, `p_puls`, `p_periodic_other`, `p_quas`, `p_sn`, `p_yso`, `prob_class` |

### 4.2 Users DB
| Model | Table | Key Columns / Constraints | Notes |
|---|---|---|---|
| `User` | `user` | `id` PK, `email` unique, `password` (hashed), `name`, `role` (default `'user'`), `email_verified`, `email_verification_token`, `reset_token`, `reset_token_expires` | Methods `has_role`, `is_admin`; decorator `role_required(role)` |
| `FavoriteGroup` | `favorite_group` | `id`, `user_id` FK→user CASCADE, `name`, `created_at`; unique `(user_id, name)` | |
| `Favorite` | `favorite` | `id`, `user_id` FK CASCADE, `locus_id`, `group_id` FK→favorite_group SET NULL, `created_at`; unique `(user_id, locus_id)` | |
| `Watchlist` | `watchlist` | `id`, `user_id` FK CASCADE, `name`, `rules_json` (Text), `sql_where` (Text), `created_at` (epoch); unique `(user_id, name)` | |
| `FilterBookmark` | `filter_bookmark` | `id`, `user_id` FK CASCADE, `name`, `query_json` (Text), `created_at` | |
| `UserObservatory` | `user_observatory` | `id`, `user_id` FK CASCADE, `name` (max 100, NOCASE), `latitude`, `longitude`, `timezone_name`, `created_at`; unique `(user_id, name)` | |
| `UserSettings` | `user_settings` | `id`, `user_id` FK CASCADE, `default_feature_plot_columns` (JSON Text), `last_selected_observatory_json` (JSON Text); unique `user_id` | |

---

## 5. API Contracts (HTTP Endpoints Only)

### 5.1 Main Blueprint (`poi_broker/app.py`)
| Method | Path | Auth | Request | Response | Errors |
|---|---|---|---|---|---|
| GET | `/` | Public | Query: `page`, `date`, `date_alert_mjd`, `alert_id`, `ztf_object_id`, `ant_passband`, `locus_id`, `locus_ra`, `locus_dec`, `magpsf`, `prob_class`, `sort__*` | HTML `main.html` (100 rows/page) | — |
| GET | `/help` | Public | — | HTML `help.html` | — |
| GET | `/contact` | Public | — | HTML `contact.html` | — |
| GET | `/profile` | Login | — | HTML `profile.html` | — |
| GET | `/download_alerts_csv` | Public | Query: `alert_id` (repeatable) | CSV | 400 missing, 404 no records, 500 error |
| GET | `/query_crossmatches` | Public | Query: `locusId` | JSON array | 400 missing, 500 error |

### 5.2 Auth Blueprint (`poi_broker/auth.py`)
| Method | Path | Auth | Request | Response | Errors |
|---|---|---|---|---|---|
| GET | `/login` | Public | `?forgot_password=true` | HTML `login.html` | — |
| POST | `/login` | Public | Form: `email`, `password`, `remember` | Redirect → `main.profile` | — |
| GET | `/signup` | Public | — | HTML `signup.html` | — |
| POST | `/signup` | Public | Form: `email`, `name`, `password` | Sends verification email | — |
| GET | `/verify-email/<token>` | Public | Token | Verifies email | — |
| POST | `/logout` | Login | — | Clears session | — |
| GET | `/forgot-password` | Public | — | HTML `forgot_password.html` | — |
| POST | `/forgot-password` | Public | Form: `email` | Sends reset email (1h) | — |
| GET | `/reset-password/<token>` | Public | Token | HTML `reset_password.html` | — |
| POST | `/reset-password/<token>` | Public | Form: new password | Sets new password | — |
| GET | `/security` | Login | — | HTML `security.html` | — |
| POST | `/change-password` | Login | Form: `current_password`, `new_password`, `new_password_confirm` | Changes password | — |

### 5.3 Favorites Blueprint (`poi_broker/routes/favorites.py`, prefix `/api`) — all Login
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/api/favorite` | Query: `locusId` | `{"fav": bool}` | — |
| GET | `/api/favorites` | Query: `groupId` (or null) | `{"favorites":[{"id","locusId"}]}` | — |
| POST | `/api/favorite` | JSON: `{"locusId","fav","groupId"}` | `{"status":"ok"}` | 400/401/404/409/500 |
| PATCH | `/api/favorite/<int:favorite_id>/group` | JSON: `{"groupId"}` | `{"status":"ok","groupId"}` | — |
| GET | `/api/favorite-groups` | — | `{"groups":[{"id","name","count"}]}` (incl. "Ungrouped") | — |
| POST | `/api/favorite-groups` | JSON: `{"name"}` | `{"status":"ok","id","name"}` 201 | — |
| DELETE | `/api/favorite-groups/<int:group_id>` | — | `{"status":"ok"}` (orphans favorites) | — |

### 5.4 Filter Bookmarks Blueprint (`poi_broker/routes/filter_bookmarks.py`, prefix `/api`) — all Login
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/api/filter-bookmarks` | — | `{"filterBookmarks":[{"id","name","params","path","created_at"}]}` | — |
| POST | `/api/filter-bookmarks` | JSON: `{"name","params"}` | `{"status":"ok",...}` 201 | 409 duplicate name |
| DELETE | `/api/filter-bookmarks/<int:bookmark_id>` | — | `{"status":"ok"}` | 404 |

> Whitelist: `ALLOWED_KEYS` (18 keys); `MAX_BOOKMARK_JSON_BYTES=16384`; `MAX_PARAM_VALUE_LEN=512`.

### 5.5 Visual Query Blueprint (`poi_broker/routes/visual_query.py`) — all Login
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/visual_query` | — | HTML `visual_query.html` | — |
| POST | `/api/preview-query` | JSON rules | `{"sql":"..."}` | 400 invalid rules |
| POST | `/api/export-query` | JSON rules | `{"count": int}` | — |
| POST | `/api/watchlist` | JSON: `{"name","rules"}` | `{"status":"ok","id","name"}` 201 | 409 dup, 406 invalid |
| GET | `/api/watchlist` | — | `{"watchlists":[{"id","name","sql_where","created_at"}]}` | — |
| DELETE | `/api/watchlist/<int:wid>` | — | `{"status":"ok"}` | 404 |

### 5.6 Lightcurve Blueprint (`poi_broker/routes/lightcurve.py`) — Public
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/query_lightcurve_data` | Query: `locusId` | Bokeh `div+script` HTML | 400/500 |
| GET | `/locus_plot_csv` | Query: `locusId` | CSV `locus_id,date_alert_mjd,ant_mag_corrected` | — |

### 5.7 Features Blueprint (`poi_broker/routes/features.py`) — Public
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/query_features` | Query: `alert_id` | JSON of all feature values | 400/404/500 |
| GET | `/query_featureplot_data` | Query: `locusId`, `features` | Bokeh `div+script` HTML | 400/500 |

### 5.8 User Observatories Blueprint (`poi_broker/routes/user_observatories.py`, prefix `/api`) — all Login
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/api/user-observatories` | — | `{"userObservatories":[{"id","name","latitude","longitude","timezone_name","created_at"}]}` | — |
| POST | `/api/user-observatories` | JSON: `{"name","latitude","longitude"}` | `{"status":"ok",...}` 201 | 409; idempotent 200 on duplicate |
| DELETE | `/api/user-observatories/<int:observatory_id>` | — | `{"status":"ok"}` | idempotent 200 `already_deleted`; may include `warning`/`fallback` |

### 5.9 User Settings Blueprint (`poi_broker/user_settings.py`) — all Login
| Method | Path | Request | Response |
|---|---|---|---|
| GET | `/settings` | — | HTML `user_settings.html` |
| POST | `/settings` | Form: `default_feature_plot_columns[]` | Saves default feature plot columns |

### 5.10 Observing Tool Blueprint (`poi_broker/observing_tool.py`) — Public
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/query_observing_plot` | Query: `obs_loc`, `obs_date`, `obs_tz`, `ra`, `dec` | Matplotlib plot (base64/HTML) | 401 for custom observatory without auth |

### 5.11 Classification Blueprint (`poi_broker/classification.py`) — Public
| Method | Path | Request | Response | Errors |
|---|---|---|---|---|
| GET | `/query_classification` | Query: `alertId` | Bokeh radar chart | — |

---

## 6. Workflows

### 6.1 Authentication Lifecycle
Signup → email verification (SHA-256 token) → login (Flask-Login) → authenticated browsing. Password reset via emailed token (1h expiry). Change password requires current password. Role-based access via `role_required`.

### 6.2 Browse & Filter Alerts
`GET /` builds a `Ztf` query with optional filters (date/MJD, alert_id prefix/full, object_id, passband, locus_id, RA/Dec ranges, magnitude, prob_class) and multiple sort keys. Uses `SearchService` + `FilterService`. Paginated 100/page.

### 6.3 Visual Query → Watchlist
`GET /visual_query` → QueryBuilder UI → `POST /api/preview-query` (SQL preview) → `POST /api/export-query` (match count) → `POST /api/watchlist` (save rules + generated SQL). Uses `QueryService` + `querybuilder_translator.py`.

### 6.4 Lightcurve / Features / Classification
- Lightcurve: `GET /query_lightcurve_data` (Bokeh) + `/locus_plot_csv`.
- Features: `GET /query_features` (values) + `/query_featureplot_data` (Bokeh).
- Classification: `GET /query_classification` (radar chart).

### 6.5 Favorites & Groups
Toggle favorite per locus → organize into groups → list with counts. Deleting a group orphans its favorites (group_id SET NULL).

### 6.6 Filter Bookmarks
Save/restore main-table filter param sets (whitelisted keys) as named bookmarks.

### 6.7 Observatories & Observing Tool
CRUD custom observatories (auto timezone via `TimezoneFinder`); last selection persisted in `UserSettings`. Observing tool computes AltAz visibility from builtin/custom observatory.

### 6.8 Edge Cases
- Duplicate names → 409 (favorite groups, watchlists, filter bookmarks, observatories).
- Group delete → favorites orphaned (not deleted).
- Observatory delete → idempotent 200 `already_deleted`; may return `warning`/`fallback`.
- Watchlists strictly private (per-user).

---

## 7. Feature-to-Spec Mapping

| Feature | Spec Requirement(s) | Primary Module(s) |
|---|---|---|
| Auth/account | REQ-001..006 | `poi_broker/auth.py` |
| Alert browsing/filtering | REQ-007 | `poi_broker/app.py`, `services/search_service.py`, `services/filter_service.py` |
| CSV download | REQ-008 | `poi_broker/app.py` |
| Crossmatch query | REQ-009 | `poi_broker/app.py` |
| Visual query / watchlists | REQ-010, REQ-011 | `routes/visual_query.py`, `services/query_service.py`, `querybuilder_translator.py` |
| Lightcurve | REQ-012 | `routes/lightcurve.py`, `services/plotting_service.py` |
| Features | REQ-013 | `routes/features.py`, `services/feature_service.py` |
| Favorites/groups | REQ-014 | `routes/favorites.py`, `services/favorites_service.py` |
| Filter bookmarks | REQ-015 | `routes/filter_bookmarks.py`, `services/filter_service.py` |
| Observatories | REQ-016 | `routes/user_observatories.py` |
| Observing tool | REQ-017 | `observing_tool.py` |
| Classification | REQ-018 | `classification.py` |
| Settings | REQ-019 | `user_settings.py` |

---

## 8. Code-to-Spec Traceability Matrix

| Spec Requirement | Code (route/function) | Test |
|---|---|---|
| REQ-001..006 | `auth.py` (login/signup/verify/logout/forgot/reset/change) | `tests/test_auth.py` |
| REQ-007 | `app.py` `start()` + `SearchService`/`FilterService` | `tests/test_app.py`, `test_alert_id_filter.py`, `test_dec_filter.py`, `test_ant_passband_filter.py` |
| REQ-008 | `app.py` `/download_alerts_csv` | `tests/test_app.py` |
| REQ-009 | `app.py` `/query_crossmatches` | `tests/test_app.py` |
| REQ-010, 011 | `routes/visual_query.py` + `QueryService` | `tests/test_visual_query*.py` |
| REQ-012 | `routes/lightcurve.py` + `PlottingService` | `tests/test_lightcurve*.py` |
| REQ-013 | `routes/features.py` + `FeatureService` | `tests/test_features*.py` |
| REQ-014 | `routes/favorites.py` + `FavoritesService` | `tests/test_favorites.py` |
| REQ-015 | `routes/filter_bookmarks.py` + `FilterService` | `tests/test_filter_bookmarks.py` |
| REQ-016 | `routes/user_observatories.py` | `tests/test_user_observatories*.py` |
| REQ-017 | `observing_tool.py` | `tests/test_observing_tool*.py` |
| REQ-018 | `classification.py` | `tests/test_classification.py` |
| REQ-019 | `user_settings.py` | `tests/test_user_settings*.py` |
| NFR-001..010 | `settings.py`, `__init__.py` (limiter, CSRF, headers, pragmas) | `tests/test_security_regressions.py`, `tests/test_smoke_routes.py` |

---

## 9. Known Gaps / Mismatches

| # | Gap | Detail | Source |
|---|---|---|---|
| G1 | No formal external spec | All requirements inferred from code; no authoritative external document | — |
| G2 | No CI workflow | No GitHub Actions / pre-commit config; validation is manual | repo root |
| G3 | No Flask-Migrate | Schema migrations are manual SQL scripts only (`tools/*.sql`) | `tools/` |
| G4 | `requirements2.txt` not guaranteed installable | Reference freeze only; may not install on all systems | `requirements2.txt` |
| G5 | Coverage targets not enforced by CI | `tests/coverage.md` states targets but no gate | `tests/coverage.md` |
| G6 | Internal service contracts undocumented | `QueryService`, `FavoritesService`, etc. excluded (per scope decision) | `services/` |

---

## 10. Assumptions and Ambiguities

| # | Item | Status | Note |
|---|---|---|---|
| A1 | All functional requirements are inferred from code, not an authoritative spec | **Inferred** | Marked in §3 |
| A2 | `role`/`role_required` behavior beyond `user`/`admin` | **Uncertain** | Only `user` default and `is_admin` observed |
| A3 | Exact error message bodies for auth failures | **Uncertain** | Behavior tested, exact strings not enumerated here |
| A4 | Watchlist `rules_json` schema | **Uncertain** | Defined by `querybuilder_translator.py`; not enumerated in this doc |
| A5 | Observing tool custom-observatory auth boundary | **Uncertain** | 401 without auth; exact UX not specified |

---

## 11. Recommended Implementation Order

1. **Auth & account** (REQ-001..006) — foundational, blocks user-scoped features.
2. **Alert browsing/filtering + CSV + crossmatch** (REQ-007..009) — core read path.
3. **Favorites/groups** (REQ-014) — first user-scoped data feature.
4. **Filter bookmarks** (REQ-015) — depends on filter infrastructure.
5. **Visual query + watchlists** (REQ-010, 011) — depends on query builder.
6. **Lightcurve / features / classification** (REQ-012, 013, 018) — plotting layer.
7. **Observatories + observing tool** (REQ-016, 017).
8. **Settings** (REQ-019) — depends on observatories + feature columns.

---

## 12. Validation Strategy (tied to spec)

1. **Per-requirement tests:** Each REQ-xxx maps to an existing test file (§8). New features must add a matching `pytest` test.
2. **Full suite:** `python -m pytest -q` (or `python -m pytest -n auto -q --tb=no` for parallel).
3. **Route registration sanity:** `python -m flask --app wsgi:app routes` with `SECRET_KEY` set.
4. **Coverage gate:** aim 50–75%/module, 75%+ total per `tests/coverage.md`.
5. **Security regressions:** `tests/test_security_regressions.py` (CSRF/auth/rate-limit) must pass on any auth/security change.
6. **Manual smoke:** browser/API sanity check for any UI/endpoint change.

---

## 13. Top 5 Spec Gaps

1. **No authoritative external spec** — requirements are code-inferred; this document is the first formal reference.
2. **No CI enforcement** — no GitHub Actions; tests/coverage are not gated automatically.
3. **No schema migration tooling** — Flask-Migrate absent; schema changes are manual SQL.
4. **Coverage targets not enforced** — `tests/coverage.md` goals have no automated gate.
5. **Internal service contracts undocumented** — service-layer APIs excluded from this spec (scope decision); may need a follow-up doc.
