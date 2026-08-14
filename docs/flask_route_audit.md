# Flask route handler audit (archived)

Reviewed 2026-08-13 · patched 2026-08-14 · remaining nits discarded 2026-08-14.

Scope: stability, validation, decoupling, request-context safety, and HTTP consistency across `favorites`, `features`, `filter_bookmarks`, `lightcurve`, and `user_observatories` (plus `favorites_service`, `feature_service`, `plotting_service`).

This file archives the Cursor canvas findings. The live canvas was removed after the audit set was closed.

## Outcome

| Bucket | Count |
|---|---|
| Patched this cycle | 8 |
| Already in codebase | 1 (bookmark `UNIQUE(user_id, name)`) |
| Closed / won't do | 11 |
| Still open | 0 |

## Patched

| Endpoint | Change |
|---|---|
| `POST /api/favorite` | Require boolean `fav`; check `groupId` ownership on insert/update; cap `locusId` at 128 |
| `GET /api/favorite` | 400 on missing/empty/overlong `locusId` |
| `GET /api/favorites` | List via `get_user_favorites()`; 400 on invalid `groupId`; `groupId=null` still ungrouped |
| `PATCH /api/favorite/<id>/group` | Body must be an object; `groupId` must be int or null (bool rejected) |
| `POST /api/favorite-groups` | Name capped at 128 after strip |
| `DELETE /api/favorite-groups/<id>` | Orphan `UPDATE` also filters `user_id` |
| Plot/CSV GETs | `READ_RATE_LIMIT_LAX` on `/query_lightcurve_data`, `/locus_plot_csv`, `/query_featureplot_data` |
| `GET /query_features` | Dropped unused `safe_serialize` / related imports (jsonify is fine for float/`None` columns) |

## Won't do (intentional)

| Item | Reason |
|---|---|
| `GET /api/favorite-groups` → 200 `[]` on exception | Contract; `profile.html` ignores status; 500 needs a UI error path |
| Bookmarks DELETE 404 vs observatories 200 `already_deleted` | Both intentional; harmonizing breaks one client |
| Observatories CRUD in the route | Same pattern as bookmarks; service extract is hygiene |
| Observatories duplicate create 200 `already_exists` | Deliberate idempotency vs groups/bookmarks 409 |
| `EarthLocation.get_site_names()` on selected-custom delete | Rare write; reuse `_get_builtin_observatory_options` later if hot |
| URL layout (`/api/*` vs `/query_*`) | Renaming breaks `main.html` |
| Silent feature-plot unknown-feature fallback | Avoids blank plots for stale saved feature lists |
| `csv.writer` on `/locus_plot_csv` | Hygiene; values from `Ztf`, not request input |
| Global JSON 500/429 handlers | Deferred; plot routes return HTML fragments |
| `current_user` decoupling in services | Callers are all `@login_required` |
| TimezoneFinder singleton on create | Rare authenticated write |

## Residual (not defects)

- `POST /api/favorite` still accepts `groupId: true` via `isinstance(..., int)` because `bool` subclasses `int`. PATCH rejects bool. First-party UI never sends it.
- `get_available_features()` remains in `feature_service` but is unused by routes.

## What stayed solid

- `/api/*` mutating routes behind `login_required`; JSON 401 for unauthenticated `/api/` paths
- Flask-WTF CSRF enforced (`tests/test_security_regressions.py`)
- Validation patterns to copy: `filter_bookmarks._sanitize_params`, `user_observatories._validate_payload`
- Observatory delete: selection fallback mutates the same session and commits once with the row delete
