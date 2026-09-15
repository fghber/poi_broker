## [json-name-type] watchlist/bookmarks/favorite-groups — non-string name caused 500
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_visual_query.py::test_watchlist_rejects_non_string_name, tests/test_filter_bookmarks.py::test_filter_bookmarks_reject_non_string_name, tests/test_favorites.py::test_favorite_group_rejects_non_string_name

## [json-locus-group-type] favorites POST — non-string locusId / bool groupId caused 500 or type-confused 404
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_favorites.py::test_favorite_rejects_non_string_locus_id, tests/test_favorites.py::test_favorite_rejects_boolean_group_id

## [whitespace-password] reset/change-password — 8 spaces passed len() unlike signup
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_security_regressions.py::test_reset_password_rejects_whitespace_only, tests/test_security_regressions.py::test_change_password_rejects_whitespace_only

## [serialize-invalid-utf8] helpers.serialize_fallback — undecodable bytes raised UnicodeDecodeError
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_helpers.py::test_safe_serialize_replaces_invalid_utf8_bytes

## [querybuilder-non-column] Filter._make_query — relationship attrs (e.g. classification) raised NotImplementedError
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_querybuilder_translator.py::test_querybuilder_rejects_non_column_attribute, tests/test_visual_query.py::test_preview_query_rejects_non_column_field

## [export-csrf-json] CSRFError handler — POST /export (JSON, not /api/) redirected instead of JSON 400
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_security_regressions.py::test_export_csrf_failure_returns_json

## [polar-empty-night] observing_tool — np.max on empty night_moon_alt during polar summer → 500
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_observing_tool.py::test_query_observing_plot_polar_summer_returns_moon_message
