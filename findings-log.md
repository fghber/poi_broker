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

## [iso-date-uncalendar] InputParser.parse_dates — ISO-looking tokens skipped calendar validation; to_mjd raised ErfaError
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_mjd_filter.py::TestMjdRoute::test_iso_date_with_invalid_month_shows_warning, tests/test_input_parser.py::TestInputParser::test_parse_dates_rejects_invalid_month

## [observing-uncalendar-date] observing_tool — obs_date split('-') without strptime; invalid calendar dates 500'd in Time()
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_observing_tool.py::test_query_observing_plot_rejects_non_calendar_date_with_numeric_coords, tests/test_observing_tool.py::test_query_observing_plot_rejects_impossible_calendar_date

## [querybuilder-malformed-rule] preview/export/watchlist — non-dict rules, non-string fields, unknown tables, unknown operators 500'd
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_visual_query.py::test_preview_query_rejects_unknown_table, tests/test_visual_query.py::test_preview_query_rejects_unsupported_operator, tests/test_visual_query.py::test_preview_query_rejects_non_dict_rule, tests/test_visual_query.py::test_preview_query_rejects_non_string_field

## [bool-latlon] POST /api/user-observatories — JSON true/false passed float() as 1.0/0.0
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_user_observatories.py::test_create_user_observatory_validation_errors (latitude True)

## [bool-custom-obs-id] last-selected observatory — isinstance(id, int) accepted JSON true as 1
- Date: 2026-09-15
- Status: fixed
- Repro: tests/test_user_observatories.py::test_save_last_observatory_api_rejects_boolean_custom_id, tests/test_observatory_context.py::TestResolveSelectedObservatory::test_custom_selection_boolean_id_returns_none
