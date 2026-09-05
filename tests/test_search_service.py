"""
Unit tests for SearchService.

Coverage targets the thin delegation methods (parse_float_input,
parse_int_input, parse_date_input, parse_mjd_input) and the guard-clause
early returns in apply_float_filter / apply_int_filter / apply_mjd_filter
that were previously uncovered.

The underlying parser/filter logic is already covered in
test_filter_service.py and test_mjd_filter.py; these tests only lock down
the delegation contracts and the empty/None ParseResult short-circuits.
"""
import pytest

from poi_broker import db
from poi_broker.models import Ztf
from poi_broker.services.input_parser import ParseResult
from poi_broker.services.search_service import SearchService


# ---------------------------------------------------------------------------
# parse_float_input / parse_int_input: delegation + allowed_range forwarding
# ---------------------------------------------------------------------------

class TestParseNumberInputs:
    def test_parse_float_input_returns_parse_result_for_valid_value(self):
        parsed = SearchService().parse_float_input("1.23")
        assert parsed is not None
        assert len(parsed.values) == 1
        assert parsed.values[0].value == "1.23"

    def test_parse_float_input_returns_none_for_invalid_text(self):
        # Invalid characters fall through InputParser.parse_numbers -> None.
        parsed = SearchService().parse_float_input("abc")
        assert parsed is None

    def test_parse_float_input_forwards_allowed_range_rejecting_out_of_range(self):
        # allowed_range must be forwarded; an out-of-range value is rejected (None).
        parsed = SearchService().parse_float_input("999.0", allowed_range=(0.0, 10.0))
        assert parsed is None

    def test_parse_float_input_forwards_allowed_range_accepting_in_range(self):
        parsed = SearchService().parse_float_input("5.0", allowed_range=(0.0, 10.0))
        assert parsed is not None
        assert parsed.values[0].value == "5.0"

    def test_parse_int_input_returns_parse_result_for_valid_value(self):
        parsed = SearchService().parse_int_input("5")
        assert parsed is not None
        assert parsed.values[0].value == "5"

    def test_parse_int_input_returns_none_for_invalid_text(self):
        parsed = SearchService().parse_int_input("abc")
        assert parsed is None

    def test_parse_int_input_forwards_allowed_range_rejecting_out_of_range(self):
        parsed = SearchService().parse_int_input("100", allowed_range=(0.0, 10.0))
        assert parsed is None


# ---------------------------------------------------------------------------
# parse_date_input / parse_mjd_input: delegation (never return None)
# ---------------------------------------------------------------------------

class TestParseDateInputs:
    def test_parse_date_input_returns_empty_result_for_invalid_text(self):
        # parse_dates never returns None; invalid input -> empty ParseResult.
        parsed = SearchService().parse_date_input("not-a-date")
        assert isinstance(parsed, ParseResult)
        assert parsed.values == []

    def test_parse_date_input_parses_iso_date(self):
        parsed = SearchService().parse_date_input("2023-10-05")
        assert len(parsed.values) == 1
        assert parsed.values[0].value == "2023-10-05"

    def test_parse_mjd_input_delegates_to_date_parser(self):
        parsed = SearchService().parse_mjd_input("2023-10-05")
        assert len(parsed.values) == 1
        assert parsed.values[0].value == "2023-10-05"

    def test_parse_mjd_input_returns_empty_result_for_invalid_text(self):
        parsed = SearchService().parse_mjd_input("garbage")
        assert isinstance(parsed, ParseResult)
        assert parsed.values == []


# ---------------------------------------------------------------------------
# apply_*_filter: guard clause returns query unchanged on empty/None parsed
# ---------------------------------------------------------------------------

class TestApplyFilterGuardClauses:
    @pytest.fixture()
    def base_query(self, app):
        with app.app_context():
            yield db.session.query(Ztf)

    def test_apply_float_filter_returns_query_unchanged_for_empty_parsed(self, base_query):
        parsed = ParseResult(values=[])
        assert SearchService().apply_float_filter(base_query, Ztf.ant_mag_corrected, parsed) is base_query

    def test_apply_float_filter_returns_query_unchanged_for_none_parsed(self, base_query):
        assert SearchService().apply_float_filter(base_query, Ztf.ant_mag_corrected, None) is base_query

    def test_apply_int_filter_returns_query_unchanged_for_empty_parsed(self, base_query):
        parsed = ParseResult(values=[])
        assert SearchService().apply_int_filter(base_query, Ztf.num_alerts, parsed) is base_query

    def test_apply_int_filter_returns_query_unchanged_for_none_parsed(self, base_query):
        assert SearchService().apply_int_filter(base_query, Ztf.num_alerts, None) is base_query

    def test_apply_mjd_filter_returns_query_unchanged_for_empty_parsed(self, base_query):
        parsed = ParseResult(values=[])
        assert SearchService().apply_mjd_filter(base_query, Ztf.date_alert_mjd, parsed) is base_query

    def test_apply_mjd_filter_returns_query_unchanged_for_none_parsed(self, base_query):
        assert SearchService().apply_mjd_filter(base_query, Ztf.date_alert_mjd, None) is base_query
