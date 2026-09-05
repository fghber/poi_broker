"""
Unit tests for FilterService / FloatFilter / IntFilter.

These tests exercise the filter helpers directly against real SQLAlchemy
columns on a seeded Ztf model (via the `app` fixture from conftest.py),
so we verify the *actual* SQL produced — not mocks of Query.filter.

Coverage:
  - FilterService.apply_filter: single exact, single with offset, single
    with >, single with <, range, range sort.
  - FloatFilter: decimals=0 exact-equal, decimals=N offset, > and <.
  - IntFilter: exact, range, >.

MJD (date -> MJD) coverage is consolidated in tests/test_mjd_filter.py.
Route-level filter coverage (warning-vs-no-warning, parser integration)
lives in tests/test_ra_filter.py, tests/test_dec_filter.py,
tests/test_mjd_date_filter.py, tests/test_date_filter_logic.py, and
tests/test_smoke_routes.py.
"""
import pytest
from astropy.time import Time

from poi_broker import db
from poi_broker.models import Ztf
from poi_broker.services.filter_service import (
    FilterService,
    FloatFilter,
    IntFilter,
)
from poi_broker.services.input_parser import ParsedValue, InputParser
from poi_broker.services.search_service import SearchService


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def seeded_ztf(app):
    """Seed three Ztf rows with known, spaced-out numeric values.

    Layout:
        Z1: ra=10.0,  dec=20.0,  mag=12.0,  num_alerts=1
        Z2: ra=50.0,  dec=0.0,   mag=18.0,  num_alerts=5
        Z3: ra=180.0, dec=-45.0, mag=20.5,  num_alerts=10
        plus three different MJD values (date-only branch).
    """
    with app.app_context():
        db.session.add_all([
            Ztf(
                alert_id="ztf_candidate:1111111111111111111",
                date_alert_mjd=Time(59945.0, format="mjd").mjd,  # 2023-01-01
                locus_id="L1", ztf_object_id="Z1",
                locus_ra=10.0, locus_dec=20.0, ant_mag_corrected=12.0,
                num_alerts=1, num_mag_values=10,
            ),
            Ztf(
                alert_id="ztf_candidate:2222222222222222222",
                date_alert_mjd=Time(60110.5, format="mjd").mjd,  # 2023-06-15 12:00
                locus_id="L2", ztf_object_id="Z2",
                locus_ra=50.0, locus_dec=0.0, ant_mag_corrected=18.0,
                num_alerts=5, num_mag_values=50,
            ),
            Ztf(
                alert_id="ztf_candidate:3333333333333333333",
                date_alert_mjd=Time(60308.0, format="mjd").mjd,  # 2023-12-31
                locus_id="L3", ztf_object_id="Z3",
                locus_ra=180.0, locus_dec=-45.0, ant_mag_corrected=20.5,
                num_alerts=10, num_mag_values=100,
            ),
        ])
        db.session.commit()
        yield


def _pv(value: str, operator: str = "") -> ParsedValue:
    """Build a ParsedValue with the given operator and clean value."""
    return ParsedValue(raw=f"{operator}{value}", operator=operator, value=value)


def _locus_ids(query) -> list[str]:
    """Return locus_id values from a Ztf query, sorted for stable comparison."""
    return sorted(r.locus_id for r in query.all())


# ---------------------------------------------------------------------------
# FilterService.apply_filter
# ---------------------------------------------------------------------------

class TestFilterServiceApplyFilter:
    def test_single_value_no_offset_uses_equality(self, app, seeded_ztf):
        """Single value, no offset → `==` (single .filter call)."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("50.0")],
                convert_callback=float, upper_offset=0.0, lower_offset=0.0,
            )
            assert _locus_ids(result) == ["L2"]

    def test_single_value_with_offset_uses_range(self, app, seeded_ztf):
        """Single value, non-zero offset → `>=` AND `<=` (two .filter calls)."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            # 50.0 ± 0.5 → only L2 (10 is too low, 180 is too high)
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("50.0")],
                convert_callback=float, upper_offset=0.5, lower_offset=0.5,
            )
            assert _locus_ids(result) == ["L2"]

    def test_greater_than_subtracts_lower_offset(self, app, seeded_ztf):
        """> uses `>= converted - lower_offset` only."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            # > 100 with lower_offset=0.1 → >= 99.9 → L2 (50) is excluded, L3 (180) included
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("100", operator=">")],
                convert_callback=float, upper_offset=0.0, lower_offset=0.1,
            )
            assert _locus_ids(result) == ["L3"]

    def test_less_than_adds_upper_offset(self, app, seeded_ztf):
        """< uses `<= converted + upper_offset` only."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            # < 50 with upper_offset=0.1 → <= 50.1 → L1 (10) included, L2 (50) included, L3 (180) excluded
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("50", operator="<")],
                convert_callback=float, upper_offset=0.1, lower_offset=0.0,
            )
            assert _locus_ids(result) == ["L1", "L2"]

    def test_range_sorts_values_regardless_of_input_order(self, app, seeded_ztf):
        """Range filter sorts the two endpoints and uses both bounds."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            # Range given in reverse order: "180 10" → 10..180 → all three
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("180"), _pv("10")],
                convert_callback=float, upper_offset=0.0, lower_offset=0.0,
            )
            assert _locus_ids(result) == ["L1", "L2", "L3"]

    def test_range_with_offset_uses_offset(self, app, seeded_ztf):
        """Range filter applies upper/lower offset to the bounds."""
        with app.app_context():
            service = FilterService()
            q = db.session.query(Ztf)
            # Range 10..50 with offset 0.5 → 9.5..50.5 → L1 (10) and L2 (50) only
            result = service.apply_filter(
                q, Ztf.locus_ra, [_pv("10"), _pv("50")],
                convert_callback=float, upper_offset=0.5, lower_offset=0.5,
            )
            assert _locus_ids(result) == ["L1", "L2"]


# ---------------------------------------------------------------------------
# FloatFilter (the previously-mutated, now-stateless class)
# ---------------------------------------------------------------------------

class TestFloatFilter:
    def test_decimals_zero_uses_equality(self, app, seeded_ztf):
        """decimals=0 → no offset → equality match against the column."""
        with app.app_context():
            ff = FloatFilter()
            q = db.session.query(Ztf)
            result = ff(q, Ztf.locus_ra, [_pv("50.0")], decimals=0)
            assert _locus_ids(result) == ["L2"]

    def test_decimals_positive_applies_offset(self, app, seeded_ztf):
        """decimals=2 → offset 0.005 → matches values within ±0.005."""
        with app.app_context():
            ff = FloatFilter()
            q = db.session.query(Ztf)
            # 50.0 ± 0.005 → L2 only
            result = ff(q, Ztf.locus_ra, [_pv("50.00")], decimals=2)
            assert _locus_ids(result) == ["L2"]

    def test_decimals_value_outside_offset_misses(self, app, seeded_ztf):
        """Value outside the decimals-derived offset must not match."""
        with app.app_context():
            ff = FloatFilter()
            q = db.session.query(Ztf)
            # decimals=5 → offset = 0.5e-5 = 5e-6. 50.01 is 0.01 away from 50.0 → miss.
            result = ff(q, Ztf.locus_ra, [_pv("50.01")], decimals=5)
            assert _locus_ids(result) == []

    def test_decimals_value_inside_offset_hits(self, app, seeded_ztf):
        """Value inside the decimals-derived offset must match."""
        with app.app_context():
            ff = FloatFilter()
            q = db.session.query(Ztf)
            # decimals=5 → offset = 5e-6. 50.000005 is 5e-6 away from 50.0 → hit.
            result = ff(q, Ztf.locus_ra, [_pv("50.000005")], decimals=5)
            assert _locus_ids(result) == ["L2"]

    def test_instance_is_reusable_across_decimal_values(self, app, seeded_ztf):
        """Regression test: a single FloatFilter instance must accept different
        `decimals` values across calls without leaking state between them.

        Before the stateless refactor, FloatFilter stored `decimals` on
        `self`, which made shared instances unsafe under threading.
        """
        with app.app_context():
            ff = FloatFilter()
            # First call with decimals=5, then decimals=0 — both should
            # be honoured independently.
            r1 = ff(db.session.query(Ztf), Ztf.locus_ra, [_pv("50.0")], decimals=5)
            r2 = ff(db.session.query(Ztf), Ztf.locus_ra, [_pv("50.0")], decimals=0)
            assert _locus_ids(r1) == ["L2"]
            assert _locus_ids(r2) == ["L2"]

    def test_range_query(self, app, seeded_ztf):
        """FloatFilter on a two-value range returns the expected subset."""
        with app.app_context():
            ff = FloatFilter()
            q = db.session.query(Ztf)
            result = ff(q, Ztf.locus_ra, [_pv("0"), _pv("100")], decimals=5)
            assert _locus_ids(result) == ["L1", "L2"]


# ---------------------------------------------------------------------------
# IntFilter
# ---------------------------------------------------------------------------

class TestIntFilter:
    def test_exact_match(self, app, seeded_ztf):
        with app.app_context():
            inf = IntFilter()
            q = db.session.query(Ztf)
            result = inf(q, Ztf.num_alerts, [_pv("5")])
            assert _locus_ids(result) == ["L2"]

    def test_range_query(self, app, seeded_ztf):
        with app.app_context():
            inf = IntFilter()
            q = db.session.query(Ztf)
            result = inf(q, Ztf.num_alerts, [_pv("3"), _pv("9")])
            assert _locus_ids(result) == ["L2"]

    def test_greater_than(self, app, seeded_ztf):
        with app.app_context():
            inf = IntFilter()
            q = db.session.query(Ztf)
            # > 5 with no offset → >= 5 → both L2 (5) and L3 (10).
            result = inf(q, Ztf.num_alerts, [_pv("5", operator=">")])
            assert _locus_ids(result) == ["L2", "L3"]


# ---------------------------------------------------------------------------
# SearchService integration (parse -> apply filter end-to-end)
# ---------------------------------------------------------------------------

class TestSearchServiceFilters:
    def test_float_filter_creates_equality_filter_for_exact_value(self, app):
        with app.app_context():
            query = db.session.query(Ztf)
            parsed = InputParser().parse_numbers('1.23')
            filtered = SearchService().apply_float_filter(query, Ztf.date_alert_mjd, parsed)
            dialect = query.session.get_bind().dialect
            compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

        assert 'date_alert_mjd' in compiled
        assert '= 1.23' in compiled or '== 1.23' in compiled

    def test_int_filter_handles_two_values_and_orders_inputs(self, app):
        with app.app_context():
            query = db.session.query(Ztf)
            parsed = InputParser().parse_numbers('5 1')
            filtered = SearchService().apply_int_filter(query, Ztf.date_alert_mjd, parsed)
            dialect = query.session.get_bind().dialect
            compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

        assert 'date_alert_mjd >= 1' in compiled
        assert 'date_alert_mjd <= 5' in compiled

    def test_mjd_filter_converts_iso_date_to_mjd_range(self, app):
        with app.app_context():
            query = db.session.query(Ztf)
            parsed = InputParser().parse_dates('2025-01-15')
            filtered = SearchService().apply_mjd_filter(query, Ztf.date_alert_mjd, parsed)
            dialect = query.session.get_bind().dialect
            compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

        assert 'date_alert_mjd' in compiled
        assert '>=' in compiled and '<=' in compiled
