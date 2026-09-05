"""
Unit tests for the Jinja template filters defined in poi_broker.app.

Covers:
- mag_filter: rounding and the 0.0 (falsy) edge case
- astro_filter: passband mapping
- format_mjd_readable: MJD -> human-readable UTC string
- epoch_utc_date: epoch-seconds -> UTC date string (User.password_changed_at)
"""

from datetime import datetime, timezone

from poi_broker.app import astro_filter, epoch_utc_date, format_mjd_readable, mag_filter


def test_mag_filter_rounds_to_three_decimals():
    assert mag_filter(12.34567) == 12.346
    assert mag_filter(12.345) == 12.345
    assert mag_filter(12.0) == 12.0


def test_mag_filter_handles_zero():
    # Regression: 0.0 is falsy, so the old `if num:` guard returned None.
    # A magnitude of 0.0 is valid and must be preserved.
    assert mag_filter(0.0) == 0.0


def test_mag_filter_handles_none():
    assert mag_filter(None) is None


def test_astro_filter_maps_passbands():
    assert astro_filter("g") == "g"
    assert astro_filter("R") == "R"
    assert astro_filter("i") == "i"


def test_astro_filter_unknown_passband_returns_empty():
    assert astro_filter("z") == ""
    assert astro_filter("") == ""


def test_format_mjd_readable_returns_utc_string():
    # MJD 58849.0 == 2020-01-01 00:00:00 UTC
    assert format_mjd_readable(58849.0) == "2020-01-01 00:00:00"


def test_format_mjd_readable_handles_none():
    assert format_mjd_readable(None) == ""


def test_format_mjd_readable_handles_invalid():
    assert format_mjd_readable("not-a-number") == ""


def test_epoch_utc_date_renders_utc_date():
    # 2026-08-27 12:00:00 UTC as epoch seconds
    epoch = int(datetime(2026, 8, 27, 12, 0, tzinfo=timezone.utc).timestamp())
    assert epoch_utc_date(epoch) == "2026-08-27"


def test_epoch_utc_date_uses_utc_not_local_tz():
    # 2026-08-27 23:30 UTC is still 2026-08-27 in UTC but already 2026-08-28
    # in UTC+1 — the filter must report the UTC date.
    epoch = int(datetime(2026, 8, 27, 23, 30, tzinfo=timezone.utc).timestamp())
    assert epoch_utc_date(epoch) == "2026-08-27"


def test_epoch_utc_date_handles_none():
    assert epoch_utc_date(None) == ""


def test_epoch_utc_date_handles_invalid():
    assert epoch_utc_date("not-a-number") == ""
    assert epoch_utc_date(float("nan")) == ""
