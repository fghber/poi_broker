"""
Test the MJD date filter logic in app.py (lines 145-149)
This tests the happy path where date_input is provided and the filter is applied,
and the else path where no date is provided or invalid input.
"""
import pytest
from poi_broker import models
from poi_broker.models import Ztf


def test_mjd_date_filter_happy_path_with_decimal(app):
    """Test when valid MJD decimal is provided - no warning (filter should be applied)."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=59190.12")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_range(app):
    """Test when MJD range is provided - no warning (filter should be applied)."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=59190 59191")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_upper_bound(app):
    """Test when MJD upper bound is provided - no warning (filter should be applied)."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=>59190")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_lower_bound(app):
    """Test when MJD lower bound is provided - no warning (filter should be applied)."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=<59190")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_no_date_provided(app):
    """Test when no date is provided - no warning message should be added."""
    client = app.test_client()
    response = client.get("/")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_invalid_date_format(app):
    """Test when invalid MJD format is provided - warning message should be added."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=invalid-date")
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'MJD filter cannot be applied' in response.data


def test_mjd_date_filter_non_numeric(app):
    """Test when non-numeric MJD is provided - warning message should be added."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=abc123")
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'MJD filter cannot be applied' in response.data


def test_mjd_date_filter_empty_string(app):
    """Test when empty string is provided - no warning (treated as no filter)."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_negative(app):
    """Test when negative MJD is provided - should be handled gracefully."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=-59190.12")
    # MJD can be negative, so this should not trigger a warning
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def test_mjd_date_filter_with_time_component(app):
    """Test when MJD with time component is provided - should be handled."""
    client = app.test_client()
    response = client.get("/?date_alert_mjd=59190.12345")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data