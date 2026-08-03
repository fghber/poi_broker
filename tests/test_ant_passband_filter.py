"""
This only tests:
1. Complete valid ant_passband (e.g., `g`, `R`, or `i`) is set - no warning added
2. Empty or missing ant_passband - no warning added
3. Invalid ant_passband (e.g., 'xyz') - warning is added
"""
from flask import request
import pytest
from poi_broker.models import Ztf

def test_ant_passband_filter_complete_valid(app):
    """Test when complete valid ant_passband is provided - no warning."""
    client = app.test_client()
    response = client.get("/?ant_passband=g")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_ant_passband_filter_invalid_value(app):
    """Test that an invalid ant_passband (e.g. 'xyz') emits a warning."""
    client = app.test_client()
    response = client.get("/?ant_passband=xyz")
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    # Optional: also assert the specific message text
    assert b'Passband filter cannot be applied' in response.data

def test_ant_passband_filter_empty_string(app):
    """Test when empty ant_passband string is provided - no warning added."""
    client = app.test_client()
    response = client.get("/?ant_passband=")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_ant_passband_filter_missing(app):
    """Test when ant_passband parameter is missing - no warning added."""
    client = app.test_client()
    response = client.get("/")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

