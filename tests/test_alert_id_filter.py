"""
Test the alert_id filter logic in app.py (lines 152-159)
This tests:
1. Complete alert_id with 18+ digits (direct query)
2. Partial alert_id with catalog prefix only (prefix matching)
3. Invalid alert_id (warning message should be shown)
"""
from flask import request
import pytest
from poi_broker.models import Ztf

def test_alert_id_filter_complete_id_ztf_candidate(app):
    """Test when complete alert_id with 'ztf_candidate' prefix is provided - no warning."""
    client = app.test_client()
    response = client.get("/?alert_id=ztf_candidate:3432493322215010047")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_alert_id_filter_complete_id_lsst(app):
    """Test when complete alert_id with 'lsst' prefix is provided - no warning."""
    client = app.test_client()
    response = client.get("/?alert_id=lsst:170094456539709554")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_alert_id_filter_prefix_only_ztf(app):
    """Test when partial alert_id with 'ztf' prefix only - no warning (prefix matching)."""
    client = app.test_client()
    response = client.get("/?alert_id=ztf")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_alert_id_filter_prefix_only_lsst(app):
    """Test when partial alert_id with 'lsst' prefix only - no warning (prefix matching)."""
    client = app.test_client()
    response = client.get("/?alert_id=lsst")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_alert_id_filter_no_alert_id_provided(app):
    """Test when no alert_id is provided - no warning message should be added."""
    client = app.test_client()
    response = client.get("/")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_alert_id_filter_partial_alert_id(app):
    """Test when partial alert_id is provided - warning message should be added."""
    client = app.test_client()
    response = client.get("/?alert_id=ztf_candidate:335155568501")
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Alert ID cannot be filter by partial IDs' in response.data

def test_alert_id_filter_invalid_format(app):
    """Test when invalid alert_id format is provided - warning message should be added."""
    client = app.test_client()
    response = client.get("/?alert_id=invalid-alert-id")
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Alert ID cannot be filter by partial IDs' in response.data

def test_alert_id_filter_empty_string(app):
    """Test when empty alert_id string is provided - no warning message should be added."""
    client = app.test_client()
    response = client.get("/?alert_id=")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data