"""
Test the ztf_object_id filter logic in app.py line 162
This only tests:
1. Complete valid ztf_object_id (e.g. `ztf09abbc`) is set - no warning added
2. Empty or missing ztf_object_id - no warning added
"""
from flask import request
import pytest
from poi_broker.models import Ztf

def test_ztf_object_id_filter_complete_valid(app):
    """Test when complete valid ztf_object_id is provided - no warning."""
    client = app.test_client()
    response = client.get("/?ztf_object_id=ztf09abbc")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_ztf_object_id_filter_partial_invalid(app):
    """Test when partial/invalid ztf_object_id is provided - nothing will be returned but no warning message should be added."""
    client = app.test_client()
    response = client.get("/?ztf_object_id=12345")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_ztf_object_id_filter_empty_string(app):
    """Test when empty ztf_object_id string is provided - no warning message should be added."""
    client = app.test_client()
    response = client.get("/?ztf_object_id=")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_ztf_object_id_filter_missing(app):
    """Test when ztf_object_id parameter is missing - no warning message should be added."""
    client = app.test_client()
    response = client.get("/")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data
