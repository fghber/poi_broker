"""
Test the locus_id filter logic in app.py line 171
This only tests:
1. Complete locus_id (e.g. `ANT2020a4izc`) is set filter is applied and no warning is added
2. Empty or missing locus_id - no warning added
"""
from flask import request
import pytest
from poi_broker.models import Ztf

def test_locus_id_filter_complete_valid(app):
    """Test when complete valid locus_id is provided - no warning."""
    client = app.test_client()
    response = client.get("/?locus_id=ANT2020a4izc")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_locus_id_filter_partial_invalid(app):
    """Test when partial/invalid locus_id is provided - nothing will be return but no warning message should be added."""
    client = app.test_client()
    response = client.get("/?locus_id=12345")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_locus_id_filter_empty_string(app):
    """Test when empty locus_id string is provided - no warning message should be added."""
    client = app.test_client()
    response = client.get("/?locus_id=")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_locus_id_filter_missing(app):
    """Test when locus_id parameter is missing - no warning message should be added."""
    client = app.test_client()
    response = client.get("/")
    assert b'<div class="alert alert-warning" role="alert">' not in response.data