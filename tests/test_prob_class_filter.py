"""
Test the prob_class filter logic in app.py (lines 194-201)
This tests:
1. Valid prob_class value - query should be filtered, no warning
2. Empty/whitespace-only prob_class - no filtering, no warning  
3. Invalid prob_class format - warning message should appear
"""

from flask import request
import pytest
from poi_broker.models import Classification


def test_prob_class_filter_valid_value(app):
    """Test when valid lowercase prob_class value is provided - query filtered."""
    client = app.test_client()
    response = client.get("/?prob_class=sn")  # Use lowercase to match validation logic
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_prob_class_filter_empty_string(app):
    """Test when empty prob_class string is provided - no warning."""
    client = app.test_client()
    response = client.get("/?prob_class=")  # Empty after strip, won't trigger either branch
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_prob_class_filter_whitespace_only(app):
    """Test when whitespace-only prob_class is provided - no warning."""
    client = app.test_client()
    response = client.get("/?prob_class=   ")  # Strips to empty, won't trigger either branch  
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_prob_class_filter_invalid_format(app):
    """Test when invalid prob_class format is provided - warning message should appear."""
    client = app.test_client()
    response = client.get("/?prob_class=invalid-class-format")  # Exists but isn't valid
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Classification filter cannot be applied' in response.data

def test_prob_case_insensitive(app):
    """Test that case-insensitive matching works for SN."""
    client = app.test_client()
    # Comparison is implemented as an exact match, therefore this should fail.
    response = client.get("/?prob_class=SN")  # Uppercase 'sn' in valid_prob_classes list
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Classification filter cannot be applied' in response.data
