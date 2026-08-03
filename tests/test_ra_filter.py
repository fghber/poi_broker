"""
Test the RA filter logic in app.py.

Right ascension is an angle in degrees and is only meaningful in the range
0° to 360°. The filter must reject out-of-range values (negatives, >360) and
malformed input, and apply the filter for valid input.
"""
import pytest


def _get_ra(client, value):
    return client.get(f"/?locus_ra={value}")


def _assert_valid(response):
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def _assert_invalid(response):
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Ra filter cannot be applied' in response.data


# --- Valid input: filter is applied, no warning ---

@pytest.mark.parametrize("value,description", [
    ("", "empty string"),
    ("118.61421", "single value"),
    ("80 90", "range"),
    ("0", "boundary zero"),
    ("360", "boundary 360"),
    ("0.001", "just inside lower"),
    ("359.999", "just inside upper"),
    ("0.0", "float zero"),
    (".56789", "leading decimal"),
    (" 123.456 ", "whitespace trimmed"),
    ("10   20", "multiple spaces in range"),
    ("1.2 3.4", "decimal range"),
    ("<80", "less than operator"),
    (">80", "greater than operator"),
    ("<80 90", "range with operators"),
    (">20 <10", "contradictory operators normalized"),
])
def test_ra_filter_valid_inputs(app, value, description):
    """Valid RA inputs should not produce a warning."""
    _assert_valid(_get_ra(app.test_client(), value))


# --- Invalid input: warning is shown ---

@pytest.mark.parametrize("value,description", [
    ("invalid-ra", "non-numeric string"),
    ("123.456ABC", "mixed letters"),
    ("123.456%", "special character"),
    ("123.456α", "unicode"),
    ("123.456\\", "backslash"),
    ("nan", "nan string"),
    ("1,2", "comma separated"),
    ("123.456,", "trailing comma"),
    (",123.456", "leading comma"),
    ("1,2,3", "multiple commas"),
    ("1e5", "exponent notation"),
    ("%2B-10%2020.5", "plus sign (URL encoded)"),
])
def test_ra_filter_invalid_format(app, value, description):
    """Invalid format inputs should produce a warning."""
    _assert_invalid(_get_ra(app.test_client(), value))


# --- Out-of-range (RA must be 0..360) ---

@pytest.mark.parametrize("value,description", [
    ("-10.5", "negative value"),
    ("-10 20", "negative range"),
    ("370", "over 360"),
    ("360.1", "just over upper"),
    ("-0.1", "just under lower"),
    ("1234567890.1", "large out of range"),
])
def test_ra_filter_out_of_range(app, value, description):
    """Out-of-range RA values should produce a warning."""
    _assert_invalid(_get_ra(app.test_client(), value))
