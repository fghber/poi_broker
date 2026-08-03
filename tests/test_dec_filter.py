"""
Test the Dec filter logic in app.py.

Declination (Dec) is an angle in degrees and is only meaningful in the range
-90° to +90° (the real range is -30° to 90°). The filter must reject out-of-range values and
malformed input, and apply the filter for valid input.
"""
import pytest


def _get_dec(client, value):
    return client.get(f"/?locus_dec={value}")


def _assert_valid(response):
    assert b'<div class="alert alert-warning" role="alert">' not in response.data


def _assert_invalid(response):
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data


# --- Valid input: filter is applied, no warning ---

@pytest.mark.parametrize("value,description", [
    ("", "empty string"),
    ("20", "single positive value"),
    ("-10.5", "single negative value"),
    ("-90 -45", "negative range"),
    ("18.8 19.4", "positive range"),
    ("-1.2 -3.4", "negative decimal range"),
    ("<-45", "less than operator"),
    (">-89", "greater than operator"),
    (">-60 <-30", "both operators normalized"),
    ("-.56789", "leading decimal point"),
    ("0", "zero"),
    ("0.0", "float zero"),
    (" -20.02131 ", "whitespace trimmed"),
    ("-20   -15", "multiple spaces"),
    ("-90 90", "maximum valid range"),
    ("-12.345678901234567890", "many decimal places (truncated)"),
    ("<-20 -15", "less than with range"),
    (">-20 <-15", "greater than less than range"),
    ("<-15 -20", "reversed range with signs"),
])
def test_dec_filter_valid_inputs(app, value, description):
    """Valid Dec inputs should not produce a warning."""
    _assert_valid(_get_dec(app.test_client(), value))


# --- Invalid input: warning is shown ---

@pytest.mark.parametrize("value,description", [
    ("invalid-dec", "non-numeric string"),
    ("95", "out of range positive"),
    ("-95", "out of range negative"),
    ("abc123", "letters with numbers at end"),
    ("123abc", "numbers with letters at end"),
    ("-20.02131%", "special character"),
    ("-20.02131ABC", "mixed case letters"),
    ("NaN", "nan-like string"),
    ("-20,-15", "comma separated"),
    ("-20.02131,", "trailing comma"),
    (",-20.02131", "leading comma"),
    ("-12345678901234567890 -1", "very long range start"),
    ("-1 12345678901234567890", "very long range end"),
    ("-20.02131α", "unicode character"),
])
def test_dec_filter_invalid_inputs(app, value, description):
    """Invalid Dec inputs should produce a warning."""
    _assert_invalid(_get_dec(app.test_client(), value))