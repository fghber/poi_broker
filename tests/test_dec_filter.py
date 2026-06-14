"""
Test the Dec filter logic in app.py (lines 180-185)
This tests the happy path where dec_input is provided and the filter is applied,
and the else path where no valid Dec input is provided.
"""
import pytest
from poi_broker import app

def test_dec_filter_logic_empty_string(app):
    """Test when dec_input is an empty string - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=")
    # Warning message SHOULD appear for empty input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_no_input_provided(app):
    """Test when no locus_dec is provided - filter_warning_message should be empty."""
    client = app.test_client()
    response = client.get("/")
    # No warning message for missing input (this path doesn't add to the string)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_happy_path_with_single_value(app):
    """Test when dec_input is provided as a single valid number - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=20")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_negative_value(app):
    """Test when dec_input is provided as a negative number - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-10.5")
    # Filter warning message should NOT appear for valid input (negative dec)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_happy_path_with_negative_range(app):
    """Test when dec_input is provided as a range - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-90 -45")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_happy_path_with_positive_range(app):
    """Test when dec_input is provided as a range - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=18.8 19.4")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_negative_decimal_numbers_in_range(app):
    """Test when dec_input has dots between range numbers - No filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-1.2 -3.4")
    # Warning message SHOULD NOT appear for valid input (decimal numbers)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_happy_path_with_less_than_sign(app):
    """Test when dec_input is provided as range with less than sign - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=<-45")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_happy_path_with_greater_than_sign(app):
    """Test when dec_input is provided as range with greater than sign - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=>-89")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_both_signs(app):
    """Test when dec_input is provided with both > and < signs - signs ignored, treated as range."""
    client = app.test_client()
    response = client.get("/?locus_dec=>-60 <-30")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_invalid_number_format(app):
    """Test when invalid number format is provided - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=invalid-dec")
    # Warning message SHOULD appear for invalid input
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_out_of_range_value(app):
    """Test when dec_input is out of valid range (-90 to +90) - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=95")
    # Warning message SHOULD appear for out of range input
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_out_of_range_negative_value(app):
    """Test when dec_input is out of valid range (-90 to +90) - negative out of range."""
    client = app.test_client()
    response = client.get("/?locus_dec=-95")
    # Warning message SHOULD appear for out of range input
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_boundary_positive(app):
    """Test when dec_input is exactly at positive boundary (+90) - should be valid."""
    client = app.test_client()
    response = client.get("/?locus_dec=90")
    # Filter warning message should NOT appear for boundary value
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_boundary_negative(app):
    """Test when dec_input is exactly at negative boundary (-90) - should be valid."""
    client = app.test_client()
    response = client.get("/?locus_dec=-90")
    # Filter warning message should NOT appear for boundary value
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_boundary_just_over(app):
    """Test when dec_input is just over positive boundary (+90.1) - should trigger warning."""
    client = app.test_client()
    response = client.get("/?locus_dec=90.1")
    # Warning message SHOULD appear for value just over boundary
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_boundary_just_under(app):
    """Test when dec_input is just under negative boundary (-90.1) - should trigger warning."""
    client = app.test_client()
    response = client.get("/?locus_dec=-90.1")
    # Warning message SHOULD appear for value just under boundary
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_range_out_of_bounds(app):
    """Test when dec_input range has one value out of bounds - should trigger warning."""
    client = app.test_client()
    response = client.get("/?locus_dec=-95 -45")
    # Warning message SHOULD appear when range is out of bounds
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_range_with_signs_out_of_bounds(app):
    """Test when dec_input range with signs has one value out of bounds."""
    client = app.test_client()
    response = client.get("/?locus_dec=<-95 >-45")
    # Warning message SHOULD appear when range is out of bounds
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_range_with_greater_than_sign(app):
    """Test when dec_input is provided as a range with greater than signs - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=>-20 >-15") # >/< signs are ignored in the logic, it's treated as a range from -20 to -15, which is valid input
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_range_with_greater_then_less_than(app):
    """Test when dec_input is provided as >-20 -15 - range from -20 to 15"""
    client = app.test_client()
    response = client.get("/?locus_dec=>-20 <-15") # >/< signs are ignored in the logic, it's treated as a range from -20 to -15, which is valid input
    # Filter warning message should NOT appear for valid input format
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_range_with_less_then_greater_than(app):
    """Test when dec_input is provided as <-20 -15 - range from negative infinity to -15."""
    client = app.test_client()
    response = client.get("/?locus_dec=<-20 -15")
    # Filter warning message should NOT appear for valid input format
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_invalid_range_format_reversed(app):
    """Test when dec_input has reversed range signs (e.g., <-20 -15 but logic expects >)."""
    client = app.test_client()
    response = client.get("/?locus_dec=<-15 -20") #this is also -20 to -15, the order of the numbers doesn't matter, the logic sorts them and the signs are ignored, so this is valid input
    # This is actually a valid format for a range, since the numbers get sorted and lesser/greater signs are ignored in the logic. So this should not trigger a warning.
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_non_numeric_string_with_numbers_at_the_end(app):
    """Test when non-numeric string is provided - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=abc123")
    # Warning message SHOULD appear for invalid input (non-numeric)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_non_numeric_string_with_numbers_at_the_beginning(app):
    """Test when non-numeric string is provided - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=123abc")
    # Warning message SHOULD appear for invalid input (non-numeric)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_leading_decimal_point(app):
    """Test when dec_input is provided with multiple decimal places - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-.56789")
    # Filter warning message should NOT appear for valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_zero(app):
    """Test when dec_input is zero - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=0")
    # Filter warning message should NOT appear for valid input (zero)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_float_zero(app):
    """Test when dec_input is float zero - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=0.0")
    # Filter warning message should NOT appear for valid input (float zero)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_whitespace(app):
    """Test when dec_input is provided with surrounding whitespace - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec= -20.02131 ")
    # Filter warning message should NOT appear for valid input (whitespace trimmed)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_special_characters(app):
    """Test when dec_input contains special characters - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20.02131%")
    # Warning message SHOULD appear for invalid input (special character)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_mixed_case(app):
    """Test when dec_input contains mixed case letters - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20.02131ABC")
    # Warning message SHOULD appear for invalid input (letters)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_nan_like_string(app):
    """Test when dec_input is 'nan' - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=NaN") # NaN is not just a string and therefore not a valid number, it should trigger a warning
    # Warning message SHOULD appear for invalid input (NaN-like string)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_comma_separated(app):
    """Test when dec_input uses comma separator - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20,-15")
    # Warning message SHOULD appear for invalid input (comma instead of space) # IDEA: Support comma as a separator in the future, but for now it should trigger a warning
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_maximum_range(app):
    """Test when dec_input is a very large range - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-90 90")
    # Filter warning message should NOT appear for valid input (large but numeric)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_multiple_spaces(app):
    """Test when dec_input has multiple spaces between numbers - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20   -15")
    # Filter warning message should NOT appear for valid input (multiple spaces)
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_unicode(app):
    """Test when dec_input contains unicode characters - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20.02131α")
    # Warning message SHOULD appear for invalid input (unicode character)
    assert b'<div class="alert alert-warning" role="alert">' in response.data

def test_dec_filter_logic_with_very_long_string_with_many_decimal_places(app):
    """Test when dec_input is a very long string - filter should be applied."""
    client = app.test_client()
    response = client.get("/?locus_dec=-12.345678901234567890") #only the first 5 decimal places are considered, the rest is ignored, still valid input
    assert b'<div class="alert alert-warning" role="alert">' not in response.data

def test_dec_filter_logic_with_trailing_comma(app):
    """Test when dec_input has trailing comma - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-20.02131,")
    # Warning message SHOULD appear for invalid input (trailing comma)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_with_leading_comma(app):
    """Test when dec_input has leading comma - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=,-20.02131")
    # Warning message SHOULD appear for invalid input (leading comma)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_starting_with_very_long_range(app):
    """Test when dec_input is a very long range - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-12345678901234567890 -1")
    # Warning message SHOULD appear for invalid input (range with different lengths)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data

def test_dec_filter_logic_ending_with_very_long_range(app):
    """Test when dec_input is a very long range - filter_warning_message should be added."""
    client = app.test_client()
    response = client.get("/?locus_dec=-1 12345678901234567890")
    # Warning message SHOULD appear for invalid input (range with different lengths)
    assert b'<div class="alert alert-warning" role="alert">' in response.data
    assert b'Dec filter cannot be applied' in response.data