from datetime import timezone
#from astropy.time import Time
from unittest.mock import patch

from poi_broker.app import _format_mjd_cached


class TestFormatMJD:
    """Test cases for _format_mjd_cached function."""

    def test_basic_mjd_conversion(self):
        """Test basic MJD to datetime conversion."""
        # MJD 60000.0 should convert to a specific date
        result = _format_mjd_cached(60000.0)
        assert isinstance(result, str)
        assert len(result) > 0
        # Verify format is correct: YYYY-MM-DD HH:MM:SS
        assert '-' in result
        assert ':' in result

    def test_mjd_zero(self):
        """Test with MJD = 0 (reference epoch)."""
        result = _format_mjd_cached(0.0)
        assert isinstance(result, str)
        # MJD 0 corresponds to 1858-11-17 00:00:00 UTC
        assert result == '1858-11-17 00:00:00'

    def test_mjd_integer(self):
        """Test with integer MJD value."""
        result = _format_mjd_cached(60000)
        assert isinstance(result, str)
        # Should match float version
        assert _format_mjd_cached(60000.0) == result

    def test_mjd_with_seconds(self):
        """Test MJD with fractional seconds."""
        result = _format_mjd_cached(60000.123456)
        assert isinstance(result, str)
        assert result == '2023-02-25 02:57:46'

    def test_caching_behavior(self):
        """Test that the function uses the cache and does not re-calculate."""
        # 1. Clear the cache to guarantee a clean state
        _format_mjd_cached.cache_clear()        
        # 2. First call: Should be a cache MISS
        _format_mjd_cached(60000.5)        
        # 3. Second call: Should be a cache HIT
        _format_mjd_cached(60000.5)        
        # 4. Verify results
        info = _format_mjd_cached.cache_info()
        assert info.misses == 1, f"Expected 1 miss, but got {info.misses}"
        assert info.hits == 1, f"Expected 1 hit, but got {info.hits}"

    def test_different_mjd_values(self):
        """Test with multiple different MJD values."""
        mjd_values = [0.0, 50000.0, 60000.0, 70000.0, 80000.0]
        results = [_format_mjd_cached(mjd) for mjd in mjd_values]
        # All should be valid strings with expected format
        for result in results:
            assert isinstance(result, str)
            assert result.count('-') == 2  # Date separator
            assert result.count(':') == 2  # Time separator