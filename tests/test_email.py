"""
Tests for the email service module.

Tests the send_email function and _format_expire_time helper with mocked SMTP.
"""

import os
import smtplib
from datetime import datetime, timezone, timedelta
from unittest.mock import MagicMock, patch, ANY

import pytest

import poi_broker.services.email_service as email_module


class TestFormatExpireTime:
    """Tests for the _format_expire_time helper function."""

    def test_none_returns_none(self):
        assert email_module._format_expire_time(None) is None

    def test_datetime_returns_iso_with_utc(self):
        dt = datetime(2026, 7, 25, 12, 34, 56, tzinfo=timezone.utc)
        result = email_module._format_expire_time(dt)
        assert result == "2026-07-25T12:34:56+00:00"

    def test_naive_datetime_assumed_utc(self):
        dt = datetime(2026, 7, 25, 12, 34, 56)  # no tzinfo
        result = email_module._format_expire_time(dt)
        assert result == "2026-07-25T12:34:56+00:00"

    def test_aware_datetime_converted_to_utc(self):
        # UTC+2
        dt = datetime(2026, 7, 25, 14, 34, 56, tzinfo=timezone(timedelta(hours=2)))
        result = email_module._format_expire_time(dt)
        assert result == "2026-07-25T12:34:56+00:00"

    def test_epoch_int_returns_iso(self):
        # 2026-07-25T17:54:56 UTC
        epoch = 1785002096
        result = email_module._format_expire_time(epoch)
        assert result == "2026-07-25T17:54:56+00:00"

    def test_epoch_string_returns_iso(self):
        epoch = "1785002096"
        result = email_module._format_expire_time(epoch)
        assert result == "2026-07-25T17:54:56+00:00"

    def test_invalid_string_returns_none(self):
        assert email_module._format_expire_time("not-a-number") is None

    def test_invalid_type_returns_none(self):
        assert email_module._format_expire_time([]) is None


class TestSendEmail:
    """Tests for the send_email function with mocked SMTP."""

    @pytest.fixture(autouse=True)
    def setup_env(self, monkeypatch):
        """Set up required environment variables for each test."""
        monkeypatch.setenv("SMTP_HOST", "smtp.example.com")
        monkeypatch.setenv("SMTP_PORT", "465")
        monkeypatch.setenv("SMTP_USER", "test@example.com")
        monkeypatch.setenv("SMTP_APP_PASSWORD", "testpass")
        monkeypatch.setenv("SMTP_FROM", "sender@example.com")
        monkeypatch.setenv("LOCAL_HOST", "localhost")

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_success(self, mock_smtp_ssl):
        """Test successful email sending."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        result = email_module.send_email(
            message="Test message",
            to_email="recipient@example.com",
            subject="Test Subject",
        )

        assert result is True
        mock_smtp_ssl.assert_called_once_with(
            "smtp.example.com", 465, local_hostname="localhost", context=ANY
        )
        mock_server.login.assert_called_once_with("test@example.com", "testpass")
        mock_server.send_message.assert_called_once()

        # Verify message headers
        sent_msg = mock_server.send_message.call_args[0][0]
        assert sent_msg["To"] == "recipient@example.com"
        assert sent_msg["From"] == "sender@example.com"
        assert sent_msg["Subject"] == "Test Subject"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_failure_returns_false(self, mock_smtp_ssl):
        """Test that SMTP failure returns False and doesn't raise."""
        mock_smtp_ssl.side_effect = smtplib.SMTPException("Connection refused")

        result = email_module.send_email(
            message="Test message",
            to_email="recipient@example.com",
        )

        assert result is False

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_requires_to_email(self, mock_smtp_ssl):
        """Test that None to_email raises ValueError."""
        with pytest.raises(ValueError, match="Recipient email address required"):
            email_module.send_email("Test message", None)

        mock_smtp_ssl.assert_not_called()

    def test_send_email_no_sender_raises_runtime_error(self, monkeypatch):
        """Test that missing sender config raises RuntimeError."""
        monkeypatch.delenv("SMTP_FROM", raising=False)
        monkeypatch.delenv("SMTP_USER", raising=False)

        with pytest.raises(RuntimeError, match="No sender address configured"):
            email_module.send_email("Test message", "recipient@example.com")

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_auto_html_includes_expire_section(self, mock_smtp_ssl):
        """Test auto-generated HTML includes Expires section when expire_time given."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        expire_dt = datetime(2026, 7, 25, 12, 34, 56, tzinfo=timezone.utc)
        email_module.send_email(
            message="Reset your password",
            to_email="recipient@example.com",
            expire_time=expire_dt,
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        html_part = sent_msg.get_body(preferencelist=("html",))
        html_content = html_part.get_content()
        assert "Expires:" in html_content
        assert "2026-07-25T12:34:56+00:00 UTC" in html_content

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_explicit_html_used_verbatim(self, mock_smtp_ssl):
        """Test that explicit html_text is used without auto expire section."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        custom_html = "<html><body><p>Custom HTML</p></body></html>"
        expire_dt = datetime(2026, 7, 25, 12, 34, 56, tzinfo=timezone.utc)
        email_module.send_email(
            message="Reset your password",
            to_email="recipient@example.com",
            html_text=custom_html,
            expire_time=expire_dt,
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        html_part = sent_msg.get_body(preferencelist=("html",))
        html_content = html_part.get_content()
        # EmailMessage adds a trailing newline to content
        assert html_content.rstrip("\n") == custom_html
        assert "Expires:" not in html_content

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_subject_defaults_to_env(self, mock_smtp_ssl):
        """Test subject falls back to SMTP_SUBJECT env var."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        assert sent_msg["Subject"] == "Notification from POI Broker"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_subject_uses_smtp_subject_env(self, mock_smtp_ssl, monkeypatch):
        """Test subject uses SMTP_SUBJECT env var when set."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server
        monkeypatch.setenv("SMTP_SUBJECT", "Custom Subject")

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        assert sent_msg["Subject"] == "Custom Subject"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_from_email_override(self, mock_smtp_ssl):
        """Test from_email parameter overrides env vars."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
            from_email="override@example.com",
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        assert sent_msg["From"] == "override@example.com"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_plain_text_part_set(self, mock_smtp_ssl):
        """Test that plain text part is set correctly."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Plain text body",
            to_email="recipient@example.com",
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        plain_part = sent_msg.get_body(preferencelist=("plain",))
        # EmailMessage adds a trailing newline
        assert plain_part.get_content().rstrip("\n") == "Plain text body"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_html_alternative_added(self, mock_smtp_ssl):
        """Test that HTML alternative is added with correct subtype."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Plain text",
            to_email="recipient@example.com",
        )

        sent_msg = mock_server.send_message.call_args[0][0]
        # Check that there are two parts (plain + html)
        parts = list(sent_msg.iter_parts())
        assert len(parts) == 2
        # First part should be plain text
        assert parts[0].get_content_type() == "text/plain"
        # Second part should be HTML
        assert parts[1].get_content_type() == "text/html"

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_login_called_when_credentials_present(self, mock_smtp_ssl):
        """Test SMTP login is called when both user and password are set."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        mock_server.login.assert_called_once_with("test@example.com", "testpass")

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_no_login_when_credentials_missing(self, mock_smtp_ssl, monkeypatch):
        """Test SMTP login is skipped when credentials are missing."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server
        monkeypatch.delenv("SMTP_USER", raising=False)
        monkeypatch.delenv("SMTP_APP_PASSWORD", raising=False)

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        mock_server.login.assert_not_called()

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_uses_smtp_host_port_from_env(self, mock_smtp_ssl, monkeypatch):
        """Test SMTP host and port are read from environment."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server
        monkeypatch.setenv("SMTP_HOST", "custom.smtp.com")
        monkeypatch.setenv("SMTP_PORT", "587")

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        mock_smtp_ssl.assert_called_once_with(
            "custom.smtp.com", 587, local_hostname="localhost", context=ANY
        )

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_uses_local_host_from_env(self, mock_smtp_ssl, monkeypatch):
        """Test LOCAL_HOST env var is used for local_hostname."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server
        monkeypatch.setenv("LOCAL_HOST", "myhost.example.com")

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        mock_smtp_ssl.assert_called_once_with(
            "smtp.example.com", 465, local_hostname="myhost.example.com", context=ANY
        )

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_ssl_context_created(self, mock_smtp_ssl):
        """Test that SSL context is created with default settings."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        email_module.send_email(
            message="Test",
            to_email="recipient@example.com",
        )

        # Verify SSL context was passed (4th positional arg or keyword arg)
        call_args = mock_smtp_ssl.call_args
        assert call_args is not None
        # context can be 4th positional or keyword argument
        if len(call_args[0]) > 3:
            assert call_args[0][3] is not None  # ssl context as positional
        else:
            assert call_args[1].get("context") is not None  # ssl context as keyword

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_logs_on_success(self, mock_smtp_ssl, caplog):
        """Test that success is logged."""
        mock_server = MagicMock()
        mock_smtp_ssl.return_value.__enter__.return_value = mock_server

        with caplog.at_level("INFO"):
            email_module.send_email(
                message="Test",
                to_email="recipient@example.com",
                subject="Test Subject",
            )

        assert "Sent email to recipient@example.com (subject=Test Subject)" in caplog.text

    @patch("poi_broker.services.email_service.smtplib.SMTP_SSL")
    def test_send_email_logs_on_failure(self, mock_smtp_ssl, caplog):
        """Test that failure is logged with exception."""
        mock_smtp_ssl.side_effect = smtplib.SMTPException("Connection refused")

        with caplog.at_level("ERROR"):
            email_module.send_email(
                message="Test",
                to_email="recipient@example.com",
            )

        assert "Failed to send email to recipient@example.com" in caplog.text
        assert "Connection refused" in caplog.text