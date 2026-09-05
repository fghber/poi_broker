"""
Email service module for sending multipart emails via SMTP.

Provides a reusable send_email function that can be used across the application
and by external tools (e.g., watchlist_digest.py).
"""

import os
import smtplib
import ssl
import logging
from email.message import EmailMessage
from datetime import datetime, timezone
from typing import Optional


logger = logging.getLogger(__name__)


def _format_expire_time(expire_time: Optional[object]) -> Optional[str]:
    """
    Format an expiration time for display in email body.
    
    Args:
        expire_time: Can be a datetime object, an epoch timestamp (int/float/str),
                     or None.
    
    Returns:
        ISO format string in UTC (e.g., "2026-07-25T12:34:56+00:00") or None if
        input is None or unparseable.
    """
    if expire_time is None:
        return None
    if isinstance(expire_time, datetime):
        if expire_time.tzinfo is None:
            expire_time = expire_time.replace(tzinfo=timezone.utc)
        else:
            expire_time = expire_time.astimezone(timezone.utc)
        return expire_time.isoformat()
    try:
        return datetime.fromtimestamp(int(expire_time), tz=timezone.utc).isoformat()
    except (TypeError, ValueError, OSError, OverflowError):
        return None


def send_email(
    message: str,
    to_email: str,
    subject: Optional[str] = None,
    html_text: Optional[str] = None,
    from_email: Optional[str] = None,
    expire_time: Optional[object] = None,
) -> bool:
    """
    Send a multipart email (plain text + optional HTML) via SMTP_SSL.
    
    Backwards-compatible: legacy calls use send_email(message, email).
    Prefer setting SMTP env vars: SMTP_HOST, SMTP_PORT, SMTP_USER, SMTP_APP_PASSWORD, SMTP_FROM.
    
    Args:
        message: Plain text message body (required).
        to_email: Recipient email address (required).
        subject: Email subject line. Defaults to SMTP_SUBJECT env or "Notification from POI Broker".
        html_text: Optional HTML version of the message. If not provided, a simple
                   HTML wrapper is generated from the plain text.
        from_email: Optional sender address. Defaults to SMTP_FROM env or SMTP_USER.
        expire_time: Optional datetime/epoch when the link/token expires. If provided
                     and no explicit html_text is given, an "Expires:" line is added
                     to the auto-generated HTML.
    
    Returns:
        True on success, False on failure (logs exception).
    
    Raises:
        ValueError: If to_email is None.
        RuntimeError: If no sender address is configured (SMTP_FROM or SMTP_USER).
    """
    # Backwards compatibility: if called as send_email(message, email)
    if to_email is None:
        raise ValueError("Recipient email address required as second argument.")

    plain_text = str(message)
    if subject is None:
        subject = os.environ.get("SMTP_SUBJECT", "Notification from POI Broker")

    # If no explicit HTML provided, create a simple HTML version
    if html_text is None:
        expire_section = ""
        expire_display = _format_expire_time(expire_time)
        if expire_display:
            expire_section = f"<p><strong>Expires:</strong> {expire_display} UTC</p>"
        
        html_text = (
            "<html><body>"
            f"<h3>{subject}</h3>"
            f"{expire_section}"
            f"<hr><pre style='white-space:pre-wrap'>{plain_text}</pre>"
            "</body></html>"
        )

    SMTP_HOST = os.environ.get("SMTP_HOST", "smtp.gmail.com")
    SMTP_PORT = int(os.environ.get("SMTP_PORT", 465))
    SMTP_USER = os.environ.get("SMTP_USER")  # required for authenticated SMTP
    SMTP_PASS = os.environ.get("SMTP_APP_PASSWORD")
    FROM = from_email or os.environ.get("SMTP_FROM") or SMTP_USER
    LOCAL_HOST = os.environ.get("LOCAL_HOST", "localhost")

    if not FROM:
        raise RuntimeError("No sender address configured (SMTP_FROM or SMTP_USER)")

    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = FROM
    msg["To"] = to_email
    msg.set_content(plain_text)
    msg.add_alternative(html_text, subtype="html")

    try:
        with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT, local_hostname=LOCAL_HOST, context=ssl.create_default_context()) as server:
            if SMTP_USER and SMTP_PASS:
                server.login(SMTP_USER, SMTP_PASS)
            server.send_message(msg)
        logger.info("Sent email to %s (subject=%s)", to_email, subject)
        return True
    except Exception as exc:
        logger.exception("Failed to send email to %s: %s", to_email, exc)
        return False