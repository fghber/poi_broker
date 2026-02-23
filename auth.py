from flask import Blueprint, render_template, request, redirect, url_for, flash
from flask_login import login_user, logout_user, login_required
from werkzeug.security import generate_password_hash, check_password_hash
from .models import User
from . import db
import secrets
from datetime import datetime, timedelta, timezone
import os
import smtplib
import socket
import logging
from email.message import EmailMessage

auth_blueprint = Blueprint('auth', __name__)

@auth_blueprint.route('/login')
def login():
    return render_template('login.html')

@auth_blueprint.route('/login', methods=['POST'])
def login_post():
    email = request.form.get('email')
    password = request.form.get('password')
    remember = True if request.form.get('remember') else False

    user = User.query.filter_by(email=email).first()

    # check if user actually exists
    # take the user supplied password, hash it, and compare it to the hashed password in database
    if not user or not check_password_hash(user.password, password):
        flash('Please check your login details and try again.')
        return redirect(url_for('auth.login')) # if user doesn't exist or password is wrong, reload the page
    
    # otherwise, we know the user has the right credentials
    login_user(user, remember=remember)
    return redirect(url_for('main.profile'))

@auth_blueprint.route('/signup')
def signup():
    return render_template('signup.html')

@auth_blueprint.route('/signup', methods=['POST'])
def signup_post():
    email = request.form.get('email')
    name = request.form.get('name')
    password = request.form.get('password')

    user = User.query.filter_by(email=email).first() # if this returns a user, then the email already exists in database

    if user: # if a user is found, we want to redirect back to signup page so user can try again  
        flash('Email address already exists')
        return redirect(url_for('auth.signup'))

    new_user = User(email=email, name=name, password=generate_password_hash(password))
    db.session.add(new_user)
    db.session.commit()

    return redirect(url_for('auth.login'))

@auth_blueprint.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('main.start'))



@auth_blueprint.route('/forgot-password')
def forgot_password():
    return render_template('forgot_password.html')

@auth_blueprint.route('/forgot-password', methods=['POST'])
def forgot_password_post():
    email = request.form.get('email')
    user = User.query.filter_by(email=email).first()
    
    if user:
        # Generate reset token
        user.reset_token = secrets.token_urlsafe(32)
        user.reset_token_expires = datetime.now(timezone.utc) + timedelta(hours=1)
        db.session.commit()
        
        # Send email with reset link (implement email sending)
        send_email(f"Reset your password using the following link: {url_for('auth.reset_password', token=user.reset_token, _external=True)}", email)
        flash('Password reset link sent to your email')
    
    return redirect(url_for('auth.login'))

@auth_blueprint.route('/reset-password/<token>')
def reset_password(token):
    user = User.query.filter_by(reset_token=token).first()
    
    if not user or user.reset_token_expires < datetime.now(timezone.utc):
        flash('Invalid or expired reset token')
        return redirect(url_for('auth.login'))
    
    return render_template('reset_password.html', token=token)

@auth_blueprint.route('/reset-password/<token>', methods=['POST'])
def reset_password_post(token):
    password = request.form.get('password')
    password_confirm = request.form.get('password_confirm')

    if not password or not password_confirm:
        flash('Please provide and confirm your new password')
        return redirect(url_for('auth.reset_password', token=token))

    if password != password_confirm:
        flash('Passwords do not match')
        return redirect(url_for('auth.reset_password', token=token))

    if len(password) < 8:
        flash('Password must be at least 8 characters')
        return redirect(url_for('auth.reset_password', token=token))

    user = User.query.filter_by(reset_token=token).first()
    if not user or not user.reset_token_expires or user.reset_token_expires < datetime.now(timezone.utc):
        flash('Invalid or expired reset token')
        return redirect(url_for('auth.login'))

    # Hash and store new password, clear the token fields
    user.password = generate_password_hash(password)
    user.reset_token = None
    user.reset_token_expires = None
    db.session.commit()

    flash('Password updated. Please log in.')
    return redirect(url_for('auth.login'))

# ...existing code...
import os
import smtplib
import socket
import logging
from email.message import EmailMessage
from datetime import datetime as _dt

logger = logging.getLogger(__name__)

def send_email(message, to_email=None, subject=None, html_text=None, from_email=None):
    """
    Send a multipart email (plain text + optional HTML).
    Backwards-compatible: legacy calls use send_email(message, email).
    Prefer setting SMTP env vars: SMTP_HOST, SMTP_PORT, SMTP_USER, SMTP_PASS, SMTP_FROM.
    """
    # Backwards compatibility: if called as send_email(message, email)
    if to_email is None:
        raise ValueError("Recipient email address required as second argument (to_email).")

    plain_text = str(message)
    if subject is None:
        subject = os.environ.get("SMTP_SUBJECT", "Notification from POI Broker")

    # If no explicit HTML provided, create a simple HTML version
    if html_text is None:
        html_text = (
            "<html><body>"
            f"<h3>{subject}</h3>"
            f"<p><strong>Time:</strong> {datetime.now(timezone.utc).isoformat()} UTC</p>"
            f"<hr><pre style='white-space:pre-wrap'>{plain_text}</pre>"
            "</body></html>"
        )

    SMTP_HOST = os.environ.get("SMTP_HOST", "smtp.gmail.com")
    SMTP_PORT = int(os.environ.get("SMTP_PORT", 465))
    SMTP_USER = os.environ.get("SMTP_USER")  # required for authenticated SMTP
    SMTP_PASS = os.environ.get("SMTP_APPP_PASSWORD")
    FROM = from_email or os.environ.get("SMTP_FROM") or SMTP_USER

    if not FROM:
        raise RuntimeError("No sender address configured (SMTP_FROM or SMTP_USER)")

    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = FROM
    msg["To"] = to_email
    msg.set_content(plain_text)
    msg.add_alternative(html_text, subtype="html")

    try:
        with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT) as server:
            if SMTP_USER and SMTP_PASS:
                server.login(SMTP_USER, SMTP_PASS)
            server.send_message(msg)
        logger.info("Sent email to %s (subject=%s)", to_email, subject)
        return True
    except Exception as exc:
        logger.exception("Failed to send email to %s: %s", to_email, exc)
        return False