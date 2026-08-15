from flask import Blueprint, render_template, request, redirect, url_for, flash, current_app, make_response
from flask_login import login_user, logout_user, login_required
from werkzeug.security import generate_password_hash, check_password_hash
from sqlalchemy.exc import IntegrityError
from .models import User
from . import db
import secrets
from datetime import datetime, timedelta, timezone
import os
import logging
from email_validator import validate_email, EmailNotValidError
from . import limiter
from .services.email_service import send_email, _format_expire_time
import hashlib


def hash_token(token: str) -> str:
    """Return a SHA-256 hex digest of a token for safe at-rest storage."""
    return hashlib.sha256(token.encode('utf-8')).hexdigest()


def _utc_now_epoch() -> int:
    return int(datetime.now(timezone.utc).timestamp())

def _is_reset_token_expired(expires_at) -> bool:
    if expires_at is None:
        return True
    if isinstance(expires_at, datetime):
        expires_epoch = int(expires_at.timestamp())
    else:
        try:
            expires_epoch = int(expires_at)
        except (TypeError, ValueError):
            return True
    return expires_epoch < _utc_now_epoch()

logger = logging.getLogger(__name__)
auth_blueprint = Blueprint('auth', __name__)

# Same copy + redirect for new signup, duplicate email, and uniqueness races
# so the response is not an account-existence oracle.
SIGNUP_GENERIC_NOTICE = (
    'If this email is available, a verification message will be sent.'
)

#IDEA: Improve emails (body and subject), add HTML version, perhaps use a proper email template, etc.


def _signup_accepted_response():
    flash(SIGNUP_GENERIC_NOTICE)
    return redirect(url_for('auth.login'))

@auth_blueprint.route('/login')
def login():
    show_reset = request.args.get('forgot_password', default=False, type=bool)
    return render_template('login.html', forgot_password=show_reset)

@auth_blueprint.route('/login', methods=['POST'])
@limiter.limit(lambda: current_app.config.get('AUTH_RATE_LIMIT_LOGIN', '10 per minute'))
def login_post():
    email_input = request.form.get('email')
    password = request.form.get('password')
    remember = True if request.form.get('remember') else False

    if email_input is None or password is None:
        flash('Please provide an email and a password and try again.')
        return redirect(url_for('auth.login')) # reload the page

    # Normalize email for lookup (no deliverability check on login)
    email = normalize_email(email_input, check_deliverability=False)
    if not email:
        flash('Please provide a valid email address and try again.')
        return redirect(url_for('auth.login'))

    user = User.query.filter_by(email=email).first()

    # check if user actually exists & provided the right password (compared to hashed password in database)
    if not user or not check_password_hash(user.password, password):
        flash('Please check your login details and try again.')
        return redirect(url_for('auth.login', forgot_password=True)) # if user doesn't exist or password is wrong, reload the page
    
    # Check if email is verified
    if not user.email_verified:
        flash('Please verify your email before logging in. Check your inbox for the verification link.')
        return redirect(url_for('auth.login'))
    
    # otherwise, we know the user has the right credentials
    login_user(user, remember=remember)
    return redirect(url_for('main.profile'))

@auth_blueprint.route('/signup')
def signup():
    return render_template('signup.html')

def normalize_email(email: str, check_deliverability: bool = False) -> str | None:
    """
    Validates and normalizes an email address for correct syntax.
    Returns the normalized email if valid, None otherwise.
    Uses email_validator for RFC 5321 compliance and domain normalization.
    
    Args:
        email: Email address to validate and normalize.
        check_deliverability: If True, verify MX records (use for account creation as per email_validator docs).
                             If False, skip MX checks (use for login/lookups).
    """
    try:
        valid = validate_email(email, check_deliverability=check_deliverability)
        return valid.normalized
    except EmailNotValidError:
        return None

@auth_blueprint.route('/signup', methods=['POST'])
@limiter.limit(lambda: current_app.config.get('AUTH_RATE_LIMIT_SIGNUP', '5 per hour'))
def signup_post():
    email_input = request.form.get('email', '').strip()
    name = request.form.get('name', '').strip()
    password = request.form.get('password', '').strip()

    # Validate and normalize email (check deliverability for new accounts)
    email = normalize_email(email_input, check_deliverability=True)
    if not email:
        flash('Invalid email address.')
        return redirect(url_for('auth.signup'))
    
    # Validate password length
    if not password or len(password) < 8:
        flash('Password must be at least 8 characters.')
        return redirect(url_for('auth.signup'))
    
    if not name:
        flash('Name is required.')
        return redirect(url_for('auth.signup'))

    user = User.query.filter_by(email=email).first() # if this returns a user, then the email already exists in database

    if user:
        return _signup_accepted_response()

    # Create user but mark as unverified
    raw_verification_token = secrets.token_urlsafe(32)
    new_user = User(
        email=email,
        name=name,
        password=generate_password_hash(password),
        email_verified=False,
        email_verification_token=hash_token(raw_verification_token)
    )

    # Send verification email
    verification_link = url_for('auth.verify_email', token=raw_verification_token, _external=True)
    result = send_email(
        message=f'Please verify your email by clicking here: {verification_link}',
        to_email=email,
        subject='Verify your POI Broker account',
        html_text=(
            '<html><body>'
            f'<h3>Welcome to POI Broker!</h3>'
            f'<p>Please verify your email by clicking the link below:</p>'
            f'<p><a href="{verification_link}">Verify Email</a></p>'
            f'<p>Or copy and paste this link: {verification_link}</p>'
            '</body></html>'
        )
    )

    if not result:
        flash('Failed to send verification email. Please try signing up again later.')
        return redirect(url_for('auth.signup'))
    
    #only add the user to the database if the email was sent successfully to avoid creating unverified accounts with invalid emails
    db.session.add(new_user)
    try:
        db.session.commit()
    except IntegrityError:
        db.session.rollback()  # uniqueness race; same response as a duplicate signup
        return _signup_accepted_response()
    except Exception as e:
        db.session.rollback()
        logger.error(f'Database error during commit: {str(e)}', exc_info=True)
        flash('Failed to create user. Please try signing up again later.') 
        return redirect(url_for('auth.signup'))

    return _signup_accepted_response()

@auth_blueprint.route('/verify-email/<token>')
def verify_email(token):
    """Verify user email via token link."""
    user = User.query.filter_by(email_verification_token=hash_token(token)).first()
    
    if not user:
        flash('Invalid or expired verification link.')
        return redirect(url_for('auth.signup'))
    
    if user.email_verified:
        flash('Email is already verified. You can now log in.')
        return redirect(url_for('auth.login'))
    
    # Mark email as verified and clear the token
    user.email_verified = True
    user.email_verification_token = None
    try:
        db.session.commit()
    except Exception as e:
        db.session.rollback()
        logger.error(f'Database error during commit: {str(e)}', exc_info=True)
        flash('Failed to verify email. Please try again later.')
        return redirect(url_for('auth.signup'))
    
    flash('Email verified successfully! You can now log in.')
    return redirect(url_for('auth.login'))


@auth_blueprint.route('/logout', methods=['POST'])
@login_required
def logout():
    logout_user()
    resp = make_response(redirect(url_for('main.start')))
    # Explicitly clear session cookie with secure flags
    resp.set_cookie(
        current_app.config.get('SESSION_COOKIE_NAME', 'session'),
        '',
        expires=0,
        httponly=True,
        secure=current_app.config.get('SESSION_COOKIE_SECURE', True),
        samesite=current_app.config.get('SESSION_COOKIE_SAMESITE', 'Lax')
    )
    return resp

@auth_blueprint.route('/forgot-password')
def forgot_password():
    return render_template('forgot_password.html')

@auth_blueprint.route('/forgot-password', methods=['POST'])
@limiter.limit(lambda: current_app.config.get('AUTH_RATE_LIMIT_FORGOT_PASSWORD', '5 per hour'))
def forgot_password_post():
    email_input = request.form.get('email')
    
    # Normalize email for lookup (no deliverability check on password reset)
    email = normalize_email(email_input, check_deliverability=False) if email_input else None
    
    user = User.query.filter_by(email=email).first() if email else None
    
    if user:
        # Generate reset token
        raw_reset_token = secrets.token_urlsafe(32)
        user.reset_token = hash_token(raw_reset_token)
        user.reset_token_expires = _utc_now_epoch() + int(timedelta(hours=1).total_seconds())
        try:
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            logger.error(f'Database error during commit: {str(e)}', exc_info=True)
            flash('Failed to generate password reset link. Please try again later.')
            return redirect(url_for('auth.forgot_password'))
        
        # Send email with reset link and expiration time
        result = send_email(
            f"Reset your password using the following link: {url_for('auth.reset_password', token=raw_reset_token, _external=True)}", 
            email,
            expire_time=user.reset_token_expires
        )
        if not result:
            flash('Failed to send password reset email. Please try again later.')
        else:
            flash('Password reset link sent to your email')
    
    return redirect(url_for('auth.login'))

@auth_blueprint.route('/reset-password/<token>')
def reset_password(token):
    user = User.query.filter_by(reset_token=hash_token(token)).first()
    
    if not user or _is_reset_token_expired(user.reset_token_expires):
        flash('Invalid or expired reset token')
        return redirect(url_for('auth.login'))
    
    return render_template('reset_password.html', token=token)

@auth_blueprint.route('/reset-password/<token>', methods=['POST'])
@limiter.limit(lambda: current_app.config.get('AUTH_RATE_LIMIT_RESET_PASSWORD', '10 per hour'))
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

    user = User.query.filter_by(reset_token=hash_token(token)).first()
    if not user or _is_reset_token_expired(user.reset_token_expires):
        flash('Invalid or expired reset token')
        return redirect(url_for('auth.login'))

    # Hash and store new password, clear the token fields
    user.password = generate_password_hash(password)
    user.reset_token = None
    user.reset_token_expires = None
    try:
        db.session.commit()
    except Exception as e:
        db.session.rollback()
        logger.error(f'Database error during commit: {str(e)}', exc_info=True)
        flash('Failed to update password. Please try again later.')
        return redirect(url_for('auth.reset_password', token=token))

    flash('Password updated. Please log in.')
    return redirect(url_for('auth.login'))


@auth_blueprint.route('/security')
@login_required
def security():
    return render_template('security.html')

@auth_blueprint.route('/change-password', methods=['POST'])
@login_required
@limiter.limit(lambda: current_app.config.get('AUTH_RATE_LIMIT_CHANGE_PASSWORD', '10 per hour'))
def change_password():
    """
    Allow authenticated users to change their password.
    Requires current password verification for security.
    """
    from flask_login import current_user
    
    current_password = request.form.get('current_password')
    new_password = request.form.get('new_password')
    new_password_confirm = request.form.get('new_password_confirm')
    
    # Validate that all fields are provided
    if not current_password or not new_password or not new_password_confirm:
        flash('Please provide current password and new password.')
        return redirect(url_for('auth.security'))
    
    # Verify current password
    if not check_password_hash(current_user.password, current_password):
        flash('Current password is incorrect.')
        return redirect(url_for('auth.security'))
    
    # Validate new password matches confirmation
    if new_password != new_password_confirm:
        flash('New passwords do not match.')
        return redirect(url_for('auth.security'))
    
    # Validate password length
    if len(new_password) < 8:
        flash('New password must be at least 8 characters.')
        return redirect(url_for('auth.security'))
    
    # Check that new password is different from current
    if check_password_hash(current_user.password, new_password):
        flash('New password must be different from current password.')
        return redirect(url_for('auth.security'))
    
    # Update password
    current_user.password = generate_password_hash(new_password)
    try:
        db.session.commit()
    except Exception as e:
        db.session.rollback()
        logger.error(f'Database error during commit: {str(e)}', exc_info=True)
        flash('Failed to change password. Please try again later.')
        return redirect(url_for('auth.security'))

    flash('Password changed successfully.')
    return redirect(url_for('auth.security'))