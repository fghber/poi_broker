from datetime import datetime, timedelta, timezone
from werkzeug.security import generate_password_hash

import poi_broker.auth as auth_module
from poi_broker.models import User


def test_normalize_email_validates_and_normalizes_address():
    normalized = auth_module.normalize_email('Test@Example.COM', check_deliverability=False)

    assert normalized == 'Test@example.com'


def test_normalize_email_rejects_invalid_address():
    assert auth_module.normalize_email('not-an-email', check_deliverability=False) is None


def test_signup_post_rejects_invalid_email(client, monkeypatch):
    monkeypatch.setattr(auth_module, 'normalize_email', lambda email, check_deliverability=True: None)

    response = client.post(
        '/signup',
        data={'email': 'invalid', 'name': 'Test User', 'password': 'Password123!'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/signup' in response.location


def test_signup_post_rejects_short_password(client, monkeypatch):
    monkeypatch.setattr(auth_module, 'normalize_email', lambda email, check_deliverability=True: 'user@example.com')

    response = client.post(
        '/signup',
        data={'email': 'user@example.com', 'name': 'Test User', 'password': 'short'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/signup' in response.location


def test_signup_post_rejects_duplicate_email(client, app, monkeypatch, user_factory):
    user_factory(email='duplicate@example.com')
    monkeypatch.setattr(auth_module, 'normalize_email', lambda email, check_deliverability=True: 'duplicate@example.com')

    response = client.post(
        '/signup',
        data={'email': 'duplicate@example.com', 'name': 'Duplicate User', 'password': 'Password123!'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/signup' in response.location


def test_signup_post_creates_user_and_sends_verification_email(client, app, monkeypatch):
    monkeypatch.setattr(auth_module, 'normalize_email', lambda email, check_deliverability=True: email.lower())
    email_calls = []

    def fake_send_email(*args, **kwargs):
        email_calls.append((args, kwargs))
        return True

    monkeypatch.setattr(auth_module, 'send_email', fake_send_email)

    response = client.post(
        '/signup',
        data={'email': 'SignupUser@example.com', 'name': 'Signup User', 'password': 'Password123!'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/login' in response.location
    assert email_calls, 'Expected send_email() to be called'

    with app.app_context():
        user = User.query.filter_by(email='signupuser@example.com').first()
        assert user is not None
        assert user.email_verified is False
        assert user.email_verification_token is not None


def test_login_post_invalid_credentials_redirects_to_forgot_password(client):
    response = client.post(
        '/login',
        data={'email': 'missing@example.com', 'password': 'wrong'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert 'forgot_password=True' in response.location


def test_login_post_unverified_user_redirects_to_login(client, user_factory):
    user_factory(email='pending@example.com', verified=False)

    response = client.post(
        '/login',
        data={'email': 'pending@example.com', 'password': 'Password123!'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/login' in response.location
    assert 'forgot_password' not in response.location


def test_login_post_success_redirects_to_profile(client, user_factory):
    user_factory(email='active@example.com', password='Password123!', verified=True)

    response = client.post(
        '/login',
        data={'email': 'active@example.com', 'password': 'Password123!'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/profile' in response.location


def test_verify_email_invalid_token_redirects_to_signup(client):
    response = client.get('/verify-email/invalid-token', follow_redirects=False)

    assert response.status_code == 302
    assert '/signup' in response.location


def test_verify_email_marks_user_verified(client, app):
    with app.app_context():
        user = User(
            email='verify@example.com',
            password=generate_password_hash('Password123!'),
            name='Verify User',
            email_verified=False,
            email_verification_token='verify-token',
        )
        from poi_broker import db
        db.session.add(user)
        db.session.commit()
        token = user.email_verification_token

    response = client.get(f'/verify-email/{token}', follow_redirects=False)

    assert response.status_code == 302
    assert '/login' in response.location

    with app.app_context():
        user = User.query.filter_by(email='verify@example.com').first()
        assert user.email_verified is True
        assert user.email_verification_token is None


def test_forgot_password_post_generates_reset_token_and_emails_user(client, app, monkeypatch, user_factory):
    user_factory(email='reset@example.com', verified=True)
    monkeypatch.setattr(auth_module, 'send_email', lambda *args, **kwargs: True)

    response = client.post(
        '/forgot-password',
        data={'email': 'reset@example.com'},
        follow_redirects=False,
    )

    assert response.status_code == 302
    assert '/login' in response.location

    with app.app_context():
        user = User.query.filter_by(email='reset@example.com').first()
        assert user.reset_token is not None
        assert user.reset_token_expires is not None
        assert user.reset_token_expires > int(datetime.now(timezone.utc).timestamp())


def test_reset_password_invalid_token_redirects_to_login(client):
    response = client.get('/reset-password/invalid-token', follow_redirects=False)

    assert response.status_code == 302
    assert '/login' in response.location


def test_is_reset_token_expired_and_format_expire_time():
    assert auth_module._is_reset_token_expired(None)
    assert auth_module._is_reset_token_expired('bad')

    future = datetime.now(timezone.utc) + timedelta(hours=1)
    assert auth_module._is_reset_token_expired(future) is False
    assert auth_module._format_expire_time(future).endswith('+00:00')
