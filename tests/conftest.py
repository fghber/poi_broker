import os
import pytest

"""
Set isolated temp SQLite paths for both alerts and users DBs via env vars
Set testing mode and disable CSRF for test client runs
Create/drop DB tables around each test app lifecycle
"""

@pytest.fixture()
def app(tmp_path, monkeypatch):
    alerts_db = tmp_path / "alerts_test.db"
    users_db = tmp_path / "users_test.db"

    monkeypatch.setenv("SECRET_KEY", "test-secret-key")
    monkeypatch.setenv("FLASK_TESTING", "1")
    monkeypatch.setenv("FLASK_DEBUG", "0")
    monkeypatch.setenv("ALERTS_DB_PATH", str(alerts_db))
    monkeypatch.setenv("USERS_DB_PATH", str(users_db))

    from poi_broker import create_app, db

    app = create_app()
    app.config.update(
        TESTING=True,
        WTF_CSRF_ENABLED=False,
    )

    with app.app_context():
        db.create_all()

    yield app

    with app.app_context():
        db.session.remove()
        db.drop_all()


@pytest.fixture()
def client(app):
    return app.test_client()


@pytest.fixture()
def auth_client(app):
    """Test client pre-logged-in as a verified user."""
    from werkzeug.security import generate_password_hash
    from poi_broker import db
    from poi_broker.models import User

    TEST_EMAIL = "smoketest@example.com"
    TEST_PASSWORD = "SmokeTest123!"

    with app.app_context():
        user = User(
            email=TEST_EMAIL,
            password=generate_password_hash(TEST_PASSWORD),
            name="Smoke Test User",
            email_verified=True,
        )
        db.session.add(user)
        db.session.commit()

    client = app.test_client()
    client.post(
        "/login",
        data={"email": TEST_EMAIL, "password": TEST_PASSWORD},
        follow_redirects=False,
    )
    return client


@pytest.fixture()
def user_factory(app):
    from werkzeug.security import generate_password_hash
    from poi_broker import db
    from poi_broker.models import User

    def factory(
        email: str = "user@example.com",
        password: str = "Password123!",
        name: str = "Test User",
        verified: bool = True,
    ):
        with app.app_context():
            user = User(
                email=email,
                password=generate_password_hash(password),
                name=name,
                email_verified=verified,
            )
            db.session.add(user)
            db.session.commit()
            return user

    return factory


@pytest.fixture()
def secure_app(tmp_path, monkeypatch):
    """App fixture with CSRF and auth rate limiting enabled for security tests."""
    alerts_db = tmp_path / "alerts_secure_test.db"
    users_db = tmp_path / "users_secure_test.db"

    monkeypatch.setenv("SECRET_KEY", "secure-test-secret-key")
    monkeypatch.setenv("FLASK_TESTING", "1")
    monkeypatch.setenv("FLASK_DEBUG", "0")
    monkeypatch.setenv("ALERTS_DB_PATH", str(alerts_db))
    monkeypatch.setenv("USERS_DB_PATH", str(users_db))
    monkeypatch.setenv("RATELIMIT_ENABLED", "1")
    monkeypatch.setenv("RATELIMIT_STORAGE_URI", "memory://")
    monkeypatch.setenv("AUTH_RATE_LIMIT_LOGIN", "3 per minute")
    monkeypatch.setenv("AUTH_RATE_LIMIT_SIGNUP", "2 per minute")
    monkeypatch.setenv("AUTH_RATE_LIMIT_FORGOT_PASSWORD", "2 per minute")
    monkeypatch.setenv("READ_RATE_LIMIT_LAX", "30 per minute")
    monkeypatch.setenv("READ_RATE_LIMIT_MEDIUM", "15 per minute")

    from poi_broker import create_app, db

    app = create_app()
    app.config.update(
        TESTING=True,
        WTF_CSRF_ENABLED=True,
    )

    with app.app_context():
        db.create_all()

    yield app

    with app.app_context():
        db.session.remove()
        db.drop_all()


@pytest.fixture()
def secure_client(secure_app):
    return secure_app.test_client()
