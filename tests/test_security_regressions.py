import re

from werkzeug.security import generate_password_hash


def _extract_hidden_csrf_token(html: str) -> str:
    match = re.search(r'name="csrf_token"\s+value="([^"]+)"', html)
    assert match, "Could not find hidden csrf_token input"
    return match.group(1)


def _extract_meta_csrf_token(html: str) -> str:
    match = re.search(r'name="csrf-token"\s+content="([^"]+)"', html)
    assert match, "Could not find csrf-token meta tag"
    return match.group(1)


def _create_verified_user(db, User, email: str, password: str, name: str):
    user = User(
        email=email,
        password=generate_password_hash(password),
        name=name,
        email_verified=True,
    )
    db.session.add(user)
    db.session.commit()
    return user.id


def _force_login(client, user_id: int) -> None:
    with client.session_transaction() as session:
        session["_user_id"] = str(user_id)
        session["_fresh"] = True


def test_api_csrf_missing_token_rejected(secure_app, secure_client):
    from poi_broker import db
    from poi_broker.models import User

    with secure_app.app_context():
        user_id = _create_verified_user(db, User, "csrf1@example.com", "Password123!", "CSRF User")

    _force_login(secure_client, user_id)

    response = secure_client.post("/api/favorite", json={"locusId": "locus-csrf-1", "fav": True, "groupId": None})
    assert response.status_code == 400
    assert response.is_json
    payload = response.get_json()
    assert payload.get("error") == "csrf validation failed"


def test_api_csrf_header_allows_authenticated_post(secure_app, secure_client):
    from poi_broker import db
    from poi_broker.models import User

    with secure_app.app_context():
        user_id = _create_verified_user(db, User, "csrf2@example.com", "Password123!", "CSRF User 2")

    _force_login(secure_client, user_id)

    html = secure_client.get("/").get_data(as_text=True)
    csrf_token = _extract_meta_csrf_token(html)

    response = secure_client.post(
        "/api/favorite",
        json={"locusId": "locus-csrf-2", "fav": True, "groupId": None},
        headers={"X-CSRFToken": csrf_token},
    )
    assert response.status_code == 200
    assert response.is_json
    assert response.get_json().get("status") == "ok"


def test_cross_user_cannot_delete_other_watchlist(secure_app, secure_client):
    from poi_broker import db
    from poi_broker.models import User, Watchlist

    with secure_app.app_context():
        owner_id = _create_verified_user(db, User, "owner@example.com", "Password123!", "Owner")
        intruder_id = _create_verified_user(db, User, "intruder@example.com", "Password123!", "Intruder")

        wl = Watchlist(
            user_id=owner_id,
            name="Owner Watchlist",
            rules_json='{"condition":"AND","rules":[{"field":"featuretable.alert_id","operator":"is_not_null"}]}',
            sql_where="featuretable.alert_id IS NOT NULL",
            created_at=1711929600,
        )
        db.session.add(wl)
        db.session.commit()
        watchlist_id = wl.id

    _force_login(secure_client, intruder_id)
    html = secure_client.get("/").get_data(as_text=True)
    csrf_token = _extract_meta_csrf_token(html)

    response = secure_client.delete(
        f"/api/watchlist/{watchlist_id}",
        headers={"X-CSRFToken": csrf_token},
    )
    assert response.status_code == 404


    with secure_app.app_context():
        still_exists = Watchlist.query.filter_by(id=watchlist_id).first()
        assert still_exists is not None


def test_login_post_rate_limited_after_retries(secure_client):
    login_page = secure_client.get("/login")
    csrf_token = _extract_hidden_csrf_token(login_page.get_data(as_text=True))

    statuses = []
    for _ in range(4):
        r = secure_client.post(
            "/login",
            data={
                "email": "nobody@example.com",
                "password": "wrong-password",
                "csrf_token": csrf_token,
            },
            follow_redirects=False,
        )
        statuses.append(r.status_code)

    assert statuses[:3] == [302, 302, 302]
    assert statuses[3] == 429
