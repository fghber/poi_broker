"""Tests for main-table filter bookmark API."""

import json

from werkzeug.security import generate_password_hash


def test_filter_bookmarks_unauthenticated(client):
    r = client.get("/api/filter-bookmarks")
    assert r.status_code == 401
    assert r.is_json
    assert r.get_json().get("error") == "authentication required"

    r = client.post("/api/filter-bookmarks", json={"name": "x", "params": {}})
    assert r.status_code == 401

    r = client.delete("/api/filter-bookmarks/1")
    assert r.status_code == 401


def test_filter_bookmarks_crud(auth_client):
    r = auth_client.post(
        "/api/filter-bookmarks",
        json={
            "name": "My locus",
            "params": {"locus_id": "ANT123", "sort__locus_ra": "asc"},
        },
    )
    assert r.status_code == 201, r.get_data(as_text=True)
    body = r.get_json()
    assert "id" in body
    assert body["name"] == "My locus"
    assert body["params"]["locus_id"] == "ANT123"
    assert "locus_id=ANT123" in body["path"]
    assert "sort__locus_ra=asc" in body["path"]
    bookmark_id = body["id"]

    r = auth_client.get("/api/filter-bookmarks")
    assert r.status_code == 200
    items = r.get_json()["filterBookmarks"]
    assert isinstance(items, list)
    assert len(items) >= 1
    match = next((x for x in items if x["id"] == bookmark_id), None)
    assert match is not None
    assert match["name"] == "My locus"

    r = auth_client.delete(f"/api/filter-bookmarks/{bookmark_id}")
    assert r.status_code == 200
    assert r.is_json
    assert r.get_json().get("status") == "ok"

    r = auth_client.get("/api/filter-bookmarks")
    assert r.status_code == 200
    ids = [x["id"] for x in r.get_json()["filterBookmarks"]]
    assert bookmark_id not in ids


def test_filter_bookmarks_validation(auth_client):
    r = auth_client.post("/api/filter-bookmarks", json={"name": "", "params": {}})
    assert r.status_code == 400
    assert "error" in r.get_json()

    r = auth_client.post("/api/filter-bookmarks", json={"params": {"locus_id": "x"}})
    assert r.status_code == 400

    r = auth_client.post(
        "/api/filter-bookmarks",
        json={"name": "bad", "params": {"not_a_real_param": "1"}},
    )
    assert r.status_code == 400
    assert "unknown" in r.get_json().get("error", "").lower()


def test_filter_bookmarks_omitted_params_empty(auth_client):
    r = auth_client.post("/api/filter-bookmarks", json={"name": "Empty filters"})
    assert r.status_code == 201
    assert r.get_json()["path"] == "/"


def test_filter_bookmarks_cross_user_delete(app, auth_client):
    from poi_broker import db
    from poi_broker.models import FilterBookmark, User

    with app.app_context():
        other = User(
            email="other_bm@example.com",
            password=generate_password_hash("OtherUser123!"),
            name="Other",
            email_verified=True,
        )
        db.session.add(other)
        db.session.commit()

        row = FilterBookmark(
            user_id=other.id,
            name="Theirs",
            query_json=json.dumps({"locus_id": "LOCUS_X"}),
        )
        db.session.add(row)
        db.session.commit()
        foreign_id = row.id

    r = auth_client.delete(f"/api/filter-bookmarks/{foreign_id}")
    assert r.status_code == 404

    with app.app_context():
        still = db.session.get(FilterBookmark, foreign_id)
        assert still is not None
