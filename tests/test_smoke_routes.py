"""
Smoke test coverage:

Public pages: /, /help, /contact
Favorites API (unauth behavior)
Visual query routes (login-required redirect behavior)
Watchlist routes (login-required redirect behavior)
Lightcurve/features/crossmatches route health and expected statuses
Test run result:

Execute command via workspace Python env: `pytest -q` or `python -m pytest -q`
"""

def test_public_pages_smoke(client):
    for path in ["/", "/help", "/contact"]:
        response = client.get(path)
        assert response.status_code == 200, f"Expected 200 for {path}, got {response.status_code}"


def test_favorites_api_smoke_unauthenticated(client):
    response = client.get("/api/favorite", query_string={"locusId": "locus-1"})
    assert response.status_code == 401
    assert response.is_json
    assert response.get_json()["error"] == "authentication required"

    response = client.get("/api/favorites")
    assert response.status_code == 401
    assert response.is_json
    assert response.get_json()["error"] == "authentication required"

    response = client.get("/api/favorite-groups")
    assert response.status_code == 401


def test_visual_query_routes_require_login(client):
    response = client.get("/visual_query")
    assert response.status_code in (301, 302)

    response = client.post("/api/preview-query", json={"rules": []})
    assert response.status_code == 401
    assert response.is_json

    response = client.post("/api/export-query", json={"rules": []})
    assert response.status_code == 401
    assert response.is_json


def test_watchlist_routes_require_login(client):
    response = client.post("/api/watchlist", json={"name": "x", "rules": {"condition": "AND", "rules": []}})
    assert response.status_code == 401
    assert response.is_json

    response = client.get("/api/watchlist")
    assert response.status_code == 401
    assert response.is_json

    response = client.delete("/api/watchlist/1")
    assert response.status_code == 401
    assert response.is_json


def test_ui_protected_route_sets_flash_message(client):
    response = client.get('/visual_query')
    assert response.status_code in (301, 302)
    with client.session_transaction() as session:
        flashes = session.get('_flashes', [])
    assert any('Log-in or Sign-Up to use this feature.' in msg for _cat, msg in flashes)


def test_authenticated_favorites_crud(auth_client):
    # Add a favorite
    r = auth_client.post("/api/favorite", json={"locusId": "locus-smoke-1", "fav": True, "groupId": None})
    assert r.status_code == 200, r.get_data(as_text=True)
    assert r.get_json()["status"] == "ok"

    # Check status
    r = auth_client.get("/api/favorite", query_string={"locusId": "locus-smoke-1"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is True

    # List all favorites
    r = auth_client.get("/api/favorites")
    assert r.status_code == 200
    ids = [f["locusId"] for f in r.get_json()["favorites"]]
    assert "locus-smoke-1" in ids

    # Create a group
    r = auth_client.post("/api/favorite-groups", json={"name": "Smoke Group"})
    assert r.status_code == 201
    group_id = r.get_json()["id"]

    # List groups
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    group_names = [g["name"] for g in r.get_json()["groups"]]
    assert "Smoke Group" in group_names

    # Delete the group
    r = auth_client.delete(f"/api/favorite-groups/{group_id}")
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"

    # Remove the favorite
    r = auth_client.post("/api/favorite", json={"locusId": "locus-smoke-1", "fav": False, "groupId": None})
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"

    r = auth_client.get("/api/favorite", query_string={"locusId": "locus-smoke-1"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is False


def test_authenticated_visual_query(auth_client):
    # The UI requires at least one valid rule before preview/save.
    minimal_rules = {
        "condition": "AND",
        "rules": [
            {
                "field": "featuretable.alert_id",
                "operator": "is_not_null",
            }
        ],
    }

    # Preview query
    r = auth_client.post("/api/preview-query", json={"rules": minimal_rules})
    assert r.status_code == 200
    assert "sql" in r.get_json()

    # Export query (match count)
    r = auth_client.post("/api/export-query", json={"rules": minimal_rules})
    assert r.status_code == 200
    assert "count" in r.get_json()

    # Visual query page renders
    r = auth_client.get("/visual_query")
    assert r.status_code == 200


def test_authenticated_visual_query_rejects_malformed_between(auth_client):
    malformed_rules = {
        "condition": "AND",
        "rules": [
            {
                "field": "featuretable.locus_ra",
                "operator": "between",
                "value": [118.0],
            }
        ],
    }

    r = auth_client.post("/api/preview-query", json={"rules": malformed_rules})
    assert r.status_code == 400
    assert r.is_json
    assert "error" in r.get_json()


def test_authenticated_visual_query_rejects_malformed_in(auth_client):
    malformed_rules = {
        "condition": "AND",
        "rules": [
            {
                "field": "featuretable.alert_id",
                "operator": "in",
                "value": "ztf_candidate:123",
            }
        ],
    }

    r = auth_client.post("/api/preview-query", json={"rules": malformed_rules})
    assert r.status_code == 400
    assert r.is_json
    assert "error" in r.get_json()


def test_authenticated_visual_query_rejects_malformed_not_in(auth_client):
    malformed_rules = {
        "condition": "AND",
        "rules": [
            {
                "field": "featuretable.alert_id",
                "operator": "not_in",
                "value": "ztf_candidate:123",
            }
        ],
    }

    r = auth_client.post("/api/preview-query", json={"rules": malformed_rules})
    assert r.status_code == 400
    assert r.is_json
    assert "error" in r.get_json()


def test_authenticated_watchlist_crud(auth_client):
    # Match current UI contract: watchlists are saved from non-empty rules only.
    minimal_rules = {
        "condition": "AND",
        "rules": [
            {
                "field": "featuretable.alert_id",
                "operator": "is_not_null",
            }
        ],
    }

    # Create watchlist
    r = auth_client.post("/api/watchlist", json={"name": "Smoke Watchlist", "rules": minimal_rules})
    assert r.status_code == 201, r.get_data(as_text=True)
    wl_id = r.get_json()["id"]

    # List watchlists
    r = auth_client.get("/api/watchlist")
    assert r.status_code == 200
    names = [w["name"] for w in r.get_json()]
    assert "Smoke Watchlist" in names

    # Delete watchlist
    r = auth_client.delete(f"/api/watchlist/{wl_id}")
    assert r.status_code == 204


def test_lightcurve_and_features_smoke(client):
    response = client.get("/query_lightcurve_data", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/locus_plot_csv", query_string={"locusId": "locus-1"})
    assert response.status_code == 200
    assert "text/csv" in response.content_type

    response = client.get("/query_features", query_string={"alert_id": "missing-alert"})
    assert response.status_code in (404, 500)

    response = client.get("/query_featureplot_data", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/query_crossmatches", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/download_alerts_csv")
    assert response.status_code == 400
