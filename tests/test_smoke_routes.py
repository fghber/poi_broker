"""
Smoke test coverage:

Public pages: /, /help, /contact
Auth-required routes (unauth behavior)
Rate limiting
Test run result:

Execute command via workspace Python env: `pytest -q` or `python -m pytest -q`
"""

from poi_broker.services.catalog_list import PAGE_SIZE

def test_public_pages_smoke(client):
    for path in ["/", "/help", "/contact"]:
        response = client.get(path)
        assert response.status_code == 200, f"Expected 200 for {path}, got {response.status_code}"


def test_catalog_serves_same_origin_bokeh(client):
    page = client.get("/")
    assert page.status_code == 200
    assert b"cdn.bokeh.org" not in page.data
    assert b"/bokeh.min.js" in page.data
    assert "cdn.bokeh.org" not in (page.headers.get("Content-Security-Policy") or "")

    js = client.get("/bokeh.min.js")
    assert js.status_code == 200
    assert js.content_type and "javascript" in js.content_type
    assert len(js.data) > 1000


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


def test_filter_bookmarks_routes_require_login(client):
    response = client.get("/api/filter-bookmarks")
    assert response.status_code == 401
    assert response.is_json

    response = client.post("/api/filter-bookmarks", json={"name": "x", "params": {}})
    assert response.status_code == 401

    response = client.delete("/api/filter-bookmarks/1")
    assert response.status_code == 401


def test_user_observatories_routes_require_login(client):
    response = client.get("/api/user-observatories")
    assert response.status_code == 401
    assert response.is_json

    response = client.post("/api/user-observatories", json={"name": "x", "latitude": 0.0, "longitude": 0.0})
    assert response.status_code == 401

    response = client.delete("/api/user-observatories/1")
    assert response.status_code == 401

    response = client.post("/api/last-observatory", json={"source": "builtin", "name": "Palomar"})
    assert response.status_code == 401
    assert response.is_json


def test_ui_protected_route_sets_flash_message(client):
    response = client.get('/visual_query')
    assert response.status_code in (301, 302)
    with client.session_transaction() as session:
        flashes = session.get('_flashes', [])
    assert any('Log-in or Sign-Up to use this feature.' in msg for _cat, msg in flashes)


def test_lightcurve_and_features_smoke(client):
    response = client.get("/query_lightcurve_data", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/locus_plot_csv", query_string={"locusId": "locus-1"})
    assert response.status_code == 200
    assert "text/csv" in response.content_type

    response = client.get("/query_features", query_string={"alert_id": "missing-alert"})
    assert response.status_code in (404, 500)
    assert response.is_json

    response = client.get("/query_featureplot_data", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/query_crossmatches", query_string={"locusId": "locus-1"})
    assert response.status_code == 200

    response = client.get("/download_alerts_csv")
    assert response.status_code == 400


def test_download_alerts_csv_caps_alert_id_count(client):
    """The endpoint serves one UI page of rows; reject larger lists server-side."""
    too_many = {"alert_id": [f"alert-{i}" for i in range(PAGE_SIZE + 1)]}
    response = client.get("/download_alerts_csv", query_string=too_many)
    assert response.status_code == 400
    assert str(PAGE_SIZE) in response.get_data(as_text=True)


def test_query_features_missing_alert_id_returns_json_error(client):
    response = client.get("/query_features")
    assert response.status_code == 400
    assert response.is_json
    assert response.get_json()["error"] == "Missing alert_id"


def test_features_routes_reject_oversized_ids(client):
    """alert_id/locusId enforce the same 128-char cap as favorites and classification."""
    long_id = "a" * 129

    response = client.get("/query_features", query_string={"alert_id": long_id})
    assert response.status_code == 400
    assert response.get_json()["error"] == "alert_id is too long"

    response = client.get("/query_featureplot_data", query_string={"locusId": long_id})
    assert response.status_code == 400
    assert response.get_json()["error"] == "locusId is too long"


def test_read_routes_rate_limited(secure_client):
    """Test that heavy read routes are rate limited."""
    # Test LAX limit (30 per minute) for main page
    for i in range(31):
        response = secure_client.get("/")
        if i < 30:
            assert response.status_code == 200
        else:
            assert response.status_code == 429  # Rate limited

    # Test LAX limit for /query_features
    for i in range(31):
        response = secure_client.get("/query_features", query_string={"alert_id": "missing-alert"})
        if i < 30:
            assert response.status_code in (404, 500)
        else:
            assert response.status_code == 429

    # Test LAX limit for /query_featureplot_data
    for i in range(31):
        response = secure_client.get("/query_featureplot_data", query_string={"locusId": "locus-1"})
        if i < 30:
            assert response.status_code == 200
        else:
            assert response.status_code == 429

    # Test LAX limit for /query_lightcurve_data
    for i in range(31):
        response = secure_client.get("/query_lightcurve_data", query_string={"locusId": "locus-1"})
        if i < 30:
            assert response.status_code == 200
        else:
            assert response.status_code == 429

    # Test LAX limit for /locus_plot_csv
    for i in range(31):
        response = secure_client.get("/locus_plot_csv", query_string={"locusId": "locus-1"})
        if i < 30:
            assert response.status_code == 200
        else:
            assert response.status_code == 429

    # Test LAX limit for /query_crossmatches
    for i in range(31):
        response = secure_client.get("/query_crossmatches", query_string={"locusId": "locus-1"})
        if i < 30:
            assert response.status_code == 200
        else:
            assert response.status_code == 429

    # Test MEDIUM limit (15 per minute) for /download_alerts_csv
    for i in range(16):
        response = secure_client.get("/download_alerts_csv")
        if i < 15:
            assert response.status_code == 400
        else:
            assert response.status_code == 429

    # Cheap 400 path only — do not render matplotlib plots in this loop.
    for i in range(31):
        response = secure_client.get("/query_observing_plot")
        if i < 30:
            assert response.status_code == 400
        else:
            assert response.status_code == 429
