"""Tests for favorites API."""

import json

from werkzeug.security import generate_password_hash


def test_favorites_unauthenticated(client):
    """Test that favorites endpoints require authentication."""
    r = client.get("/api/favorite", query_string={"locusId": "test"})
    assert r.status_code == 401
    assert r.is_json
    assert r.get_json().get("error") == "authentication required"

    r = client.get("/api/favorites")
    assert r.status_code == 401

    r = client.post("/api/favorite", json={"locusId": "test", "fav": True})
    assert r.status_code == 401

    r = client.patch("/api/favorite/1/group", json={"groupId": None})
    assert r.status_code == 401

    r = client.get("/api/favorite-groups")
    assert r.status_code == 401

    r = client.post("/api/favorite-groups", json={"name": "Test Group"})
    assert r.status_code == 401

    r = client.delete("/api/favorite-groups/1")
    assert r.status_code == 401


def test_favorites_crud(auth_client):
    """Test basic favorites CRUD operations."""
    # Initially not favorited
    r = auth_client.get("/api/favorite", query_string={"locusId": "test-locus-1"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is False

    # Add favorite
    r = auth_client.post("/api/favorite", json={"locusId": "test-locus-1", "fav": True})
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"

    # Now should be favorited
    r = auth_client.get("/api/favorite", query_string={"locusId": "test-locus-1"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is True

    # List all favorites
    r = auth_client.get("/api/favorites")
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    assert isinstance(favs, list)
    assert len(favs) == 1
    assert favs[0]["locusId"] == "test-locus-1"

    # Add another favorite
    r = auth_client.post("/api/favorite", json={"locusId": "test-locus-2", "fav": True})
    assert r.status_code == 200

    # List should have 2
    r = auth_client.get("/api/favorites")
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    assert len(favs) == 2
    locus_ids = [f["locusId"] for f in favs]
    assert "test-locus-1" in locus_ids
    assert "test-locus-2" in locus_ids

    # Remove first favorite
    r = auth_client.post("/api/favorite", json={"locusId": "test-locus-1", "fav": False})
    assert r.status_code == 200

    # Should not be favorited anymore
    r = auth_client.get("/api/favorite", query_string={"locusId": "test-locus-1"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is False

    # List should have 1
    r = auth_client.get("/api/favorites")
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    assert len(favs) == 1
    assert favs[0]["locusId"] == "test-locus-2"

    # Clean up
    r = auth_client.post("/api/favorite", json={"locusId": "test-locus-2", "fav": False})
    assert r.status_code == 200


def test_favorites_groups_crud(auth_client):
    """Test favorite groups CRUD operations."""
    # Create a group
    r = auth_client.post("/api/favorite-groups", json={"name": "Test Group"})
    assert r.status_code == 201
    group_data = r.get_json()
    assert "id" in group_data
    assert group_data["name"] == "Test Group"
    group_id = group_data["id"]

    # List groups - should have Ungrouped + Test Group
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    assert isinstance(groups, list)
    assert len(groups) == 2
    group_names = [g["name"] for g in groups]
    assert "Test Group" in group_names
    assert "Ungrouped" in group_names

    # All should have count 0
    for g in groups:
        assert g["count"] == 0

    # Add favorite to group
    r = auth_client.post("/api/favorite", json={"locusId": "grouped-fav", "fav": True, "groupId": group_id})
    assert r.status_code == 200

    # Group count should update
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    test_group = next(g for g in groups if g["name"] == "Test Group")
    assert test_group["count"] == 1

    # List favorites in group
    r = auth_client.get("/api/favorites", query_string={"groupId": group_id})
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    assert len(favs) == 1
    assert favs[0]["locusId"] == "grouped-fav"

    # Delete group
    r = auth_client.delete(f"/api/favorite-groups/{group_id}")
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"

    # Group should be gone
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    group_names = [g["name"] for g in groups]
    assert "Test Group" not in group_names
    assert "Ungrouped" in group_names

    # Favorite should be orphaned to Ungrouped
    ungrouped = next(g for g in groups if g["name"] == "Ungrouped")
    assert ungrouped["count"] == 1

    # Clean up
    r = auth_client.post("/api/favorite", json={"locusId": "grouped-fav", "fav": False})
    assert r.status_code == 200


def test_favorites_validation(auth_client):
    """Test favorites input validation."""
    # Missing locusId
    r = auth_client.post("/api/favorite", json={"fav": True})
    assert r.status_code == 400
    assert "Missing locusId" in r.get_json()["status"]

    # Invalid JSON
    r = auth_client.post("/api/favorite", data="invalid json")
    assert r.status_code == 400
    assert "Invalid or missing JSON body" in r.get_json()["status"]

    # Empty group name
    r = auth_client.post("/api/favorite-groups", json={"name": ""})
    assert r.status_code == 400
    assert "name required" in r.get_json()["error"]

    # Missing name
    r = auth_client.post("/api/favorite-groups", json={})
    assert r.status_code == 400
    assert "name required" in r.get_json()["error"]

    # Invalid favorite ID for group update
    r = auth_client.patch("/api/favorite/99999/group", json={"groupId": None})
    assert r.status_code == 404
    assert "favorite not found" in r.get_json()["error"]

    # Invalid group ID for group update
    # First create a favorite
    r = auth_client.post("/api/favorite", json={"locusId": "validation-test", "fav": True})
    assert r.status_code == 200

    # Get the favorite ID
    r = auth_client.get("/api/favorites")
    fav_id = r.get_json()["favorites"][0]["id"]

    # Try to move to non-existent group
    r = auth_client.patch(f"/api/favorite/{fav_id}/group", json={"groupId": 99999})
    assert r.status_code == 404
    assert "group not found" in r.get_json()["error"]

    # Clean up
    r = auth_client.post("/api/favorite", json={"locusId": "validation-test", "fav": False})
    assert r.status_code == 200


def test_favorites_cross_user_isolation(app, auth_client):
    """Test that users cannot access other users' favorites."""
    from poi_broker import db
    from poi_broker.models import Favorite, FavoriteGroup, User

    with app.app_context():
        # Create another user
        other_user = User(
            email="other_fav@example.com",
            password=generate_password_hash("OtherUser123!"),
            name="Other User",
            email_verified=True,
        )
        db.session.add(other_user)
        db.session.commit()

        # Create favorite and group for other user
        other_fav = Favorite(
            user_id=other_user.id,
            locus_id="other-locus",
        )
        db.session.add(other_fav)

        other_group = FavoriteGroup(
            user_id=other_user.id,
            name="Other Group",
        )
        db.session.add(other_group)
        db.session.commit()

        other_fav_id = other_fav.id
        other_group_id = other_group.id

    # Current user should not see other's favorite
    r = auth_client.get("/api/favorite", query_string={"locusId": "other-locus"})
    assert r.status_code == 200
    assert r.get_json()["fav"] is False

    # Should not see in favorites list
    r = auth_client.get("/api/favorites")
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    locus_ids = [f["locusId"] for f in favs]
    assert "other-locus" not in locus_ids

    # Should not see other's group
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    group_names = [g["name"] for g in groups]
    assert "Other Group" not in group_names

    # Should not be able to delete other's favorite (but we don't have direct delete, so test group delete)
    r = auth_client.delete(f"/api/favorite-groups/{other_group_id}")
    assert r.status_code == 404
    assert "group not found" in r.get_json()["error"]

    # Verify other user's data still exists
    with app.app_context():
        still_fav = db.session.get(Favorite, other_fav_id)
        assert still_fav is not None
        still_group = db.session.get(FavoriteGroup, other_group_id)
        assert still_group is not None


def test_favorites_group_operations(auth_client):
    """Test advanced group operations like moving favorites."""
    # Create two groups
    r = auth_client.post("/api/favorite-groups", json={"name": "Group A"})
    assert r.status_code == 201
    group_a_id = r.get_json()["id"]

    r = auth_client.post("/api/favorite-groups", json={"name": "Group B"})
    assert r.status_code == 201
    group_b_id = r.get_json()["id"]

    # Add favorites to Group A
    r = auth_client.post("/api/favorite", json={"locusId": "move-test-1", "fav": True, "groupId": group_a_id})
    assert r.status_code == 200

    r = auth_client.post("/api/favorite", json={"locusId": "move-test-2", "fav": True, "groupId": group_a_id})
    assert r.status_code == 200

    # Check counts
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    group_a = next(g for g in groups if g["id"] == group_a_id)
    group_b = next(g for g in groups if g["id"] == group_b_id)
    assert group_a["count"] == 2
    assert group_b["count"] == 0

    # Get favorite ID to move
    r = auth_client.get("/api/favorites", query_string={"groupId": group_a_id})
    assert r.status_code == 200
    favs = r.get_json()["favorites"]
    fav_id = favs[0]["id"]

    # Move favorite from Group A to Group B
    r = auth_client.patch(f"/api/favorite/{fav_id}/group", json={"groupId": group_b_id})
    assert r.status_code == 200

    # Check updated counts
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    group_a = next(g for g in groups if g["id"] == group_a_id)
    group_b = next(g for g in groups if g["id"] == group_b_id)
    assert group_a["count"] == 1
    assert group_b["count"] == 1

    # Move back to ungrouped
    r = auth_client.patch(f"/api/favorite/{fav_id}/group", json={"groupId": None})
    assert r.status_code == 200

    # Check final counts
    r = auth_client.get("/api/favorite-groups")
    assert r.status_code == 200
    groups = r.get_json()["groups"]
    group_a = next(g for g in groups if g["id"] == group_a_id)
    group_b = next(g for g in groups if g["id"] == group_b_id)
    ungrouped = next(g for g in groups if g["name"] == "Ungrouped")
    assert group_a["count"] == 1
    assert group_b["count"] == 0
    assert ungrouped["count"] == 1

    # Clean up
    r = auth_client.post("/api/favorite", json={"locusId": "move-test-1", "fav": False})
    assert r.status_code == 200
    r = auth_client.post("/api/favorite", json={"locusId": "move-test-2", "fav": False})
    assert r.status_code == 200
    r = auth_client.delete(f"/api/favorite-groups/{group_a_id}")
    assert r.status_code == 200
    r = auth_client.delete(f"/api/favorite-groups/{group_b_id}")
    assert r.status_code == 200