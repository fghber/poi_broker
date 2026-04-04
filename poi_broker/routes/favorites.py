"""Favorites API routes blueprint."""

import logging
from flask import Blueprint, jsonify, request
from flask_login import login_required, current_user
from ..services.favorites_service import (
    get_favorite_status,
    get_user_favorites,
    toggle_favorite,
    update_favorite_group,
    get_favorite_groups,
    create_favorite_group,
    delete_favorite_group,
)

logger = logging.getLogger(__name__)

favorites_bp = Blueprint('favorites', __name__, url_prefix='/api')


@favorites_bp.route('/favorite', methods=['GET'])
@login_required
def api_favorite_get():
    """
    GET /api/favorite?locusId=...  -> returns {"fav": true/false}
    """
    locus_id = request.args.get('locusId')
    fav = get_favorite_status(locus_id)
    return jsonify({'fav': fav})


@favorites_bp.route('/favorites', methods=['GET'])
@login_required
def api_favorites_get():
    """GET /api/favorites?groupId=<id> -> return favorites with IDs, optionally filtered by group."""
    group_id = request.args.get('groupId', type=int)
    
    # If groupId=null explicitly, show only ungrouped
    if 'groupId' in request.args and request.args.get('groupId') == 'null':
        group_id = None
    elif group_id is None and 'groupId' not in request.args:
        # No group filter specified, show all
        group_id = 'all'
    
    if group_id == 'all':
        from .. import db
        from ..models import Favorite
        favs = [{'id': r.id, 'locusId': r.locus_id} for r in Favorite.query.filter_by(user_id=current_user.id).all()]
    else:
        favs = get_user_favorites(group_id)
    
    return jsonify({'favorites': favs})


@favorites_bp.route('/favorite', methods=['POST'])
@login_required
def api_favorite_post():
    """
    POST /api/favorite
    JSON body: { "locusId": "...", "fav": true, "groupId": null }
    """
    data = request.get_json(silent=True)
    if data is None:
        return jsonify({'status': 'Invalid or missing JSON body'}), 400
    
    locus_id = data.get('locusId')
    if not locus_id:
        return jsonify({'status': 'Missing locusId'}), 400
    
    fav_flag = bool(data.get('fav'))
    group_id = data.get('groupId')
    
    result, status_code = toggle_favorite(locus_id, fav_flag, group_id)
    return jsonify(result), status_code


@favorites_bp.route('/favorite/<int:favorite_id>/group', methods=['PATCH'])
@login_required
def api_favorite_update_group(favorite_id):
    """PATCH /api/favorite/<id>/group -> move favorite to a group. Body: {"groupId": <id> or null}"""
    data = request.get_json(silent=True)
    if data is None:
        return jsonify({'error': 'Invalid JSON'}), 400
    
    group_id = data.get('groupId')
    result, status_code = update_favorite_group(favorite_id, group_id)
    return jsonify(result), status_code


@favorites_bp.route('/favorite-groups', methods=['GET'])
@login_required
def api_favorite_groups_get():
    """GET /api/favorite-groups -> return all user's groups with favorites counts."""
    groups = get_favorite_groups()
    return jsonify({'groups': groups})


@favorites_bp.route('/favorite-groups', methods=['POST'])
@login_required
def api_favorite_groups_post():
    """POST /api/favorite-groups -> create a new group. Body: {"name": "Group A"}"""
    data = request.get_json(silent=True)
    if not data or not data.get('name'):
        return jsonify({'error': 'name required'}), 400
    
    result, status_code = create_favorite_group(data['name'])
    return jsonify(result), status_code


@favorites_bp.route('/favorite-groups/<int:group_id>', methods=['DELETE'])
@login_required
def api_favorite_groups_delete(group_id):
    """DELETE /api/favorite-groups/<id> -> delete a group (orphans its favorites)."""
    result, status_code = delete_favorite_group(group_id)
    return jsonify(result), status_code
