"""Visual query builder routes blueprint."""

import logging
import json
import time
from flask import Blueprint, jsonify, request, render_template
from flask_login import login_required, current_user
from .. import db
from ..models import Watchlist
from ..services.query_service import get_preview_sql, get_query_match_count, build_query_from_rules

logger = logging.getLogger(__name__)

visual_query_bp = Blueprint('visual_query', __name__)


@visual_query_bp.route('/visual_query')
@login_required
def visual_query():
    """Show visual query interface."""
    return render_template('visual_query.html')


@visual_query_bp.route('/api/preview-query', methods=['POST'])
@login_required
def preview_query():
    """Get SQL preview for a query builder payload."""
    try:
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            return jsonify({'error': 'Invalid or missing JSON'}), 400

        # Accept either direct rules object or wrapped payload {'rules': {...}}.
        rules_payload = data
        if isinstance(data.get('rules'), dict) and 'condition' in data.get('rules', {}):
            rules_payload = data['rules']

        sql_preview = get_preview_sql(rules_payload)
        return jsonify({'sql': sql_preview})
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        logger.error('Error in preview_query', exc_info=True)
        return jsonify({'error': str(e)}), 500


@visual_query_bp.route('/api/export-query', methods=['POST'])
@login_required
def export_query():
    """Get match count for a query builder payload."""
    try:
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            return jsonify({'error': 'Invalid or missing JSON'}), 400

        # Accept either direct rules object or wrapped payload {'rules': {...}}.
        rules_payload = data
        if isinstance(data.get('rules'), dict) and 'condition' in data.get('rules', {}):
            rules_payload = data['rules']

        count = get_query_match_count(rules_payload)
        return jsonify({'count': count})
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        logger.error('Error in export_query', exc_info=True)
        return jsonify({'error': str(e)}), 500


@visual_query_bp.route('/api/watchlist', methods=['POST'])
@login_required
def save_watchlist():
    """Save a watchlist from query builder rules."""
    try:
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            return jsonify({'error': 'Invalid or missing JSON'}), 400

        name = data.get('name', '').strip()
        if not name:
            return jsonify({'error': 'Watchlist name is required'}), 400
        if len(name) > 128:
            return jsonify({'error': 'Watchlist name must be 128 characters or fewer'}), 400

        rules_payload = data.get('rules')
        if isinstance(rules_payload, dict) and 'condition' in rules_payload:
            pass
        elif isinstance(data.get('rules'), dict):
            rules_payload = data['rules']
        else:
            return jsonify({'error': 'Invalid querybuilder rules payload'}), 400

        # Build query and get SQL
        _, sql_where = build_query_from_rules(rules_payload)

        now_epoch = int(time.time())

        watchlist = Watchlist(
            user_id=current_user.id,
            name=name,
            rules_json=json.dumps(rules_payload),
            sql_where=sql_where,
            created_at=now_epoch,
        )
        db.session.add(watchlist)
        try:
            db.session.commit()
        except Exception:
            db.session.rollback()
            return jsonify({'error': f"A watchlist named '{name}' already exists"}), 409

        return jsonify({'id': watchlist.id, 'name': watchlist.name}), 201
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        logger.error('Error in save_watchlist', exc_info=True)
        return jsonify({'error': str(e)}), 500


@visual_query_bp.route('/api/watchlist', methods=['GET'])
@login_required
def list_watchlists():
    """List all watchlists for the current user."""
    try:
        rows = Watchlist.query.filter_by(user_id=current_user.id).order_by(Watchlist.created_at.desc()).all()
        result = [
            {'id': w.id, 'name': w.name, 'sql_where': w.sql_where, 'created_at': w.created_at}
            for w in rows
        ]
        return jsonify(result)
    except Exception as e:
        logger.error('Error in list_watchlists', exc_info=True)
        return jsonify({'error': str(e)}), 500


@visual_query_bp.route('/api/watchlist/<int:wid>', methods=['DELETE'])
@login_required
def delete_watchlist(wid):
    """Delete a watchlist by ID."""
    try:
        watchlist = Watchlist.query.filter_by(id=wid, user_id=current_user.id).first()
        if watchlist is None:
            return jsonify({'error': 'Not found'}), 404
        db.session.delete(watchlist)
        db.session.commit()
        return '', 204
    except Exception as e:
        logger.error('Error in delete_watchlist', exc_info=True)
        return jsonify({'error': str(e)}), 500
