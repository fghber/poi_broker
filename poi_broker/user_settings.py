import json
from flask import Blueprint, render_template, request, redirect, url_for, flash
from flask_login import login_required, current_user
import logging
from sqlalchemy.exc import IntegrityError
from . import db
from .constants.features import FEATURE_COLUMNS, default_feature_plot_columns
from .models import UserSettings

logger = logging.getLogger(__name__)
user_settings_bp = Blueprint('user_settings', __name__)

_MAX_DEFAULT_FEATURE_PLOT_COLUMNS = 10


def _normalize_feature_columns(columns):
    if not isinstance(columns, list):
        return default_feature_plot_columns()

    normalized = []
    for col in columns:
        if isinstance(col, str) and col in FEATURE_COLUMNS and col not in normalized:
            normalized.append(col)
        if len(normalized) >= _MAX_DEFAULT_FEATURE_PLOT_COLUMNS:
            break

    return normalized if normalized else default_feature_plot_columns()


def get_saved_feature_plot_columns(user_id):
    settings = UserSettings.query.filter_by(user_id=user_id).first()
    if not settings or not settings.default_feature_plot_columns:
        return default_feature_plot_columns()

    try:
        columns = json.loads(settings.default_feature_plot_columns)
    except (TypeError, ValueError):
        return default_feature_plot_columns()

    return _normalize_feature_columns(columns)

@user_settings_bp.route('/settings', methods=['GET'])
@login_required
def settings():
    selected_features = get_saved_feature_plot_columns(current_user.id)
    return render_template(
        'user_settings.html',
        selected_features=selected_features,
        available_features=list(FEATURE_COLUMNS.keys()),
    )

@user_settings_bp.route('/settings', methods=['POST'])
@login_required
def save_settings():
    # Accept both naming conventions used by forms/clients:
    # - name="default_feature_plot_columns[]"
    # - name="default_feature_plot_columns"
    selected_features = request.form.getlist('default_feature_plot_columns')
    if not selected_features:
        selected_features = request.form.getlist('default_feature_plot_columns[]')

    if not selected_features:
        flash('No features selected. Please select at least one feature.', 'error')
        return redirect(url_for('user_settings.settings'))

    settings = UserSettings.query.filter_by(user_id=current_user.id).first()
    if settings is None:
        settings = UserSettings(user_id=current_user.id)
    settings.default_feature_plot_columns = json.dumps(_normalize_feature_columns(selected_features))
    db.session.add(settings)

    try:
        db.session.commit()
    except IntegrityError:
        db.session.rollback()
        flash('Failed to save feature plot defaults. Please try if you have any unsaved changes.', 'warning')
        return redirect(url_for('user_settings.settings'))
    except Exception as e:
        db.session.rollback()
        logger.error(f'Database error during commit: {str(e)}', exc_info=True)
        flash('Failed to save feature plot defaults. Please try again later.', 'danger')
        return redirect(url_for('user_settings.settings'))

    flash('Your default feature plot columns have been saved.', 'success')
    return redirect(url_for('user_settings.settings'))