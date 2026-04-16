import json

from flask import Blueprint, render_template, request, redirect, url_for, flash
from flask_login import login_required, current_user

from . import db
from .constants.features import FEATURE_COLUMNS, default_feature_plot_columns

user_settings_bp = Blueprint('user_settings', __name__)

_MAX_DEFAULT_FEATURE_PLOT_COLUMNS = 10

"""
CREATE TABLE user_settings (
    id INTEGER PRIMARY KEY,
    user_id INTEGER NOT NULL,
    default_feature_plot_columns TEXT,
    FOREIGN KEY(user_id) REFERENCES user(id) ON DELETE CASCADE,
    UNIQUE(user_id)
);
"""


class UserSettings(db.Model):
    __bind_key__ = 'users'
    __tablename__ = 'user_settings'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', ondelete='CASCADE'), nullable=False, index=True)
    default_feature_plot_columns = db.Column(db.Text, nullable=True)

    __table_args__ = (db.UniqueConstraint('user_id', name='uix_user_settings_user_id'),)

    def __repr__(self):
        return f"<UserSettings user_id={self.user_id}>"


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


def save_feature_plot_columns_for_user(user_id, feature_columns):
    settings = UserSettings.query.filter_by(user_id=user_id).first()
    if settings is None:
        settings = UserSettings(user_id=user_id)
    settings.default_feature_plot_columns = json.dumps(_normalize_feature_columns(feature_columns))
    db.session.add(settings)
    db.session.commit()
    return settings


@user_settings_bp.route('/settings', methods=['GET'])
@login_required
def settings():
    selected_features = get_saved_feature_plot_columns(current_user.id)
    return render_template(
        'user_settings.html',
        selected_features=selected_features,
        available_features=list(FEATURE_COLUMNS.keys()),
    )


def _submitted_feature_plot_columns():
    cols = request.form.getlist('default_feature_plot_columns[]')
    if cols:
        return cols
    return request.form.getlist('default_feature_plot_columns')


@user_settings_bp.route('/settings', methods=['POST'])
@login_required
def save_settings():
    selected_features = _submitted_feature_plot_columns()
    saved = save_feature_plot_columns_for_user(current_user.id, selected_features)
    flash('Your default feature plot columns have been saved.', 'success')
    return redirect(url_for('user_settings.settings'))
