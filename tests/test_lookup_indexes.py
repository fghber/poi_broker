"""Lookup indexes declared on alerts/users models must exist after create_all()."""

import sqlalchemy as sa

from poi_broker import db

EXPECTED_ALERTS_INDEXES = {
    "idx_featuretable_alert_id",
    "idx_featuretable_locus_id",
    "idx_featuretable_ztf_object_id",
    "idx_featuretable_ant_passband",
    "idx_featuretable_locus_ra",
    "idx_featuretable_locus_dec",
    "idx_featuretable_ant_mag_corrected",
    "idx_featuretable_date_alert_mjd",
    "idx_crossmatches_locus_id",
    "idx_classification_prob_class",
}

EXPECTED_USERS_INDEXES = {
    "uix_filter_bookmark_user_name",
    "idx_export_task_status_updated_at",
    "idx_export_task_status_created_at",
    "uix_export_task_one_active_per_user",
}


def test_alerts_lookup_indexes_exist_after_create_all(app):
    with app.app_context():
        rows = db.session.execute(
            sa.text("SELECT name FROM sqlite_master WHERE type='index'")
        ).fetchall()

    names = {row[0] for row in rows}
    missing = EXPECTED_ALERTS_INDEXES - names
    assert not missing, f"missing lookup indexes: {sorted(missing)}"


def test_users_housekeeping_indexes_exist_after_create_all(app):
    with app.app_context():
        rows = db.session.execute(
            sa.text("SELECT name FROM sqlite_master WHERE type='index'"),
            bind_arguments={"bind": db.engines["users"]},
        ).fetchall()

    names = {row[0] for row in rows}
    missing = EXPECTED_USERS_INDEXES - names
    assert not missing, f"missing users indexes: {sorted(missing)}"
