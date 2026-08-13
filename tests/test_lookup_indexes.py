"""Lookup indexes declared on alerts models must exist after create_all()."""

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


def test_alerts_lookup_indexes_exist_after_create_all(app):
    with app.app_context():
        rows = db.session.execute(
            sa.text("SELECT name FROM sqlite_master WHERE type='index'")
        ).fetchall()

    names = {row[0] for row in rows}
    missing = EXPECTED_ALERTS_INDEXES - names
    assert not missing, f"missing lookup indexes: {sorted(missing)}"
