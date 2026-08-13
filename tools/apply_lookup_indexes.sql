-- One-shot lookup indexes for an existing alerts SQLite file.
-- Safe to re-run (IF NOT EXISTS). Do not use db.create_all() on a populated
-- production DB, and do not run tools/alertsdb_schema.sql CREATE TABLE against it.
--
-- SQLite will lock the file while each index builds (~minutes on ~1.4M wide
-- rows). Apply during a maintenance window, then ANALYZE.
--
--   sqlite3 /path/to/alerts.db < tools/apply_lookup_indexes.sql
--
-- Names match poi_broker.models and tools/alertsdb_schema.sql. A DESC index on
-- date_alert_mjd is omitted: SQLite uses this B-tree in both directions, and
-- idx_featuretable_date_alert_mjd already exists on the measured live file.
-- idx_classification_alert_id is omitted: alert_id is already the PK.

CREATE INDEX IF NOT EXISTS idx_featuretable_date_alert_mjd ON featuretable (date_alert_mjd);
CREATE INDEX IF NOT EXISTS idx_featuretable_alert_id ON featuretable (alert_id);
CREATE INDEX IF NOT EXISTS idx_featuretable_ztf_object_id ON featuretable (ztf_object_id);
CREATE INDEX IF NOT EXISTS idx_featuretable_locus_id ON featuretable (locus_id);
CREATE INDEX IF NOT EXISTS idx_featuretable_ant_passband ON featuretable (ant_passband);
CREATE INDEX IF NOT EXISTS idx_featuretable_locus_ra ON featuretable (locus_ra);
CREATE INDEX IF NOT EXISTS idx_featuretable_locus_dec ON featuretable (locus_dec);
CREATE INDEX IF NOT EXISTS idx_featuretable_ant_mag_corrected ON featuretable (ant_mag_corrected);
CREATE INDEX IF NOT EXISTS idx_crossmatches_locus_id ON crossmatches (locus_id);
CREATE INDEX IF NOT EXISTS idx_classification_prob_class ON classification (prob_class);

ANALYZE;
