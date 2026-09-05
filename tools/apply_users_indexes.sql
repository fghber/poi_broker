-- One-shot indexes / constraints for an existing users SQLite file.
-- Safe to re-run (IF NOT EXISTS) when there are no duplicate bookmark names.
-- Do not use db.create_all() on a populated production DB, and do not re-run
-- tools/usersdb_schema.sql CREATE TABLE against it.
--
--   sqlite3 /path/to/users.db < tools/apply_users_indexes.sql
--
-- Names match poi_broker.models and tools/usersdb_schema.sql.
-- Non-unique export indexes are created first so a duplicate-bookmark abort
-- does not skip them.

-- O6: ExportTask housekeeping composites
CREATE INDEX IF NOT EXISTS idx_export_task_status_updated_at
    ON export_task (status, updated_at);
CREATE INDEX IF NOT EXISTS idx_export_task_status_created_at
    ON export_task (status, created_at);

-- O5: abort if duplicate (user_id, name) rows would make the unique index fail.
-- sqlite3 prints any duplicate groups; the CHECK insert then stops the script.
SELECT user_id, name, COUNT(*) AS duplicate_count
FROM filter_bookmark
GROUP BY user_id, name
HAVING COUNT(*) > 1;

CREATE TEMP TABLE _filter_bookmark_unique_guard (
    ok INTEGER NOT NULL CHECK (ok = 1)
);
INSERT INTO _filter_bookmark_unique_guard (ok)
SELECT CASE
    WHEN EXISTS (
        SELECT 1
        FROM filter_bookmark
        GROUP BY user_id, name
        HAVING COUNT(*) > 1
    ) THEN 0
    ELSE 1
END;
DROP TABLE _filter_bookmark_unique_guard;

CREATE UNIQUE INDEX IF NOT EXISTS uix_filter_bookmark_user_name
    ON filter_bookmark (user_id, name);

ANALYZE;
