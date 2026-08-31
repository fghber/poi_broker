-- One-shot: add password-changed epoch on an existing users SQLite file.
-- Required for session invalidation (see docs/password_reset/deployment.md).
-- The application never alters the database itself. Do not re-run if the
-- column already exists (SQLite ADD COLUMN will fail).
--
--   sqlite3 /path/to/users.db < tools/apply_password_changed_at.sql

ALTER TABLE user ADD COLUMN password_changed_at INTEGER;
