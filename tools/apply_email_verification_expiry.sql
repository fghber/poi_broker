-- One-shot: add verification-token expiry on an existing users SQLite file.
-- Schema upgrades are manual: the application never alters the database
-- (see docs/password_reset/deployment.md). Do not re-run if the column
-- already exists (SQLite ADD COLUMN will fail).
--
--   sqlite3 /path/to/users.db < tools/apply_email_verification_expiry.sql

ALTER TABLE user ADD COLUMN email_verification_token_expires INTEGER;
