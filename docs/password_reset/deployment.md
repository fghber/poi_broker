# Deployment: session invalidation (`user.password_changed_at`)

> **Historical record — this migration has been applied.** Both one-shot
> scripts referenced below were applied to the production `users.db` and
> removed from the repository in 2026-10
> (`docs/adr/adr_db_schema_migration_strategy.md`). Fresh installations get
> both columns from `tools/usersdb_schema.sql`. Do not try to run the
> scripts; keep this page as the deployment record.

Password resets and password changes now end **every** existing session and
remembered login for that account. This is enforced by a new epoch-seconds
column, `password_changed_at`, on the `user` table of the **users database**
(`users.db` only — the alerts database is untouched), plus
`email_verification_token_expires` which replaced a legacy boot-time ALTER
helper.

There was no auto-migration: **the application never alters the database.**
The two schema upgrades were applied manually with one-shot SQL scripts
(`tools/apply_password_changed_at.sql`, `tools/apply_email_verification_expiry.sql`),
run before restarting the web service — once the app boots against the new
ORM, any authenticated request hits a SELECT referencing
`password_changed_at`, which errors out with `no such column` until the
column exists.

## What changed for users

- Nobody was logged out by the deployment itself. Accounts untouched by their
  owner (`password_changed_at` still NULL) behave exactly as before,
  remember-me included.
- The first time each user changes their password (or uses a reset link),
  every session and remember-me login for the account ends — including the
  browser performing the change. They sign back in with the new password;
  its cookie was cleared during logout.
- Remember-me cookies left behind on *other* devices remain in those browsers
  but are inert: the server rejects them because they carry no valid login
  timestamp against the newer `password_changed_at`. Browsers drop them at
  cookie expiry or the next sign-in on that device.

## Rollback (if ever needed)

The column is additive and nullable, and old code ignores it. To fully revert,
remove it with `ALTER TABLE user DROP COLUMN password_changed_at;`
(SQLite ≥ 3.35) and restart on the previous release.
