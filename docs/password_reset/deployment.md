# Deployment: session invalidation (`user.password_changed_at`)

Password resets and password changes now end **every** existing session and
remembered login for that account. This is enforced by a new epoch-seconds
column, `password_changed_at`, on the `user` table of the **users database**
(`users.db` only — the alerts database is untouched).

There is no auto-migration: **the application never alters the database.**
Schema upgrades are applied manually with one-shot SQL scripts in `tools/`
(before v3.6 the boot-time helper `_ensure_email_verification_expiry_column`
auto-added one column; that helper is removed, hence the conditional step
below). Both changes are single-column `ALTER TABLE`s and are safe to apply
while the app is running.

## 1. Locate the users DB and check current columns

The path comes from `USERS_DB_PATH` in the app's environment (`.env`), which
defaults to a sibling `_broker_db/users.db` outside the repository.

```bash
sqlite3 "$USERS_DB_PATH" "PRAGMA table_info(user);"
```

Look at the `name` column of the output:

| Output contains          | Means                        | Action             |
| ------------------------ | ---------------------------- | ------------------ |
| `password_changed_at`    | Already migrated             | Skip step 2a       |
| `email_verification_token_expires` | Boot helper ran its ALTER on this site previously | Skip step 2b |
| Neither                  | Pre-v3.5 schema              | Run both 2a and 2b |

## 2. Apply the one-shot SQL script(s)

**Do not re-run a script against an already-migrated database** — SQLite's
`ADD COLUMN` fails if the column exists (see Troubleshooting).

### 2a. Required for this release

```bash
sqlite3 "$USERS_DB_PATH" < tools/apply_password_changed_at.sql
```

### 2b. Only when upgrading across versions that shipped the boot-time helper

```bash
sqlite3 "$USERS_DB_PATH" < tools/apply_email_verification_expiry.sql
```

## 3. Restart

Restart the Gunicorn service hosting the web app as usual
(`<your-poi-broker-web-unit>` is whatever unit runs `wsgi:app`; the Huey
worker needs no restart, but restarting it does no harm):

```bash
sudo systemctl restart <your-poi-broker-web-unit>
```

Apply the script(s) **before** restarting: once the app boots against the new
ORM, any authenticated request hits a SELECT referencing
`password_changed_at`, which errors out with `no such column` until the
column exists. The reverse order (script first) has no such window.

## What changes for users

- Nobody is logged out by the deployment itself. Accounts untouched by their
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

## Rollback

No rollback is required for a mistaken upgrade: the column is additive and
nullable, and old code ignores it. To fully revert, remove it with
`ALTER TABLE user DROP COLUMN password_changed_at;` (SQLite ≥ 3.35) and
restart on the previous release.

## Troubleshooting

### `duplicate column name: password_changed_at`

The script ran twice, or the column pre-exists. Nothing is damaged — verify
with `PRAGMA table_info(user);` and move on.

### `no such table: user`

Wrong database file: re-check `USERS_DB_PATH`. Note there are two SQLite
files in play (`ztf_alerts_stream.db`, `users.db`); both migration scripts
belong to the users one.
