# ADR: Epoch-based session invalidation on password change/reset

- **Status:** Accepted
- **Date:** 2026-08
- **Feature:** Auth session lifecycle (`poi_broker/auth.py`, `user_loader`)

## Context

Flask-Login sessions are stateless signed cookies validated by the
`user_loader` callback; the app keeps no server-side session store. Before
this change, neither a password reset (`POST /reset-password/<token>`) nor a
password change (`POST /change-password`) touched login state, so a session
cookie stolen before the credential change kept full access after the victim
rotated their password (insufficient session expiration on credential
change; remember-me cookies likewise survived, since flask-login restores
them through the same `user_loader`).

## Decision

Add an epoch-seconds watermark column `user.password_changed_at` (nullable,
users DB) and a `login_at` epoch marker inside the signed session, written by
`POST /login`. The `user_loader` — the single choke point through which both
session cookies and remember-cookie restores authenticate — rejects the
request unless the session carries an integer `login_at > password_changed_at`.
Both password-write routes stamp the watermark and call `logout_user()`,
which also expires the acting browser's remember-me cookie.

Epoch-second stamps cannot order events inside the watermark second, so that
second belongs to the old credential, hence the strict comparison. A
successful `POST /login` necessarily presented the current credential, so the
route stamps `login_at` just past the watermark
(`max(now, password_changed_at + 1)`) rather than bouncing a legitimate login
made in the same second as the write.

One uniform rule results: **any password write ends every session and
remembered login for that account, including the browser performing the
change.** `password_changed_at IS NULL` (account never rotated its password)
skips the check entirely, so deployment logs nobody out and pre-existing
remembered logins keep working until that account's first rotation.

## Alternatives considered

### Server-side session store (Flask-Session + SQLAlchemy/Redis) (rejected)

Full session enumeration and per-session revocation, but a new storage
dependency and a store read on every request. Same footprint logic as the
Huey ADR (`docs/async_export/adr_huey_sqlite.md`): no new networked service
for a small team.

### Flask-Login freshness / `SESSION_PROTECTION='strong'` (rejected)

Freshness marks (`_fresh`) say when the login happened, not whether
credentials changed since; "strong" protection keys on request fingerprints
(UA/IP), which is orthogonal. Neither mechanism binds sessions to credential
changes.

### Rotating the user's identity on password change (rejected)

Bumping `user.id` (or a token the remember cookie must match) would
invalidate sessions only via the same loader check we built anyway — while
`user.id` is referenced by half the users-DB foreign keys and embedded in the
remember cookie. The watermark column delivers the identical effect with no
FK churn.

## Consequences

- Stolen cookies cannot outlive a password reset/change; there is nothing
  client-side to keep them valid.
- The acting browser is signed out on password change and redirected to
  login (flash explains why). Accepted UX cost of the uniform rule.
- For accounts that have rotated their password at least once, remembered
  logins stop resurrecting sessions once the browser-side session expires
  (default 24h): a bare remember token carries no `login_at` and proves
  nothing about credential freshness. Deliberate, documented tradeoff.
- Deploy order matters: apply `tools/apply_password_changed_at.sql` **before**
  restarting into the new ORM, or authenticated pages 500 on the missing
  column. Runbook: `docs/password_reset/deployment.md`.
- The application never writes DDL. This ADR also removed the legacy
  boot-time helper (`_ensure_email_verification_expiry_column`); all schema
  upgrades are manual one-shot SQL scripts.

## Invariants (for implementers and coding agents)

1. Every code path that writes `user.password` MUST also stamp
   `user.password_changed_at` (currently `reset_password_post`,
   `change_password`).
2. Every successful `POST /login` MUST set `session['login_at']` to an epoch
   integer; the loader treats a missing or non-int marker as absent
   (reject-when-armed).
3. The check lives in the `user_loader` only — it must not be duplicated into
   `before_request` hooks; the loader is the single point covering remember
   cookie restores too.
4. Comparison is strict (`login_at < changed_at` rejects). Keep it that way:
   equality surviving is what would let a same-second acting session stay
   valid if the sign-out-everywhere rule is ever relaxed.
5. Never re-introduce boot-time `ALTER TABLE` (auto-migration). Manual
   `tools/*.sql` only; pinned by `tests/test_security_regressions.py::
   test_boot_does_not_alter_legacy_users_db`.
6. Automated tests: `tests/test_security_regressions.py` (invalidation set:
   reset kills sessions, change signs out everywhere, remember cookie
   cleared, boot does not alter legacy DB).

## Related code

| Path | Role |
|------|------|
| `poi_broker/__init__.py` | `user_loader` epoch check |
| `poi_broker/auth.py` | `login_at` stamp; watermark + sign-out on reset/change |
| `poi_broker/models.py` | `User.password_changed_at` column |
| `poi_broker/app.py` | `epoch_utc_date` Jinja filter (Security page display) |
| `tools/apply_password_changed_at.sql` | One-shot upgrade script |
| `tools/usersdb_schema.sql` | Fresh-install schema (column included) |
| `docs/password_reset/deployment.md` | Deploy/upgrade runbook |
| `tests/test_security_regressions.py` | Behavior guards |
