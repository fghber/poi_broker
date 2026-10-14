# ADR: SQLite schema migration strategy (one-shot apply scripts, no framework)

> Completion/retirement rules and status notes below are historical, superseded
> by [Migration completion and current status](adr_migration_completion_status.md).
> Use that consolidated ledger for open work. Fresh-install-only schema files
> and the no-runtime-DDL decision remain in force.

- **Status:** Accepted
- **Date:** 2026-10
- **Scope:** all database schema changes (`users.db`, alerts DB); deployment docs
- **Disposition of:** auditSBA.md finding 4 (resolved by cleanup)

## Context

auditSBA.md finding 4 flagged that `tools/usersdb_schema.sql` is not
idempotent (bare `CREATE TABLE`s) while `docs/async_export/ubuntu_add_worker.md`
instructed applying it to an existing production `users.db`. The team's model
was already different: the full schema files are for **fresh installations**,
and schema changes for existing deployments travel as **one-shot
`apply_*.sql` scripts**. Production databases are fully migrated as of
2026-10 (maintainer statement; the parent repo `ninahernitschek/poi_broker_frontend`
`main` carries the same schema tooling), which makes four of the five then-existing
apply scripts spent artifacts.

## Decision

1. **The application never alters the database.** No auto-migration at boot —
   pinned by `tests/test_security_regressions.py::test_boot_does_not_alter_legacy_users_db`.
2. **Fresh installations** are built from `tools/usersdb_schema.sql` and
   `tools/alertsdb_schema.sql`. Those files are **never** run against an
   existing database (they are intentionally not idempotent — plain DDL catches
   mistakes on a fresh DB instead of silently skipping).
3. **Existing deployments** are migrated by one-shot `apply_*.sql` scripts.
   Each script's lifecycle:
   - **created** together with the fix (schema-file update + script + tests),
   - **applied** to production when the fix deploys,
   - **removed** from the repo once the fix is committed/merged into the
     parent `main` **and** applied in production.
4. Current state (2026-10): `apply_email_verification_expiry.sql`,
   `apply_password_changed_at.sql`, `apply_export_snapshot_mjd.sql`, and
   `apply_lookup_indexes.sql` were removed under this rule.
   **`tools/apply_users_indexes.sql` is retained**: it carries the
   one-active-export-per-user partial unique index
   (`uix_export_task_one_active_per_user`), a fix that was uncommitted at
   decision time — the parent's copy of the script does not contain it, so
   production has not received it yet. Remove the script after that fix is
   committed, merged to the parent, and applied in production.
5. **Caveat:** the parent repo also tracks these files, so an upstream merge
   may re-add deleted scripts — re-delete them.

## Alternatives considered

- **Migration framework** (Alembic etc.): rejected — against the project's
  SQLite-operations model (see `docs/async_export/adr_huey_sqlite.md`): two
  plain SQLite files, one web process, no migration state table to drift from
  reality.
- **Make the schema files idempotent** (`IF NOT EXISTS` everywhere) so they
  double as migrations: rejected — it blurs "fresh install DDL" (where plain
  DDL catches typos) with "migration" (where re-runnability matters), and
  hides partial-application states instead of failing loudly.
- **Delete all five apply scripts**: rejected — `apply_users_indexes.sql`
  carried a fix that was not yet committed anywhere; production had not
  received its index. Deleting it would have orphaned the fix.

## Consequences

- Removing an apply script is safe only when the fix it carries is committed,
  in the parent `main`, and applied in production — check all three, not just
  the file's presence in the parent.
- Historical script contents remain recoverable from git history.
- Deployment docs must describe *pending* migrations only; applied ones are
  removed from docs together with their scripts (historical mentions become
  past tense with a pointer here).

## Invariants (for implementers and coding agents)

1. Never run `usersdb_schema.sql` / `alertsdb_schema.sql` against an existing
   database.
2. Never make the app alter the database at boot or at runtime.
3. A DB change is not "done" until the schema files, the apply script, and the
   doc note exist together; the apply script is deleted only under rule 3
   above.
4. `tools/apply_users_indexes.sql` must not be deleted until
   `uix_export_task_one_active_per_user` is confirmed applied in production
   (its uncommitted fix must also be committed).

## Related

| Path | Role |
|------|------|
| `tools/usersdb_schema.sql` | users.db schema — fresh installs |
| `tools/alertsdb_schema.sql` | alerts DB schema — fresh installs |
| `tools/apply_users_indexes.sql` | retained: pending one-active-export-per-user migration |
| `tools/apply_fix_bad_magnitudes.sql` | retained: pending 83 corrupt-magnitude repair |
| `tools/apply_catalog_sort_indexes.sql` | retained: pending composite sort indexes |
| `docs/catalog_sort_migration.md` | pending runbook for the two scripts above |
| `tests/test_apply_users_indexes.py` | pins the retained script's invariants |
| `docs/password_reset/deployment.md` | historical runbook, kept as record |
| `tests/test_security_regressions.py::test_boot_does_not_alter_legacy_users_db` | pins "the app never alters the DB" |
