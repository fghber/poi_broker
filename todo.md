# Fixes

- [x] change order of Summary: Watchlist <-> Bookmarks
- [x] Rename "Saved Filters" to Bookmarked Filters

# Cleanup

- [x] align jQuery/Bootstrap versions and integrity tags
- JSON responses differs on error types and escaping: Some endpoints use jsonify, others build Response/current_app.response_class with json.dumps or safe_serialize (query_features, query_crossmatches).
- [x] request.query_string.decode('ascii'): Non-ASCII query strings can raise; utf-8 with errors policy is safer.
- [x] implement `get_flashed_messages(with_categories=True)` globally for all categoires used: 'danger', 'success', 'info',  'warning'
- [x] apply db.session.commit() pattern when
  * Violating a unique constraint (email column)
  * Violating a foreign key constraint
  * Violating a NOT NULL constraint

# Considerations

- [x] consider extending rate-limiting to heavy read routes/endpoints (/query_features, /query_crossmatches, /download_alerts_csv, main /)
- [x] use a single grouped query to get_favorite_groups (favorites_service.py): For each group it runs Favorite.query.filter_by(group_id=g.id).count() — classic N+1. 
- [x] Consider refactoring the spinner into a reusable JavaScript module.
- Update main table data via AJAX/API calls instead of page loads/GET

# New Features

- [x] Save table filters (URL) as bookmark (My Search/Filter)
- [x] Allow users selecting (up to 10) default features to plot
- [x] Allow users to create custom observatory coordinates for the observing planning tool
- [x] Add Default Observatory Coordinates -> Last-used becomes the default for the next session
- [x] Document new features in the README.md
- [x] Add Bulk Export based on Visual Query (Top 1M/Preview or All)
  - [x] Create CSV fully async, inform user when ready

# Future

- Migrate to more capable DB (PostgreSQL)
- Change column type: ant_magband REAL -> TEXT
  ```
  ALTER TABLE featuretable ADD COLUMN ant_passband_str TEXT
  UPDATE featuretable SET ant_passband_str = CAST(ant_passband as TEXT)
  ALTER TABLE featuretable DROP COLUMN ant_passband;
  ALTER TABLE featuretable RENAME COLUMN ant_passband_str TO ant_passband;

  Recommended (safer, canonical SQLite way)
  CREATE TABLE new_featuretable (
        -- same schema, but:
        ant_passband TEXT,
        ...
    );
    INSERT INTO new_featuretable (...)
    SELECT
        CAST(ant_passband AS TEXT),
        ...
    FROM featuretable;
    DROP TABLE featuretable;
    ALTER TABLE new_featuretable RENAME TO featuretable;
    HOWEVER: SQLite happily stores TEXT in a REAL column (Type Affinity). This only matters when moving to PosgreSQL/MySQL
  ```