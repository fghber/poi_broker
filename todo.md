# Fixes
- [x] change order of Summary: Watchlist <-> Bookmarks
- [x] Rename "Saved Filters" to Bookmarked Filters

# Cleanup
- align jQuery/Bootstrap versions and integrity tags
- JSON responses differs on error types and escaping: Some endpoints use jsonify, others build Response/current_app.response_class with json.dumps or safe_serialize (query_features, query_crossmatches).
- request.query_string.decode('ascii'): Non-ASCII query strings can raise; utf-8 with errors policy is safer.
- implement `get_flashed_messages(with_categories=True)` globally for all categoires used: 'danger', 'success', 'info',  'warning'
- apply db.session.commit() pattern when
  - Violating a unique constraint (email column)
  - Violating a foreign key constraint
  - Violating a NOT NULL constraint
  ```
    try:
        db.session.commit()
    except IntegrityError:
        db.session.rollback()
        flash('Failed to generate password reset link. Please try again later.')
        return redirect(url_for('auth.forgot_password'))
    except Exception as e:
        db.session.rollback()
        logger.error('Database error during commit', exc_info=True)
        flash('Failed to generate password reset link. Please try again later.')
        return redirect(url_for('auth.forgot_password'))
  ```

# Considerations

- [x] consider extending rate-limiting to heavy read routes/endpoints (/query_features, /query_crossmatches, /download_alerts_csv, main /)
- [x] use a single grouped query to get_favorite_groups (favorites_service.py): For each group it runs Favorite.query.filter_by(group_id=g.id).count() — classic N+1. 
- Spinner handling logic is repeated across multiple files. Consider refactoring this into a reusable JavaScript module.
- Update main table data via AJAX/API calls instead of page loads/GET

# New Features
- Allow users to create custom observatory coordinates for the observing planning tool
[x] Save table filters (URL) as bookmark (My Search/Filter)
[x] Allow users selecting (up to 10) default features to plot
- Add Bulk Export based on Visual Query (Top 1000/Preview or All)
  - Create full CSV async, inform user when ready

# Future
- Indepentent Python Client API Export Package 
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