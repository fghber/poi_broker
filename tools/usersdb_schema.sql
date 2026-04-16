---users.db---
CREATE TABLE
    user (
        id INTEGER NOT NULL,
        email VARCHAR(100),
        password VARCHAR(100),
        name VARCHAR(1000),
        role VARCHAR(20),
        reset_token VARCHAR(128),
        reset_token_expires INTEGER,
        email_verified INTEGER DEFAULT 0,
        email_verification_token VARCHAR(128),
        PRIMARY KEY (id),
        UNIQUE (email)
    );
CREATE INDEX ix_user_email_verification_token ON user (email_verification_token);

CREATE TABLE
    favorite (
        id INTEGER PRIMARY KEY,
        user_id INTEGER NOT NULL,
        locus_id TEXT NOT NULL,
        created_at INTEGER,
        group_id INTEGER REFERENCES favorite_group (id) ON DELETE SET NULL,
        UNIQUE (user_id, locus_id),
        FOREIGN KEY (user_id) REFERENCES user (id) ON DELETE CASCADE
    );
CREATE INDEX ix_favorite_group_id ON favorite (group_id);
CREATE INDEX ix_favorite_user_id ON favorite (user_id);
CREATE INDEX ix_favorite_locus_id ON favorite (locus_id);

CREATE TABLE
    favorite_group (
        id INTEGER PRIMARY KEY,
        user_id INTEGER NOT NULL,
        name TEXT NOT NULL,
        created_at TEXT NOT NULL,
        UNIQUE (user_id, name),
        FOREIGN KEY (user_id) REFERENCES user (id) ON DELETE CASCADE
    );
CREATE INDEX ix_favorite_group_user_id ON favorite_group (user_id);

CREATE TABLE
    watchlist (
        id INTEGER PRIMARY KEY,
        user_id INTEGER NOT NULL,
        name TEXT NOT NULL,
        rules_json TEXT NOT NULL,
        sql_where TEXT NOT NULL,
        created_at INTEGER NOT NULL DEFAULT (strftime ('%s', 'now')),
        UNIQUE (user_id, name),
        FOREIGN KEY (user_id) REFERENCES user (id) ON DELETE CASCADE
    );
CREATE INDEX ix_watchlist_user_id ON watchlist (user_id);

CREATE TABLE 
    filter_bookmark (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        user_id INTEGER NOT NULL,
        name VARCHAR(200) NOT NULL,
        query_json TEXT NOT NULL,
        created_at TEXT NOT NULL,
        FOREIGN KEY (user_id) REFERENCES user(id) ON DELETE CASCADE
	);
CREATE INDEX ix_filter_bookmark_user_id ON filter_bookmark (user_id);

CREATE TABLE 
    user_settings (
        id INTEGER PRIMARY KEY,
        user_id INTEGER NOT NULL,
        default_feature_plot_columns TEXT,
        FOREIGN KEY (user_id) REFERENCES user(id) ON DELETE CASCADE,
        UNIQUE (user_id)
	);