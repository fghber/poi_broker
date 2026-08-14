"""apply_users_indexes.sql must detect duplicate bookmark names before UNIQUE."""

from pathlib import Path

import pytest

APPLY_SQL = Path(__file__).resolve().parents[1] / 'tools' / 'apply_users_indexes.sql'


def _schema(conn) -> None:
    conn.executescript(
        """
        CREATE TABLE filter_bookmark (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            name VARCHAR(200) NOT NULL,
            query_json TEXT NOT NULL,
            created_at TEXT NOT NULL
        );
        CREATE TABLE export_task (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            status VARCHAR(20) NOT NULL DEFAULT 'PENDING',
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        );
        """
    )


def test_apply_users_indexes_aborts_on_duplicate_bookmark_names(tmp_path):
    import sqlite3

    db_path = tmp_path / 'users.db'
    conn = sqlite3.connect(db_path)
    _schema(conn)
    conn.execute(
        "INSERT INTO filter_bookmark (user_id, name, query_json, created_at) "
        "VALUES (1, 'dup', '{}', '2020-01-01')"
    )
    conn.execute(
        "INSERT INTO filter_bookmark (user_id, name, query_json, created_at) "
        "VALUES (1, 'dup', '{}', '2020-01-02')"
    )
    conn.commit()
    conn.close()

    script = APPLY_SQL.read_text(encoding='utf-8')
    conn = sqlite3.connect(db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.executescript(script)

    names = {
        row[0]
        for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index'"
        )
    }
    assert 'idx_export_task_status_updated_at' in names
    assert 'idx_export_task_status_created_at' in names
    assert 'uix_filter_bookmark_user_name' not in names
    conn.close()


def test_apply_users_indexes_creates_unique_index_when_names_unique(tmp_path):
    import sqlite3

    db_path = tmp_path / 'users.db'
    conn = sqlite3.connect(db_path)
    _schema(conn)
    conn.execute(
        "INSERT INTO filter_bookmark (user_id, name, query_json, created_at) "
        "VALUES (1, 'a', '{}', '2020-01-01')"
    )
    conn.execute(
        "INSERT INTO filter_bookmark (user_id, name, query_json, created_at) "
        "VALUES (1, 'b', '{}', '2020-01-02')"
    )
    conn.commit()
    conn.close()

    script = APPLY_SQL.read_text(encoding='utf-8')
    conn = sqlite3.connect(db_path)
    conn.executescript(script)
    names = {
        row[0]
        for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index'"
        )
    }
    assert 'uix_filter_bookmark_user_name' in names
    conn.close()
