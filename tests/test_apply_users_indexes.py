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
    assert 'uix_export_task_one_active_per_user' in names
    conn.close()


def test_apply_users_indexes_enforces_one_active_export_per_user(tmp_path):
    import sqlite3

    db_path = tmp_path / 'users.db'
    conn = sqlite3.connect(db_path)
    _schema(conn)
    conn.executescript(APPLY_SQL.read_text(encoding='utf-8'))

    conn.execute(
        "INSERT INTO export_task (user_id, status, created_at, updated_at) "
        "VALUES (1, 'PENDING', '2020-01-01', '2020-01-01')"
    )
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO export_task (user_id, status, created_at, updated_at) "
            "VALUES (1, 'RUNNING', '2020-01-02', '2020-01-02')"
        )
    # Completed tasks do not count against the single-active-task invariant.
    conn.execute(
        "INSERT INTO export_task (user_id, status, created_at, updated_at) "
        "VALUES (1, 'SUCCESS', '2020-01-03', '2020-01-03')"
    )
    conn.commit()
    conn.close()


def test_apply_users_indexes_aborts_on_duplicate_active_export_tasks(tmp_path):
    """Mixed active statuses (e.g. 1 PENDING + 1 RUNNING) must fail the partial unique index."""
    import sqlite3

    db_path = tmp_path / 'users.db'
    conn = sqlite3.connect(db_path)
    _schema(conn)
    conn.execute(
        "INSERT INTO export_task (user_id, status, created_at, updated_at) "
        "VALUES (1, 'PENDING', '2020-01-01', '2020-01-01')"
    )
    conn.execute(
        "INSERT INTO export_task (user_id, status, created_at, updated_at) "
        "VALUES (1, 'RUNNING', '2020-01-02', '2020-01-02')"
    )
    conn.commit()
    conn.close()

    script = APPLY_SQL.read_text(encoding='utf-8')
    conn = sqlite3.connect(db_path)
    # The pre-check query finds the duplicate (user_id=1, count=2)
    precheck = conn.execute(
        """
        SELECT user_id, COUNT(*) AS duplicate_count
        FROM export_task
        WHERE status IN ('PENDING', 'RUNNING')
        GROUP BY user_id
        HAVING COUNT(*) > 1
        """
    ).fetchall()
    assert precheck == [(1, 2)]

    # The script as a whole must reject the table due to the unique constraint
    with pytest.raises(sqlite3.IntegrityError):
        conn.executescript(script)
    conn.close()


def test_active_export_statuses_lockstep_across_code_and_schemas():
    """ExportTask.ACTIVE_STATUSES must match the partial index predicates in SQL schemas."""
    from poi_broker.models import ExportTask

    expected_predicate = f"status IN ({', '.join(repr(s) for s in ExportTask.ACTIVE_STATUSES)})"
    apply_sql = APPLY_SQL.read_text(encoding='utf-8')
    schema_sql = (Path(__file__).resolve().parents[1] / 'tools' / 'usersdb_schema.sql').read_text(encoding='utf-8')

    assert expected_predicate in apply_sql, f"Missing {expected_predicate} in apply_users_indexes.sql"
    assert expected_predicate in schema_sql, f"Missing {expected_predicate} in usersdb_schema.sql"

