"""
Integration test for Huey task queue in both development and production modes.

This test demonstrates:
1. Development mode (MemoryHuey, immediate=True, synchronous)
2. Production mode (SqliteHuey, async with background worker)

Usage:
    # Test development mode (default)
    python -m pytest tests/test_huey_tasks.py -v
    
    # Test production mode with SQLite backend
    export HUEY_BACKEND=sqlite
    export HUEY_SQLITE_PATH=./test_huey.db
    python -m pytest tests/test_huey_tasks.py -v
"""

import os
from datetime import datetime, timezone

import pytest

# Set test environment before imports
os.environ.setdefault('SECRET_KEY', 'test-secret-key')
os.environ.setdefault('HUEY_BACKEND', 'memory')


def test_huey_imports():
    """Test that Huey can be imported and configured."""
    from poi_broker.huey_config import create_huey, get_huey_backend
    
    backend = get_huey_backend()
    assert backend in ('memory', 'sqlite'), f"Invalid backend: {backend}"
    
    huey = create_huey()
    assert huey is not None, "Failed to create Huey instance"


def test_app_with_memory_huey(tmp_path):
    """Test app creation with MemoryHuey backend (development mode)."""
    os.environ['HUEY_BACKEND'] = 'memory'
    
    # Need to reimport to pick up new environment variable
    import importlib

    import poi_broker.extensions
    import poi_broker.huey_config
    importlib.reload(poi_broker.huey_config)
    importlib.reload(poi_broker.extensions)
    
    from poi_broker import create_app
    app = create_app()
    
    assert app is not None
    with app.app_context():
        from poi_broker.extensions import huey
        assert hasattr(huey, 'immediate'), "MemoryHuey should have 'immediate' attribute"


def test_app_with_sqlite_huey(tmp_path, monkeypatch):
    """Test SqliteHuey backend configuration through create_huey()."""
    # Regression: the sqlite branch previously passed utc_time=True (invalid) and
    # crashed at import. This test goes through create_huey() so the configured
    # path (huey_config) is exercised, not just a hand-built SqliteHuey.
    monkeypatch.setenv('HUEY_BACKEND', 'sqlite')
    monkeypatch.setenv('HUEY_SQLITE_PATH', str(tmp_path / 'test_huey.db'))

    from poi_broker.huey_config import create_huey
    test_huey = create_huey()

    assert test_huey is not None
    assert test_huey.name == 'poi_broker'
    # Database directory should exist (parent of db_path)
    assert (tmp_path / 'test_huey.db').parent.exists()
    # utc should be enabled at the Huey level
    assert test_huey.utc is True


def test_create_huey_requires_sqlite_path_without_app_context(tmp_path, monkeypatch):
    """Regression: worker (no Flask context) must not silently fall back to cwd."""
    monkeypatch.setenv('HUEY_BACKEND', 'sqlite')
    monkeypatch.delenv('HUEY_SQLITE_PATH', raising=False)

    from poi_broker.huey_config import get_huey_sqlite_path
    with pytest.raises(RuntimeError):
        get_huey_sqlite_path()


def test_worker_entrypoint_registers_tasks():
    """Regression: poi_broker.worker.huey must register create_export_file.

    huey_consumer only imports the module of the object passed to it. If the
    consumer targets poi_broker.extensions.huey, tasks are never imported and the
    worker runs but never executes anything. The worker entrypoint must import
    poi_broker.tasks so the registry is populated.

    Sibling tests reload `poi_broker.extensions`, which creates a *new* Huey
    instance whose registry is empty until tasks are re-registered against it.
    We re-import tasks + worker here so the assertion reflects the current shared
    instance rather than a stale one captured at collection time.
    """
    import importlib

    import poi_broker.extensions
    import poi_broker.tasks
    import poi_broker.worker

    # Reload extensions first so the registry is fresh (sibling tests may have
    # created a new Huey instance whose registry still holds the old task).
    importlib.reload(poi_broker.extensions)
    # Re-import so tasks register against the *current* shared Huey instance.
    importlib.reload(poi_broker.tasks)
    importlib.reload(poi_broker.worker)

    huey = poi_broker.extensions.huey
    # worker must be backed by the same instance as extensions
    assert poi_broker.worker.huey is huey
    task = huey._registry.string_to_task('poi_broker.tasks.create_export_file')
    assert task is not None, (
        "create_export_file not registered on the worker's Huey instance"
    )


def test_get_worker_app_is_process_singleton_under_concurrency(app, monkeypatch):
    """L2: concurrent Huey threads must share one Flask app, not N create_app()s."""
    import threading

    import poi_broker.tasks as tasks_mod

    create_calls = []
    entered = threading.Event()
    release = threading.Event()

    def _slow_create_app():
        create_calls.append(1)
        entered.set()
        assert release.wait(timeout=2), 'create_app was not released'
        return app

    monkeypatch.setattr(tasks_mod, 'create_app', _slow_create_app)

    n_threads = 8
    results: list = []
    errors: list = []

    def _worker():
        try:
            results.append(tasks_mod._get_worker_app())
        except Exception as exc:
            errors.append(exc)

    threads = [threading.Thread(target=_worker) for _ in range(n_threads)]
    for thread in threads:
        thread.start()

    assert entered.wait(timeout=2), 'create_app was never entered'
    release.set()
    for thread in threads:
        thread.join(timeout=2)
        assert not thread.is_alive()

    assert errors == []
    assert len(create_calls) == 1
    assert len(results) == n_threads
    assert all(got is app for got in results)


def test_reset_stale_export_tasks(app):
    """Stale PENDING/RUNNING exports are reset to FAILED."""
    from datetime import timedelta

    from poi_broker import db
    from poi_broker.models import ExportTask, User
    from poi_broker.tasks import reset_stale_export_tasks

    with app.app_context():
        ExportTask.query.delete()
        User.query.filter_by(email='stale@example.com').delete()
        User.query.filter_by(email='stale2@example.com').delete()
        User.query.filter_by(email='fresh@example.com').delete()
        db.session.commit()

        user1 = User(email='stale@example.com', password='hashed', name='Stale', email_verified=True)
        user2 = User(email='stale2@example.com', password='hashed', name='Stale2', email_verified=True)
        user3 = User(email='fresh@example.com', password='hashed', name='Fresh', email_verified=True)
        db.session.add_all([user1, user2, user3])
        db.session.commit()

        def _make(user, status, age_seconds):
            t = ExportTask(user_id=user.id, status=status)
            t.updated_at = datetime.now(timezone.utc) - timedelta(seconds=age_seconds)
            db.session.add(t)
            return t

        stale_pending = _make(user1, 'PENDING', 60 * 60)   # 1h old -> stale
        stale_running = _make(user2, 'RUNNING', 60 * 60)   # 1h old -> stale
        fresh_pending = _make(user3, 'PENDING', 60)        # 1m old -> not stale
        done = _make(user1, 'SUCCESS', 60 * 60)            # old but terminal -> untouched
        db.session.commit()

        n = reset_stale_export_tasks(max_age_seconds=30 * 60)
        db.session.refresh(stale_pending)
        db.session.refresh(stale_running)
        db.session.refresh(fresh_pending)
        db.session.refresh(done)

        assert n == 2
        assert stale_pending.status == 'FAILED'
        assert stale_running.status == 'FAILED'
        assert fresh_pending.status == 'PENDING'
        assert done.status == 'SUCCESS'

        ExportTask.query.delete()
        User.query.delete()
        db.session.commit()


def test_export_task_creation(app):
    """Test creating an ExportTask and enqueuing a background job."""
    from poi_broker import db
    from poi_broker.models import ExportTask, User
    
    with app.app_context():
        # Clean up any existing test data
        ExportTask.query.delete()
        User.query.filter_by(email='test@example.com').delete()
        db.session.commit()
        
        # Create a test user
        test_user = User(
            email='test@example.com',
            password='hashed_password',
            name='Test User',
            email_verified=True
        )
        db.session.add(test_user)
        db.session.commit()
        
        # Create an export task
        task = ExportTask(
            user_id=test_user.id,
            status='PENDING'
        )
        db.session.add(task)
        db.session.commit()
        
        # Verify task was created
        retrieved_task = ExportTask.query.filter_by(id=task.id).first()
        assert retrieved_task is not None
        assert retrieved_task.status == 'PENDING'
        assert retrieved_task.user_id == test_user.id
        
        # Clean up
        ExportTask.query.delete()
        User.query.delete()
        db.session.commit()


def test_huey_task_execution_memory_backend(app):
    """Test Huey task execution in development mode (immediate/synchronous)."""
    from poi_broker.extensions import huey
    
    with app.app_context():
        # Test that tasks can be created and registered with Huey
        
        # Create a simple test task
        @huey.task()
        def test_task(x):
            return x * 2
        
        # Execute the task - it should return immediately in memory mode
        result = test_task(21)
        
        # Verify the task was executed and returned a Result object
        # (even in immediate mode, Huey returns Result objects)
        assert result is not None
        # If immediate mode is enabled, we can get the actual value
        if hasattr(result, 'id'):
            # This is a Result object, which is expected for Huey tasks
            assert result.id is not None
        else:
            # Or it returned the direct value
            assert result == 42


def test_export_routes_exist():
    """Test that export routes are registered."""
    from poi_broker import create_app
    
    app = create_app()
    
    with app.app_context():
        rules = [str(rule) for rule in app.url_map.iter_rules()]
        export_routes = [r for r in rules if 'export' in r]
        
        assert '/export' in export_routes, "GET /export route not found"
        assert '/export/download/<int:task_id>' in export_routes, "Download route not found"
        assert '/export/status/<int:task_id>' in export_routes, "Status route not found"


def test_export_task_status_transitions(app):
    """Test ExportTask status transitions through its lifecycle."""
    from poi_broker import db
    from poi_broker.models import ExportTask, User
    
    with app.app_context():
        # Clean up
        ExportTask.query.delete()
        User.query.filter_by(email='status_test@example.com').delete()
        db.session.commit()
        
        # Create user
        user = User(
            email='status_test@example.com',
            password='hashed',
            name='Status Test',
            email_verified=True
        )
        db.session.add(user)
        db.session.commit()
        
        # Create task and verify status transitions
        task = ExportTask(user_id=user.id, status='PENDING')
        db.session.add(task)
        db.session.commit()
        
        # PENDING -> RUNNING
        task.status = 'RUNNING'
        db.session.commit()
        assert ExportTask.query.get(task.id).status == 'RUNNING'
        
        # RUNNING -> SUCCESS
        task.status = 'SUCCESS'
        task.file_path = '/tmp/export_123.csv'
        db.session.commit()
        assert ExportTask.query.get(task.id).status == 'SUCCESS'
        
        # Verify timestamps updated
        assert task.updated_at is not None
        assert task.created_at is not None
        
        # Clean up
        ExportTask.query.delete()
        User.query.delete()
        db.session.commit()


def test_create_export_file_cleans_partial_file_on_failure(app, monkeypatch):
    """Regression (#5): a mid-write failure removes the partial CSV file and
    marks the ExportTask FAILED instead of leaving an orphan file on disk.

    We make build_export_query_from_rules succeed but the row iteration raise,
    so the file has been opened (header written) before the exception — exactly
    the "disk full / row serialization" scenario. The task's error handler must
    unlink the partial file.
    """
    from pathlib import Path

    import poi_broker.tasks as tasks_mod
    from poi_broker import db
    from poi_broker.models import ExportTask, User

    with app.app_context():
        user = User(
            email='export_fail_cleanup@example.com',
            password='hashed',
            name='Fail Cleanup',
            email_verified=True,
        )
        db.session.add(user)
        db.session.commit()
        task = ExportTask(user_id=user.id, status='PENDING')
        db.session.add(task)
        db.session.commit()
        task_id = task.id
        user_id = user.id

    # build_export_query_from_rules succeeds, but iterating the result raises
    # inside the CSV-write loop (after the file is opened and the header written).
    def _boom_iter():
        raise RuntimeError('simulated write failure')
        yield  # pragma: no cover -- makes this a generator so iteration raises

    class _BoomQuery:
        def yield_per(self, _n):
            return iter(_boom_iter())

    def _fake_build(*_a, **_kw):
        return _BoomQuery(), 'fake where'

    monkeypatch.setattr(tasks_mod, 'build_export_query_from_rules', _fake_build)
    monkeypatch.setattr(tasks_mod, '_get_worker_app', lambda: app)

    # Snapshot the exports dir before the run. The fixture app's instance path
    # is the shared repo instance/ dir, which may already contain files from
    # earlier manual exports — so we assert that the task leaves no NEW file.
    exports_dir = Path(app.instance_path) / 'exports'
    exports_dir.mkdir(parents=True, exist_ok=True)
    before = {p.name for p in exports_dir.iterdir()}

    # Run the task synchronously via .func so we don't depend on Huey backend mode.
    tasks_mod.create_export_file.func(
        query_params={'rules': [{'x': 1}]}, user_id=user_id, task_id=task_id,
    )

    with app.app_context():
        task = db.session.get(ExportTask, task_id)
        assert task is not None
        assert task.status == 'FAILED'
        assert task.error_message == 'Export failed. Please try again.'
        assert 'simulated write failure' not in (task.error_message or '')
        assert task.file_path is None

    # No NEW partial CSV may be left behind by the failed task.
    after = {p.name for p in exports_dir.iterdir()}
    assert after == before, f'failed task left a partial file: {after - before}'


def test_create_export_file_does_not_resurrect_failed_task(app, monkeypatch):
    """F6: if the stale guard FAILED the row mid-write, the worker must not
    flip it back to SUCCESS and must discard the CSV.
    """
    from pathlib import Path

    import poi_broker.tasks as tasks_mod
    from poi_broker import db
    from poi_broker.models import ExportTask, User

    with app.app_context():
        user = User(
            email='export_resurrect@example.com',
            password='hashed',
            name='Resurrect',
            email_verified=True,
        )
        db.session.add(user)
        db.session.commit()
        task = ExportTask(user_id=user.id, status='PENDING')
        db.session.add(task)
        db.session.commit()
        task_id = task.id
        user_id = user.id

    # First yield_per iteration: simulate another process failing the task
    # (stale guard), then yield one row so the worker would otherwise succeed.
    failed_once = {'done': False}

    def _iter():
        if not failed_once['done']:
            failed_once['done'] = True
            with app.app_context():
                row = db.session.get(ExportTask, task_id)
                row.status = 'FAILED'
                row.error_message = 'Export aborted: no progress for over 30 minute(s).'
                row.updated_at = datetime.now(timezone.utc)
                db.session.commit()
        # One fake row so build_export_row is exercised if the worker continues.
        yield type('Row', (), {})()

    class _Query:
        def yield_per(self, _n):
            return iter(_iter())

    monkeypatch.setattr(
        tasks_mod,
        'build_export_query_from_rules',
        lambda *_a, **_kw: (_Query(), 'fake where'),
    )
    monkeypatch.setattr(
        tasks_mod,
        'build_export_row',
        lambda _row: {col: '' for col in tasks_mod.get_export_columns()},
    )
    monkeypatch.setattr(tasks_mod, '_get_worker_app', lambda: app)

    exports_dir = Path(app.instance_path) / 'exports'
    exports_dir.mkdir(parents=True, exist_ok=True)
    before = {p.name for p in exports_dir.iterdir()}

    tasks_mod.create_export_file.func(
        query_params={'rules': [{'x': 1}]}, user_id=user_id, task_id=task_id,
    )

    with app.app_context():
        task = db.session.get(ExportTask, task_id)
        assert task is not None
        assert task.status == 'FAILED'
        assert task.file_path is None
        assert 'aborted' in (task.error_message or '')

    after = {p.name for p in exports_dir.iterdir()}
    assert after == before, f'resurrected export left a file: {after - before}'


def test_export_heartbeat_keeps_running_task_fresh(app, monkeypatch):
    """Heartbeat refreshes updated_at so reset_stale_export_tasks leaves a live export alone."""
    from datetime import timedelta
    from pathlib import Path

    import poi_broker.tasks as tasks_mod
    from poi_broker import db
    from poi_broker.models import ExportTask, User
    from poi_broker.tasks import reset_stale_export_tasks

    with app.app_context():
        user = User(
            email='export_heartbeat@example.com',
            password='hashed',
            name='Heartbeat',
            email_verified=True,
        )
        db.session.add(user)
        db.session.commit()
        task = ExportTask(user_id=user.id, status='PENDING')
        db.session.add(task)
        db.session.commit()
        task_id = task.id
        user_id = user.id

    # Force a heartbeat on every row by setting the interval to 0.
    monkeypatch.setattr(tasks_mod, 'EXPORT_HEARTBEAT_SECONDS', 0)

    row_count = {'n': 0}

    def _iter():
        # Yield two rows so the heartbeat runs at least once inside the loop.
        yield type('Row', (), {})()
        row_count['n'] += 1
        yield type('Row', (), {})()
        row_count['n'] += 1

    class _Query:
        def yield_per(self, _n):
            return iter(_iter())

    monkeypatch.setattr(
        tasks_mod,
        'build_export_query_from_rules',
        lambda *_a, **_kw: (_Query(), 'fake where'),
    )
    monkeypatch.setattr(
        tasks_mod,
        'build_export_row',
        lambda _row: {col: '' for col in tasks_mod.get_export_columns()},
    )
    monkeypatch.setattr(tasks_mod, '_get_worker_app', lambda: app)

    heartbeats = []
    real_transition = tasks_mod._transition_export_task

    def _tracking_transition(task_id_arg, from_statuses, **values):
        ok = real_transition(task_id_arg, from_statuses, **values)
        if ok and from_statuses == ['RUNNING'] and values.get('status') == 'RUNNING':
            heartbeats.append(True)
        return ok

    monkeypatch.setattr(tasks_mod, '_transition_export_task', _tracking_transition)

    tasks_mod.create_export_file.func(
        query_params={'rules': [{'x': 1}]}, user_id=user_id, task_id=task_id,
    )

    assert heartbeats, 'expected at least one RUNNING heartbeat during write'
    assert row_count['n'] == 2

    with app.app_context():
        task = db.session.get(ExportTask, task_id)
        assert task is not None
        assert task.status == 'SUCCESS'
        file_path = task.file_path
        # A freshly heartbeated RUNNING row must not be stale-failed.
        task.status = 'RUNNING'
        task.updated_at = datetime.now(timezone.utc) - timedelta(seconds=10)
        db.session.commit()
        n = reset_stale_export_tasks(max_age_seconds=30 * 60)
        db.session.refresh(task)
        assert n == 0
        assert task.status == 'RUNNING'
        if file_path:
            Path(file_path).unlink(missing_ok=True)
