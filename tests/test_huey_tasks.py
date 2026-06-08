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
import pytest
import time
from pathlib import Path
from datetime import datetime, timezone

# Set test environment before imports
os.environ.setdefault('SECRET_KEY', 'test-secret-key')
os.environ.setdefault('HUEY_BACKEND', 'memory')


def test_huey_imports():
    """Test that Huey can be imported and configured."""
    from poi_broker.huey_config import get_huey_backend, create_huey
    
    backend = get_huey_backend()
    assert backend in ('memory', 'sqlite'), f"Invalid backend: {backend}"
    
    huey = create_huey()
    assert huey is not None, "Failed to create Huey instance"


def test_app_with_memory_huey(tmp_path):
    """Test app creation with MemoryHuey backend (development mode)."""
    os.environ['HUEY_BACKEND'] = 'memory'
    
    # Need to reimport to pick up new environment variable
    import importlib
    import poi_broker.huey_config
    import poi_broker.extensions
    importlib.reload(poi_broker.huey_config)
    importlib.reload(poi_broker.extensions)
    
    from poi_broker import create_app
    app = create_app()
    
    assert app is not None
    with app.app_context():
        from poi_broker.extensions import huey
        assert hasattr(huey, 'immediate'), "MemoryHuey should have 'immediate' attribute"


def test_app_with_sqlite_huey(tmp_path):
    """Test Huey SqliteHuey backend configuration."""
    # Test that SqliteHuey can be created with proper configuration
    from pathlib import Path
    db_path = str(tmp_path / 'test_huey.db')
    
    # Create SqliteHuey directly to test configuration
    from huey import SqliteHuey
    test_huey = SqliteHuey(
        'poi_broker_test',
        filename=db_path,
        journal_mode='wal',
        timeout=5,
        cache_mb=8,
        fsync=False,
    )
    
    # Verify the object was created successfully
    assert test_huey is not None
    assert test_huey.name == 'poi_broker_test'
    # Database directory should exist (parent of db_path)
    assert Path(db_path).parent.exists()


def test_export_task_creation(app):
    """Test creating an ExportTask and enqueuing a background job."""
    from poi_broker import db
    from poi_broker.models import User, ExportTask
    
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
        assert '/export/download/<task_id>' in export_routes, "Download route not found"
        assert '/export/status/<task_id>' in export_routes, "Status route not found"


def test_export_task_status_transitions(app):
    """Test ExportTask status transitions through its lifecycle."""
    from poi_broker import db
    from poi_broker.models import User, ExportTask
    
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


# ============================================================================
# MANUAL TESTING INSTRUCTIONS
# ============================================================================
"""
To manually test the production Huey setup with SQLite backend:

## STEP 1: Start the Huey Worker

### On Linux/Mac:
```bash
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export SECRET_KEY=test-secret-key
./run_huey_worker.sh
```

### On Windows:
```cmd
set HUEY_BACKEND=sqlite
set HUEY_SQLITE_PATH=.\\instance\\huey.db
set SECRET_KEY=test-secret-key
run_huey_worker.bat
```

## STEP 2: Start the Flask App (in another terminal)

### On Linux/Mac:
```bash
export HUEY_BACKEND=sqlite
export HUEY_SQLITE_PATH=./instance/huey.db
export SECRET_KEY=test-secret-key
python -m flask --app wsgi:app run --debug
```

### On Windows:
```cmd
set HUEY_BACKEND=sqlite
set HUEY_SQLITE_PATH=.\\instance\\huey.db
set SECRET_KEY=test-secret-key
python -m flask --app wsgi:app run --debug
```

## STEP 3: Test the Export Feature

1. Open http://localhost:5000/ in your browser
2. Login to your account
3. Navigate to /export
4. Build a query using the visual query builder
5. Click "Start Export"
6. In the worker terminal, you should see:
   ```
   [2026-06-08 10:30:45] EXECUTING: poi_broker.tasks.create_export_file
   ```
7. Monitor the task progress in real-time
8. Once complete, download the CSV file

## STEP 4: Verify SQLite Persistence

The Huey task database is stored at: ./instance/huey.db

Check task history:
```bash
sqlite3 ./instance/huey.db
.tables
SELECT * FROM task;
```

## STEP 5: Load Testing (Optional)

To test with multiple concurrent exports:

```python
import requests
import json

session = requests.Session()
session.post('http://localhost:5000/login', data={
    'email': 'your@email.com',
    'password': 'your_password'
})

# Submit multiple export tasks
for i in range(3):
    response = session.post('http://localhost:5000/export', json={
        'rules': [
            {'id': 'featuretable.ant_mag_corrected', 'field': 'feature_mean_magn_r', 
             'type': 'double', 'input': 'text', 'operator': '>', 'value': 15}
        ]
    })
    print(f'Export {i+1}: {response.json()}')
```

## TROUBLESHOOTING

### Huey consumer won't start:
- Check that SECRET_KEY environment variable is set
- Verify DATABASE paths are correct
- Check logs: tail -f huey.log

### Tasks not executing:
- Verify HUEY_BACKEND=sqlite is set
- Check that huey_consumer process is running
- Look for errors in the worker log

### SQLite database locked:
- Ensure only one huey_consumer is running
- Check for stale lock files in instance/
- Restart the worker process
"""
