"""
Huey configuration for both development and production environments.

This module provides a flexible Huey setup that supports:
- Development: MemoryHuey with immediate execution (synchronous, no external dependencies)
- Production: SqliteHuey with persistent queue (async, requires huey_consumer worker)

Usage:
    from poi_broker.huey_config import get_huey
    huey = get_huey()

Environment Variables:
    HUEY_BACKEND: 'memory' (default) or 'sqlite' for production
    HUEY_SQLITE_PATH: Path to huey.db (default: instance/huey.db)
    HUEY_IMMEDIATE: 'true' or 'false' - run tasks synchronously (for development)
"""

import os
from pathlib import Path


def get_huey_backend():
    """
    Determine which Huey backend to use based on environment.
    
    Returns:
        str: 'memory' (development) or 'sqlite' (production)
    """
    backend = os.environ.get('HUEY_BACKEND', 'memory').lower()
    if backend not in ('memory', 'sqlite'):
        raise ValueError(f"Invalid HUEY_BACKEND: {backend}. Must be 'memory' or 'sqlite'")
    return backend


def get_huey_sqlite_path():
    """
    Get the SQLite database path for Huey.
    
    Returns:
        str: Absolute path to huey.db
    """
    # First check if explicitly provided
    explicit_path = os.environ.get('HUEY_SQLITE_PATH')
    if explicit_path:
        return explicit_path
    
    # Otherwise use instance folder (default for Flask)
    # For worker processes without Flask app, fall back to current directory
    try:
        from flask import current_app
        instance_path = current_app.instance_path
    except (RuntimeError, ImportError):
        # No Flask app context (e.g., in huey_consumer worker)
        instance_path = os.path.join(os.getcwd(), 'instance')
    
    Path(instance_path).mkdir(parents=True, exist_ok=True)
    return str(Path(instance_path) / 'huey.db')


def create_huey():
    """
    Create and configure Huey instance based on environment.
    
    Returns:
        Huey: Configured Huey instance (MemoryHuey or SqliteHuey)
    """
    backend_type = get_huey_backend()
    
    if backend_type == 'sqlite':
        from huey import SqliteHuey
        db_path = get_huey_sqlite_path()
        return SqliteHuey(
            'poi_broker',
            filename=db_path,
            journal_mode='wal',
            timeout=5,
            cache_mb=8,
            fsync=False,
            utc_time=True
        )
    else:
        # Default: MemoryHuey for development
        from huey import MemoryHuey
        immediate = os.environ.get('HUEY_IMMEDIATE', 'true').lower() == 'true'
        return MemoryHuey('poi_broker', immediate=immediate)
