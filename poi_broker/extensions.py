"""
Extension instances for Flask app.

These are initialized without an app instance to avoid circular imports.
They are configured in the app factory (create_app).
"""

from .huey_config import create_huey


# Main Huey instance - created based on HUEY_BACKEND environment variable
# Development (default): MemoryHuey with immediate=True (synchronous)
# Production: SqliteHuey with persistent queue
huey = create_huey()
