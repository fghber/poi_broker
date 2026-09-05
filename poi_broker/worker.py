"""
Huey consumer entrypoint.

The ``huey_consumer`` CLI only imports the single object passed to it (via
``load_class``), e.g. ``poi_broker.extensions.huey``. Importing ``extensions``
alone does NOT register the application's task functions on the Huey instance,
because the ``@huey.task`` decorators live in other modules (``poi_broker.tasks``).
A worker that never imports those modules starts up cleanly but silently never
executes any queued task.

This module exists so the consumer targets a single importable path that:
  1. triggers creation of the shared ``huey`` instance (``extensions``),
  2. imports ``tasks`` to register every ``@huey.task`` on that instance.

Point ``huey_consumer`` at ``poi_broker.worker.huey`` instead of
``poi_broker.extensions.huey``.

Example:
    python -m huey.bin.huey_consumer poi_broker.worker.huey --workers=2 --worker-type=thread
"""

from .extensions import huey # noqa: I001 (intentionally imports huey from extensions before tasks to ensure proper initialization order)
from . import tasks  # noqa: F401 (registers create_export_file & friends on `huey`)

__all__ = ['huey']