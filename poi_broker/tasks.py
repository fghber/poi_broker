"""
Background tasks for data export using Huey.
"""

import csv
import logging
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

from huey import crontab

from . import create_app, db
from .extensions import huey
from .models import ExportTask
from .services.query_service import (
    build_export_row,
    build_query_from_rules,
    get_export_columns,
)

logger = logging.getLogger(__name__)


def _env_int(name: str, default: int) -> int:
    """Read an integer env var, falling back to ``default`` on missing/invalid."""
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


# Seconds a task may stay PENDING/RUNNING before it is considered abandoned.
STALE_TASK_MAX_AGE = _env_int('EXPORT_STALE_MAX_AGE_SECONDS', 30 * 60)

# Exports (CSV files + ExportTask rows) older than this are removed by the
# periodic housekeeping task.
EXPORT_RETENTION_DAYS = _env_int('EXPORT_RETENTION_DAYS', 10)


@huey.task()
def create_export_file(query_params: dict, user_id: int, task_id: int):
    """
    Background task to export data to CSV.
    
    Args:
        query_params: Dictionary with 'rules' key containing filter rules
        user_id: User ID for the export task
        task_id: Integer ID of the ExportTask record
    
    Updates:
        ExportTask status from PENDING to RUNNING then to SUCCESS/FAILED
        Sets file_path on SUCCESS or error_message on FAILED
    """
    app = create_app()
    
    with app.app_context():
        export_task = db.session.get(ExportTask, task_id)
        if not export_task:
            logger.exception(f'ExportTask {task_id} not found')
            return
        
        # Track any partially-written file so we can clean it up on failure.
        file_path: Path | None = None
        
        try:
            # Update status to RUNNING
            export_task.status = 'RUNNING'
            export_task.updated_at = datetime.now(timezone.utc)
            db.session.commit()
            logger.info(f'ExportTask {task_id} set to RUNNING')
            
            # Build and execute query
            filtered_query, where_clause = build_query_from_rules(query_params)
            logger.info(f'ExportTask {task_id} query: {where_clause}')
            
            # Create exports directory in instance path
            exports_dir = Path(app.instance_path) / 'exports'
            exports_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate filename with timestamp
            timestamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
            filename = f'export_{user_id}_{timestamp}.csv'
            file_path = exports_dir / filename
            
            # Write CSV file
            all_columns = get_export_columns()

            row_count = 0
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=all_columns)
                writer.writeheader()

                # Stream in batches instead of .all(): keeps memory bounded for
                # large exports (the Python sqlite3 driver still buffers the raw
                # DBAPI rows, but ORM objects are created/de-allocated per batch).
                for ztf_row, classification_row in filtered_query.yield_per(1000):
                    writer.writerow(build_export_row(ztf_row, classification_row))
                    row_count += 1
            
            logger.info(f'CSV file created at {file_path} with {row_count} rows')
            
            # Update task with SUCCESS status and file path
            export_task.status = 'SUCCESS'
            export_task.file_path = str(file_path)
            export_task.updated_at = datetime.now(timezone.utc)
            db.session.commit()
            logger.info(f'ExportTask {task_id} completed successfully')
            
        except Exception as e:
            logger.exception(f'ExportTask {task_id} failed with an error.')
            # Remove any partially-written CSV so failed exports don't leave
            # orphan files behind (cleanup_expired_exports only cleans SUCCESS).
            if file_path is not None:
                try:
                    file_path.unlink(missing_ok=True)
                except OSError:
                    logger.warning(f'Could not remove partial export file {file_path}')
            export_task.status = 'FAILED'
            export_task.error_message = str(e)
            export_task.updated_at = datetime.now(timezone.utc)
            try:
                db.session.commit()
            except Exception:
                db.session.rollback()
                logger.exception(f'Could not mark ExportTask {task_id} as FAILED')


def reset_stale_export_tasks(max_age_seconds: int = STALE_TASK_MAX_AGE) -> int:
    """
    Fail export tasks that are stuck in PENDING or RUNNING for too long.

    A crashed or restarted worker leaves tasks in PENDING/RUNNING forever, which
    permanently blocks the user's active-task slot (the API returns 409 while one
    exists). This resets such tasks to FAILED so the user can retry.

    Returns:
        int: number of tasks reset to FAILED
    """
    cutoff = datetime.now(timezone.utc) - timedelta(seconds=max_age_seconds)
    from .models import ExportTask
    stale = ExportTask.query.filter(
        ExportTask.status.in_(['PENDING', 'RUNNING']),
        ExportTask.updated_at < cutoff
    ).all()
    for task in stale:
        task.status = 'FAILED'
        task.error_message = (
            f"Export aborted: no progress for over {max_age_seconds // 60} minute(s). "
            "Please try again."
        )
        task.updated_at = datetime.now(timezone.utc)
    if stale:
        db.session.commit()
        logger.info(f'Marked {len(stale)} stale export task(s) as FAILED')
    return len(stale)


def cleanup_expired_exports(max_age_days: int = EXPORT_RETENTION_DAYS) -> int:
    """
    Remove exports older than ``max_age_days``: delete the CSV file, then the
    ExportTask DB row. Files are removed first so a crash mid-cleanup cannot
    leave a DB row pointing at a missing file (the download route already guards
    against missing files, but this keeps the DB consistent).

    Returns:
        int: number of exports removed.
    """
    cutoff = datetime.now(timezone.utc) - timedelta(days=max_age_days)
    from .models import ExportTask

    expired = ExportTask.query.filter(
        ExportTask.created_at < cutoff,
        ExportTask.status == 'SUCCESS',
    ).all()

    removed = 0
    for task in expired:
        # 1) Remove the CSV file first (best-effort).
        if task.file_path:
            try:
                Path(task.file_path).unlink(missing_ok=True)
            except OSError:
                logger.warning(f'Could not remove export file {task.file_path}')
        # 2) Remove the DB row.
        db.session.delete(task)
        removed += 1

    if removed:
        db.session.commit()
        logger.info(f'Removed {removed} expired export(s) older than {max_age_days} day(s)')
    return removed


@huey.periodic_task(crontab(minute='*/5'))
def cleanup_stale_export_tasks() -> int:
    """
    Periodic housekeeping task (every 5 minutes):
      - fail stale PENDING/RUNNING exports (reset_stale_export_tasks)
      - remove exports older than EXPORT_RETENTION_DAYS (cleanup_expired_exports)

    In memory/development mode this only executes when the consumer runs with
    periodic enabled; the startup reset in :func:`reset_stale_export_tasks`
    covers the app side.
    """
    app = create_app()
    with app.app_context():
        reset_stale_export_tasks()
        return cleanup_expired_exports()


@huey.periodic_task(crontab(hour='*/6'))
def prune_huey_results() -> int:
    """
    Periodic housekeeping for the Huey SQLite queue (every 6 hours).

    Huey stores task results in the ``taskresult`` table and never prunes them
    by default, so ``huey.db`` grows over time. This flushes stored results
    (and any stale schedule entries) to bound the queue database size.

    Returns:
        int: number of result rows flushed.
    """
    try:
        flushed = huey.storage.flush_results()
        # Also drop stale periodic schedule entries to keep the schedule table small.
        huey.storage.flush_schedule()
        logger.info(f'Pruned {flushed} huey result(s)')
        return int(flushed or 0)
    except Exception:
        logger.exception('Failed to prune huey results')
        return 0
