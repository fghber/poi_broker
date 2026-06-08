"""
Background tasks for data export using Huey.
"""

import csv
import logging
from pathlib import Path
from datetime import datetime, timezone

from .extensions import huey
from . import db, create_app
from .models import ExportTask, Ztf, Classification
from .services.query_service import build_query_from_rules

logger = logging.getLogger(__name__)


@huey.task()
def create_export_file(query_params: dict, user_id: int, task_id: str):
    """
    Background task to export data to CSV.
    
    Args:
        query_params: Dictionary with 'rules' key containing filter rules
        user_id: User ID for the export task
        task_id: UUID of the ExportTask record
    
    Updates:
        ExportTask status from PENDING to RUNNING then to SUCCESS/FAILED
        Sets file_path on SUCCESS or error_message on FAILED
    """
    app = create_app()
    
    with app.app_context():
        export_task = db.session.get(ExportTask, task_id, synchronize_session=False)
        if not export_task:
            logger.error(f'ExportTask {task_id} not found')
            return
        
        try:
            # Update status to RUNNING
            export_task.status = 'RUNNING'
            export_task.updated_at = datetime.now(timezone.utc)
            db.session.commit()
            logger.info(f'ExportTask {task_id} set to RUNNING')
            
            # Build and execute query
            filtered_query, where_clause = build_query_from_rules(query_params)
            results = filtered_query.all()
            logger.info(f'Query returned {len(results)} results for task {task_id}')
            
            # Create exports directory in instance path
            exports_dir = Path(app.instance_path) / 'exports'
            exports_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate filename with timestamp
            timestamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
            filename = f'export_{user_id}_{timestamp}.csv'
            file_path = exports_dir / filename
            
            # Write CSV file
            if results:
                # Extract column names from first result (Ztf object)
                ztf_obj = results[0][0]
                column_names = [col.name for col in Ztf.__table__.columns]
                classification_cols = ['classification_id', 'classification_type', 'classification_date']
                all_columns = column_names + classification_cols
                
                with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=all_columns)
                    writer.writeheader()
                    
                    for ztf_row, classification_row in results:
                        row_dict = {col.name: getattr(ztf_row, col.name) for col in Ztf.__table__.columns}
                        if classification_row:
                            row_dict['classification_id'] = classification_row.id
                            row_dict['classification_type'] = classification_row.classification_type
                            row_dict['classification_date'] = classification_row.classification_date
                        else:
                            row_dict['classification_id'] = None
                            row_dict['classification_type'] = None
                            row_dict['classification_date'] = None
                        writer.writerow(row_dict)
            else:
                # Write empty CSV with headers
                column_names = [col.name for col in Ztf.__table__.columns]
                classification_cols = ['classification_id', 'classification_type', 'classification_date']
                all_columns = column_names + classification_cols
                
                with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=all_columns)
                    writer.writeheader()
            
            logger.info(f'CSV file created at {file_path}')
            
            # Update task with SUCCESS status and file path
            export_task.status = 'SUCCESS'
            export_task.file_path = str(file_path)
            export_task.updated_at = datetime.now(timezone.utc)
            db.session.commit()
            logger.info(f'ExportTask {task_id} completed successfully')
            
        except Exception as e:
            logger.error(f'ExportTask {task_id} failed with error: {str(e)}', exc_info=True)
            export_task.status = 'FAILED'
            export_task.error_message = str(e)
            export_task.updated_at = datetime.now(timezone.utc)
            db.session.commit()
