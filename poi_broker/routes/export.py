"""Export data API routes blueprint."""

import logging
import json
from flask import Blueprint, render_template, jsonify, request, send_file, flash, redirect, url_for
from flask_login import login_required, current_user
from pathlib import Path

from .. import db
from ..models import ExportTask
from ..tasks import create_export_file

logger = logging.getLogger(__name__)

export_bp = Blueprint('export', __name__, url_prefix='/export')


@export_bp.route('', methods=['GET'])
@login_required
def export_page():
    """
    Render the export page with query builder UI and current export status.
    
    GET /export -> HTML page with query builder and export status
    """
    # Get the most recent export task for the user
    recent_task = ExportTask.query.filter_by(user_id=current_user.id).order_by(
        ExportTask.created_at.desc()
    ).first()
    
    return render_template(
        'export.html',
        recent_task=recent_task
    )


@export_bp.route('', methods=['POST'])
@login_required
def export_submit():
    """
    Handle export form submission. Creates an ExportTask and enqueues the background job.
    
    POST /export with JSON: {rules: [...]} -> redirects to GET /export with message
    """
    data = request.get_json() if request.is_json else {}
    query_params = data.get('rules')
    
    if not query_params:
        return jsonify({'error': 'No query parameters provided'}), 400
    
    # Check if user already has an active (PENDING or RUNNING) export task
    active_task = ExportTask.query.filter_by(user_id=current_user.id).filter(
        ExportTask.status.in_(['PENDING', 'RUNNING'])
    ).first()
    
    if active_task:
        return jsonify({
            'error': 'You already have an active export task. Please wait for it to complete or download the file.',
            'task_id': active_task.id
        }), 409
    
    try:
        # Create new ExportTask record
        export_task = ExportTask(
            user_id=current_user.id,
            status='PENDING'
        )
        db.session.add(export_task)
        db.session.commit()
        logger.info(f'Created ExportTask {export_task.id} for user {current_user.id}')
        
        # Enqueue the background task
        create_export_file(
            query_params={'rules': query_params},
            user_id=current_user.id,
            task_id=export_task.id
        )
        logger.info(f'Enqueued background task for ExportTask {export_task.id}')
        
        return jsonify({
            'success': True,
            'message': 'Export started. This may take a few moments.',
            'task_id': export_task.id
        }), 202
    
    except Exception as e:
        logger.error(f'Failed to create export task: {str(e)}', exc_info=True)
        return jsonify({'error': 'Failed to start export task'}), 500


@export_bp.route('/download/<task_id>', methods=['GET'])
@login_required
def export_download(task_id: str):
    """
    Download an exported CSV file.
    
    GET /export/download/<task_id> -> CSV file download
    """
    export_task = db.session.get(ExportTask, task_id, synchronize_session=False)
    
    if not export_task:
        logger.warning(f'Export task {task_id} not found')
        flash('Export task not found', 'danger')
        return redirect(url_for('export.export_page'))
    
    if export_task.user_id != current_user.id:
        logger.warning(f'User {current_user.id} attempted to download task {task_id} from user {export_task.user_id}')
        flash('Unauthorized', 'danger')
        return redirect(url_for('export.export_page'))
    
    if export_task.status != 'SUCCESS':
        logger.warning(f'Export task {task_id} has status {export_task.status}, not ready for download')
        flash(f'Export task is not ready for download (status: {export_task.status})', 'warning')
        return redirect(url_for('export.export_page'))
    
    if not export_task.file_path or not Path(export_task.file_path).exists():
        logger.error(f'Export file not found at {export_task.file_path}')
        flash('Export file not found on disk', 'danger')
        return redirect(url_for('export.export_page'))
    
    try:
        return send_file(
            export_task.file_path,
            as_attachment=True,
            download_name=f'export_{task_id}.csv'
        )
    except Exception as e:
        logger.error(f'Failed to serve export file {export_task.file_path}: {str(e)}', exc_info=True)
        flash('Failed to download file', 'danger')
        return redirect(url_for('export.export_page'))


@export_bp.route('/status/<task_id>', methods=['GET'])
@login_required
def export_status(task_id: str):
    """
    Get the status of an export task.
    
    GET /api/export/status/<task_id> -> JSON with task status
    """
    export_task = db.session.get(ExportTask, task_id, synchronize_session=False)
    
    if not export_task:
        return jsonify({'error': 'Task not found'}), 404
    
    if export_task.user_id != current_user.id:
        return jsonify({'error': 'Unauthorized'}), 403
    
    return jsonify({
        'task_id': export_task.id,
        'status': export_task.status,
        'created_at': export_task.created_at.isoformat(),
        'updated_at': export_task.updated_at.isoformat(),
        'file_path': export_task.file_path,
        'error_message': export_task.error_message
    })
