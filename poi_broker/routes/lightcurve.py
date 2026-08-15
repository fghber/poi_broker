"""Lightcurve visualization routes blueprint."""

import logging
from flask import Blueprint, Response, current_app, jsonify, request
from sqlalchemy.orm import Query
from ..services.plotting_service import bokeh_json_payload, create_bokeh_lightcurve_figure
from .. import db, limiter
from ..models import Ztf

logger = logging.getLogger(__name__)

lightcurve_bp = Blueprint('lightcurve', __name__)


def build_lightcurve_plot_query(locus_id: str) -> Query:
    """Column-projected lightcurve plot query for a locus."""
    return (
        db.session.query(Ztf.date_alert_mjd, Ztf.ant_mag_corrected, Ztf.ant_passband)
        .filter(Ztf.locus_id == locus_id)
    )


def build_lightcurve_csv_query(locus_id: str) -> Query:
    """Column-projected lightcurve CSV query for a locus."""
    return (
        db.session.query(Ztf.locus_id, Ztf.date_alert_mjd, Ztf.ant_mag_corrected)
        .filter(Ztf.locus_id == locus_id)
    )


@lightcurve_bp.route('/query_lightcurve_data', methods=['GET'])
@limiter.limit(lambda: current_app.config.get('READ_RATE_LIMIT_LAX', '30 per minute'))
def query_lightcurve_data():
    """
    Get lightcurve plot data for a locus ID.
    Returns JSON ``{div, script}`` Bokeh components.
    """
    locusId = request.args.get('locusId')
    if not locusId:
        return jsonify({'error': 'Missing locusId'}), 400

    try:
        data = build_lightcurve_plot_query(locusId).all()

        # Create Bokeh plot
        div, script = create_bokeh_lightcurve_figure(data)
        return jsonify(bokeh_json_payload(div, script))
    except Exception as e:
        logger.error(f'Error querying lightcurve data: {str(e)}', exc_info=True)
        return jsonify({'error': 'Error querying lightcurve data'}), 500


@lightcurve_bp.route('/locus_plot_csv', methods=['GET'])
@limiter.limit(lambda: current_app.config.get('READ_RATE_LIMIT_LAX', '30 per minute'))
def get_locus_plot():
    """
    Export lightcurve data as CSV for a locus ID.
    """
    locusId = request.args.get('locusId')
    if not locusId:
        return Response('Missing locusId', status=400)

    try:
        data = build_lightcurve_csv_query(locusId).all()

        csv = 'locus_id,date_alert_mjd,ant_mag_corrected\n'

        # Build CSV from query results
        for row in data:
            if row.date_alert_mjd is not None and row.ant_mag_corrected is not None:
                csv += f'{row.locus_id},{row.date_alert_mjd},{row.ant_mag_corrected}\n'

        return Response(
            csv,
            mimetype="text/csv",
            headers={"Content-disposition": "attachment; filename=myplot.csv"}
        )
    except Exception as e:
        logger.error(f'Error exporting lightcurve CSV: {str(e)}', exc_info=True)
        return Response('Error exporting lightcurve CSV', status=500)
