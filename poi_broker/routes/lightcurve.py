"""Lightcurve visualization routes blueprint."""

import logging
from flask import Blueprint, Response, request
from ..services.plotting_service import create_bokeh_lightcurve_figure
from .. import db
from ..models import Ztf

logger = logging.getLogger(__name__)

lightcurve_bp = Blueprint('lightcurve', __name__)


@lightcurve_bp.route('/query_lightcurve_data', methods=['GET'])
def query_lightcurve_data():
    """
    Get lightcurve plot data for a locus ID.
    Returns Bokeh HTML/script components.
    """
    locusId = request.args.get('locusId')
    if not locusId:
        return Response('Missing locusId', status=400)

    try:
        # Query lightcurve data for this locus_id
        lightcurve_query = db.session.query(Ztf)
        lightcurve_query = lightcurve_query.filter(Ztf.locus_id == locusId)
        lightcurve_query = lightcurve_query.options(
            db.load_only(Ztf.date_alert_mjd, Ztf.ant_mag_corrected, Ztf.ant_passband)
        )
        
        data = lightcurve_query.all()

        # Create Bokeh plot
        div, script = create_bokeh_lightcurve_figure(data)
        
        # Return the components to the HTML template
        return f'{div}{script}'
    except Exception as e:
        logger.error(f'Error querying lightcurve data: {str(e)}', exc_info=True)
        return Response(f'Error: {str(e)}', status=500)


@lightcurve_bp.route('/locus_plot_csv', methods=['GET'])
def get_locus_plot():
    """
    Export lightcurve data as CSV for a locus ID.
    """
    locusId = request.args.get('locusId')
    if not locusId:
        return Response('Missing locusId', status=400)

    try:
        lightcurve_query = db.session.query(Ztf)
        lightcurve_query = lightcurve_query.filter(Ztf.locus_id == locusId)
        lightcurve_query = lightcurve_query.options(
            db.load_only(Ztf.locus_id, Ztf.date_alert_mjd, Ztf.ant_mag_corrected)
        )
        data = lightcurve_query.all()

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
        return Response(f'Error: {str(e)}', status=500)
