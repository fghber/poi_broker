"""Feature query routes blueprint."""

import logging
import json
from flask import Blueprint, Response, request, jsonify, current_app
from ..services.feature_service import query_features_by_alert_id, query_feature_plot_data, get_available_features
from ..services.plotting_service import create_bokeh_feature_plot
from ..helpers import safe_serialize

logger = logging.getLogger(__name__)

features_bp = Blueprint('features', __name__)


@features_bp.route('/query_features', methods=['GET'])
def query_features():
    """
    Get all features for a given alert_id.
    Returns JSON with all feature values.
    """
    alert_id = request.args.get('alert_id')
    if not alert_id:
        return Response('Missing alert_id', status=400)

    try:
        data = query_features_by_alert_id(alert_id)
        if data is None:
            return jsonify({'error': 'No feature record found for alert_id'}), 404

        response = current_app.response_class(
            response=json.dumps(data),
            status=200,
            mimetype='application/json'
        )
        return response
    except Exception as e:
        logger.error(f'Error querying features: {str(e)}', exc_info=True)
        return jsonify({'error': str(e)}), 500


@features_bp.route('/query_featureplot_data', methods=['GET'])
def query_featureplot_data():
    """
    Get feature plot data for a locus ID.
    Returns Bokeh HTML/script components with selected features plotted.
    """
    locusId = request.args.get('locusId')
    if not locusId:
        return Response('Missing locusId', status=400)

    selected_features = request.args.get('features')
    
    try:
        # Parse selected features if provided
        feature_list = None
        if selected_features:
            feature_list = [f.strip() for f in selected_features.split(',') if f.strip()]
            if not feature_list:
                feature_list = None
        
        # Query feature plot data
        data, used_features = query_feature_plot_data(locusId, feature_list)

        # Create Bokeh plot
        div, script = create_bokeh_feature_plot(data, used_features)
        
        # Return the components to the HTML template
        return f'{div}{script}'
    except Exception as e:
        logger.error(f'Error querying feature plot data: {str(e)}', exc_info=True)
        return Response(f'Error: {str(e)}', status=500)
