"""Feature query service for handling ZTF feature data."""

import logging
import json
from flask import current_app
from .. import db
from ..models import Ztf
from ..helpers import object_as_dict, safe_serialize
from ..constants.features import FEATURE_COLUMN_LIST, default_feature_plot_columns

logger = logging.getLogger(__name__)


def get_available_features():
    """
    Return the list of available feature column names.
    
    Returns:
        list: Feature column names (strings)
    """
    return FEATURE_COLUMN_LIST


def query_features_by_alert_id(alert_id):
    """
    Query all features for a given alert_id.
    
    Args:
        alert_id: The alert_id to query
    
    Returns:
        dict: Feature data as a dictionary, or None if not found
    """
    try:
        # Select only the feature columns to avoid returning unrelated fields like date_log
        selected_columns = [getattr(Ztf, col) for col in FEATURE_COLUMN_LIST]
        feature_query = db.session.query(*selected_columns)
        feature_query = feature_query.filter(Ztf.alert_id == alert_id)

        row = feature_query.first()
        if row is None:
            return None

        # Row is a SQLAlchemy Row object with keys matching FEATURE_COLUMN_LIST
        return {col: getattr(row, col) for col in FEATURE_COLUMN_LIST}
    except Exception as e:
        logger.error(f'Error querying features for alert_id {alert_id}: {str(e)}', exc_info=True)
        raise


def query_feature_plot_data(locus_id, selected_features=None):
    """
    Query feature plot data for a given locus_id.
    
    Args:
        locus_id: The locus_id to query
        selected_features: List of feature column names to load (optional, defaults via default_feature_plot_columns)
    
    Returns:
        list: Ztf rows with requested feature columns loaded
    """
    try:
        if not selected_features:
            feature_list = default_feature_plot_columns()
        else:
            # Limit to 10 features for plotting
            feature_list = selected_features[:10]
            # Validate against known features
            feature_list = [f for f in feature_list if f in FEATURE_COLUMN_LIST]
            if not feature_list:
                feature_list = default_feature_plot_columns()
        
        # Build column list: always include date and mag, plus selected features
        columns_to_load = [Ztf.date_alert_mjd, Ztf.ant_mag_corrected] + [getattr(Ztf, col) for col in feature_list]
        
        featureplot_query = db.session.query(Ztf)
        featureplot_query = featureplot_query.filter(Ztf.locus_id == locus_id)
        featureplot_query = featureplot_query.options(db.load_only(*columns_to_load))
        
        data = featureplot_query.all()
        return data, feature_list
    except Exception as e:
        logger.error(f'Error querying feature plot data for locus_id {locus_id}: {str(e)}', exc_info=True)
        raise
