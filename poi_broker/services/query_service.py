"""Query building and execution service for visual query builder."""
import logging

from sqlalchemy.dialects import sqlite as _sqlite_dialect

from .. import db
from ..models import Classification, Ztf
from ..querybuilder_translator import Filter

logger = logging.getLogger(__name__)


# Full set of Ztf (object) columns. Used for filtering, watchlists and the bulk
# export, so it lives once here as the authoritative list.
ZTF_COLUMNS = [col.name for col in Ztf.__table__.columns] # type: ignore

# Classification columns surfaced in the bulk export. Data for these is pulled
# off a joined Classification row (which the query returns as the second element).
CLASSIFICATION_COLUMNS = ['classification_alert_id', 'classification_prob_class']


def get_export_columns():
    """
    Return the column names used when writing an export CSV row.

    Returns:
        list[str]: Ztf object columns followed by classification columns.
    """
    return ZTF_COLUMNS + CLASSIFICATION_COLUMNS


def build_export_row(ztf_row, classification_row):
    """
    Build a single CSV row dict (field name -> value) for a bulk export.

    The query built by build_query_from_rules selects (Ztf, Classification) via
    an outer join, so each result is a tuple (ztf_row, classification_or_none).

    Args:
        ztf_row: A Ztf instance.
        classification_row: A Classification instance or None.

    Returns:
        dict[str, Any]: Field names from get_export_columns() to values.
    """
    row_dict = {col: getattr(ztf_row, col) for col in ZTF_COLUMNS}
    if classification_row:
        row_dict['classification_alert_id'] = classification_row.alert_id
        row_dict['classification_prob_class'] = classification_row.prob_class
    else:
        row_dict['classification_alert_id'] = None
        row_dict['classification_prob_class'] = None
    return row_dict


def build_query_from_rules(rules_payload):
    """
    Build and execute a query from querybuilder rules.
    
    Args:
        rules_payload: Dict with 'rules' key containing filter rules
    
    Returns:
        tuple: (filtered_query, where_clause_str) or raises Exception if invalid
    
    Raises:
        ValueError: If rules payload is invalid
    """
    if not isinstance(rules_payload, dict) or 'rules' not in rules_payload:
        raise ValueError('Invalid querybuilder rules payload')
    
    if not isinstance(rules_payload.get('rules'), list) or len(rules_payload['rules']) == 0:
        raise ValueError('At least one filter rule is required')
    
    base_query = db.session.query(Ztf, Classification).outerjoin(
        Classification,
        Ztf.alert_id == Classification.alert_id
    )
    
    models_dict = {'featuretable': Ztf, 'classification': Classification}
    myfilter = Filter(models_dict, base_query)
    filtered_query = myfilter.querybuilder(rules_payload)
    
    where_clause = filtered_query.whereclause
    if where_clause is None:
        raise ValueError('No filter conditions could be built from the rules')
    
    compiled = where_clause.compile(
        dialect=_sqlite_dialect.dialect(),
        compile_kwargs={'literal_binds': True}
    )
    
    return filtered_query, str(compiled)


def get_preview_sql(rules_payload):
    """
    Get SQL preview string from rules payload.
    
    Args:
        rules_payload: Dict with 'rules' key containing filter rules
    
    Returns:
        str: SQL WHERE clause as string
    """
    _, where_clause_str = build_query_from_rules(rules_payload)
    return where_clause_str


def get_query_match_count(rules_payload):
    """
    Get number of matching records for a query.
    
    Args:
        rules_payload: Dict with 'rules' key containing filter rules
    
    Returns:
        int: Number of matching records
    """
    filtered_query, _ = build_query_from_rules(rules_payload)
    match_count = filtered_query.order_by(None).with_entities(db.func.count()).scalar()
    return int(match_count or 0)
