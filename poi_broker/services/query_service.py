"""Query building and execution service for visual query builder."""

import logging
from sqlalchemy.dialects import sqlite as _sqlite_dialect
from .. import db
from ..models import Ztf, Classification
from ..querybuilder_translator import Filter

logger = logging.getLogger(__name__)


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
    
    base_query = db.session.query(Ztf.alert_id, Classification).outerjoin(
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
