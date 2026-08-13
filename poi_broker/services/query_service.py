"""Query building and execution service for visual query builder."""
from __future__ import annotations

import logging
from typing import Any

from sqlalchemy.dialects import sqlite as _sqlite_dialect
from sqlalchemy.engine import Row
from sqlalchemy.orm import Query

from .. import db
from ..models import Classification, Ztf
from ..querybuilder_translator import Filter

logger = logging.getLogger(__name__)


# Full set of Ztf (object) columns. Used for filtering, watchlists and the bulk
# export, so it lives once here as the authoritative list.
ZTF_COLUMNS = [col.name for col in Ztf.__table__.columns] # type: ignore

# Classification columns surfaced in the bulk export CSV header.
CLASSIFICATION_COLUMNS = ['classification_alert_id', 'classification_prob_class']

# Core column projection for bulk CSV export. Selecting columns (not entities)
# skips ORM identity-map hydration. Classification is limited to the two CSV
# fields instead of the full mapped table. Ztf still exports every mapped
# column so the CSV schema is unchanged.
EXPORT_SELECT_COLS = [getattr(Ztf, name) for name in ZTF_COLUMNS] + [
    Classification.alert_id.label('classification_alert_id'),
    Classification.prob_class.label('classification_prob_class'),
]

EXPORT_YIELD_PER = 1000


def get_export_columns() -> list[str]:
    """
    Return the column names used when writing an export CSV row.

    Returns:
        list[str]: Ztf object columns followed by classification columns.
    """
    return ZTF_COLUMNS + CLASSIFICATION_COLUMNS


def build_export_row(row: Row) -> dict[str, Any]:
    """
    Build a single CSV row dict from a Core column Row.

    Args:
        row: A SQLAlchemy Row from :func:`build_export_query_from_rules`.

    Returns:
        dict[str, Any]: Field names from get_export_columns() to values.
    """
    mapping = row._mapping
    return {name: mapping[name] for name in get_export_columns()}


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


def build_export_query_from_rules(rules_payload: dict) -> tuple[Query, str]:
    """
    Build a column-only export query from querybuilder rules.

    Same filters as :func:`build_query_from_rules`, but the SELECT list is Core
    columns so iteration does not construct Ztf/Classification instances or
    populate the identity map. ``yield_per`` batches row construction; the
    sqlite3 driver still buffers DBAPI results, so this is not a true stream.

    Args:
        rules_payload: Dict with 'rules' key containing filter rules.

    Returns:
        tuple: (export_query, where_clause_str)
    """
    filtered_query, where_clause = build_query_from_rules(rules_payload)
    export_query = (
        filtered_query.enable_eagerloads(False)
        .with_entities(*EXPORT_SELECT_COLS)
        .execution_options(yield_per=EXPORT_YIELD_PER, stream_results=True)
    )
    return export_query, where_clause


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
