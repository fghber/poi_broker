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


def build_query_from_rules(rules_payload: dict) -> tuple[Query, dict]:
    """
    Build an ORM query from querybuilder rules.

    Filter values stay as bind parameters. Do not compile this query with
    ``literal_binds=True`` for execution; persist ``rules_payload`` and rebuild.

    Args:
        rules_payload: Dict with 'rules' key containing filter rules.

    Returns:
        tuple: (filtered_query, rules_payload)

    Raises:
        ValueError: If rules payload is invalid or builds no WHERE clause.
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

    if filtered_query.whereclause is None:
        raise ValueError('No filter conditions could be built from the rules')

    return filtered_query, rules_payload


def build_export_query_from_rules(rules_payload: dict) -> tuple[Query, dict]:
    """
    Build a column-only export query from querybuilder rules.

    Same filters as :func:`build_query_from_rules`, but the SELECT list is Core
    columns so iteration does not construct Ztf/Classification instances or
    populate the identity map. ``yield_per`` batches row construction; the
    sqlite3 driver still buffers DBAPI results, so this is not a true stream.

    Args:
        rules_payload: Dict with 'rules' key containing filter rules.

    Returns:
        tuple: (export_query, rules_payload)
    """
    filtered_query, rules_payload = build_query_from_rules(rules_payload)
    export_query = (
        filtered_query.enable_eagerloads(False)
        .with_entities(*EXPORT_SELECT_COLS)
        .execution_options(yield_per=EXPORT_YIELD_PER, stream_results=True)
    )
    return export_query, rules_payload


def get_preview_sql(rules_payload: dict) -> str:
    """
    Compile a display-only WHERE preview from rules.

    Inlines bind values for the UI / ``Watchlist.sql_where`` column. Never
    execute the returned string; re-run :func:`build_query_from_rules` instead.

    Args:
        rules_payload: Dict with 'rules' key containing filter rules.

    Returns:
        str: SQL WHERE clause as a human-readable preview.
    """
    filtered_query, _ = build_query_from_rules(rules_payload)
    compiled = filtered_query.whereclause.compile(
        dialect=_sqlite_dialect.dialect(),
        compile_kwargs={'literal_binds': True},
    )
    return str(compiled)


def get_query_match_count(rules_payload: dict) -> int:
    """
    Get number of matching records for a query.

    Args:
        rules_payload: Dict with 'rules' key containing filter rules.

    Returns:
        int: Number of matching records.
    """
    filtered_query, _ = build_query_from_rules(rules_payload)
    match_count = filtered_query.order_by(None).with_entities(db.func.count()).scalar()
    return int(match_count or 0)


def execute_watchlist(
    rules_payload: dict,
    start_mjd: float,
    end_mjd: float,
    limit: int,
) -> list[str]:
    """
    Return alert IDs matching rules in ``[start_mjd, end_mjd)``.

    Rebuilds the ORM filter from ``rules_payload``. Never executes stored SQL.

    Args:
        rules_payload: Querybuilder rules dict persisted as ``rules_json``.
        start_mjd: Inclusive MJD lower bound.
        end_mjd: Exclusive MJD upper bound.
        limit: Maximum number of alert IDs to return.

    Returns:
        list[str]: Matching ``alert_id`` values, newest first.
    """
    query, _ = build_query_from_rules(rules_payload)
    rows = (
        query.with_entities(Ztf.alert_id)
        .filter(Ztf.date_alert_mjd >= start_mjd, Ztf.date_alert_mjd < end_mjd)
        .order_by(Ztf.date_alert_mjd.desc())
        .limit(limit)
        .all()
    )
    return [alert_id for (alert_id,) in rows if alert_id]
