"""Hybrid catalog pagination: keyset for date-only sort, OFFSET+1 otherwise."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Any

from sqlalchemy import func, select, tuple_
from sqlalchemy.orm import Query

from poi_broker import db
from poi_broker.models import Ztf
from poi_broker.services.catalog_query import CatalogQueryBuild

PAGE_SIZE = 100
CATALOG_TOTAL_TTL_SECONDS = 300.0

_catalog_total_cache: dict[str, float | int | None] = {
    'value': None,
    'expires': 0.0,
}


@dataclass(frozen=True)
class KeysetCursor:
    date_alert_mjd: float
    alert_id: str
    locus_id: str

    def as_before_params(self) -> dict[str, str]:
        return {
            'before_mjd': _format_mjd_cursor(self.date_alert_mjd),
            'before_alert_id': self.alert_id,
            'before_locus_id': self.locus_id,
        }

    def as_after_params(self) -> dict[str, str]:
        return {
            'after_mjd': _format_mjd_cursor(self.date_alert_mjd),
            'after_alert_id': self.alert_id,
            'after_locus_id': self.locus_id,
        }


def _format_mjd_cursor(value: float) -> str:
    return format(float(value), '.16g')


def clear_catalog_total_cache() -> None:
    _catalog_total_cache['value'] = None
    _catalog_total_cache['expires'] = 0.0


def cached_catalog_total(ttl_seconds: float = CATALOG_TOTAL_TTL_SECONDS) -> int:
    """Unfiltered featuretable COUNT(*), cached per process."""
    now = time.monotonic()
    cached_value = _catalog_total_cache['value']
    expires = _catalog_total_cache['expires']
    if cached_value is not None and now < float(expires):
        return int(cached_value)

    total = db.session.execute(select(func.count()).select_from(Ztf)).scalar_one()
    _catalog_total_cache['value'] = int(total)
    _catalog_total_cache['expires'] = now + ttl_seconds
    return int(total)


def count_matches(query: Query) -> int:
    """COUNT(*) on an already-filtered query; strips ORDER BY."""
    match_count = query.order_by(None).with_entities(func.count()).scalar()
    return int(match_count or 0)


def fetch_page(query: Query, page: int, page_size: int | None = None) -> tuple[list[Any], bool]:
    """OFFSET page without a COUNT. has_next comes from the extra row."""
    if page_size is None:
        page_size = PAGE_SIZE
    if page < 1:
        raise ValueError('page must be >= 1')
    items = query.limit(page_size + 1).offset((page - 1) * page_size).all()
    has_next = len(items) > page_size
    return items[:page_size], has_next


def _keyset_order(sort_desc: bool):
    if sort_desc:
        return (Ztf.date_alert_mjd.desc(), Ztf.alert_id.desc(), Ztf.locus_id.desc())
    return (Ztf.date_alert_mjd.asc(), Ztf.alert_id.asc(), Ztf.locus_id.asc())


def _keyset_tuple_filter(cursor: KeysetCursor, toward_smaller: bool):
    columns = tuple_(Ztf.date_alert_mjd, Ztf.alert_id, Ztf.locus_id)
    values = tuple_(cursor.date_alert_mjd, cursor.alert_id, cursor.locus_id)
    if toward_smaller:
        return columns < values
    return columns > values


def fetch_keyset(
    query: Query,
    *,
    cursor: KeysetCursor | None,
    direction: str,
    sort_desc: bool,
    page_size: int | None = None,
) -> tuple[list[Any], bool, bool]:
    """Return (items, has_next, has_prev) for date-only keyset pagination.

    direction: first | before | after | last
    """
    if page_size is None:
        page_size = PAGE_SIZE
    base = query.order_by(None)
    reverse_after_fetch = False
    has_next = False
    has_prev = False

    if direction == 'first':
        page_query = base.order_by(*_keyset_order(sort_desc))
        has_prev = False
    elif direction == 'last':
        page_query = base.order_by(*_keyset_order(not sort_desc))
        reverse_after_fetch = True
        has_next = False
    elif direction == 'before':
        if cursor is None:
            raise ValueError('before direction requires a cursor')
        toward_smaller = sort_desc
        page_query = base.filter(_keyset_tuple_filter(cursor, toward_smaller)).order_by(
            *_keyset_order(sort_desc)
        )
        has_prev = True
    elif direction == 'after':
        if cursor is None:
            raise ValueError('after direction requires a cursor')
        toward_smaller = not sort_desc
        page_query = base.filter(_keyset_tuple_filter(cursor, toward_smaller)).order_by(
            *_keyset_order(not sort_desc)
        )
        reverse_after_fetch = True
        has_next = True
    else:
        raise ValueError(f'unknown keyset direction: {direction}')

    items = page_query.limit(page_size + 1).all()
    has_extra = len(items) > page_size
    items = items[:page_size]
    if reverse_after_fetch:
        items = list(reversed(items))

    if direction == 'first':
        has_next = has_extra
    elif direction == 'last':
        has_prev = has_extra
    elif direction == 'before':
        has_next = has_extra
    elif direction == 'after':
        has_prev = has_extra

    return items, has_next, has_prev


def cursor_from_row(row: Any) -> KeysetCursor:
    return KeysetCursor(
        date_alert_mjd=float(row.date_alert_mjd),
        alert_id=row.alert_id,
        locus_id=row.locus_id,
    )


def parse_keyset_cursor(args, prefix: str) -> KeysetCursor | None:
    mjd_text = args.get(f'{prefix}_mjd')
    alert_id = args.get(f'{prefix}_alert_id')
    locus_id = args.get(f'{prefix}_locus_id')
    if not mjd_text or not alert_id or not locus_id:
        return None
    try:
        mjd_value = float(mjd_text)
    except (TypeError, ValueError):
        return None
    return KeysetCursor(date_alert_mjd=mjd_value, alert_id=alert_id, locus_id=locus_id)


def parse_keyset_request(args) -> tuple[str, KeysetCursor | None]:
    """Return (direction, cursor) for a date-only request."""
    if args.get('last') == '1':
        return 'last', None
    before = parse_keyset_cursor(args, 'before')
    if before is not None:
        return 'before', before
    after = parse_keyset_cursor(args, 'after')
    if after is not None:
        return 'after', after
    return 'first', None


def _request_page(args) -> int | None:
    if hasattr(args, 'get'):
        try:
            page = args.get('page', default=None, type=int)
            if page is not None:
                return page
        except TypeError:
            pass
        raw_page = args.get('page')
    else:
        raw_page = None
    if raw_page in (None, ''):
        return None
    try:
        return int(raw_page)
    except (TypeError, ValueError):
        return None


def should_use_keyset(is_date_only_sort: bool, args) -> bool:
    """Keyset unless an extra sort is active or this is a bare page jump."""
    if not is_date_only_sort:
        return False
    direction, _cursor = parse_keyset_request(args)
    if direction != 'first':
        return True
    page = _request_page(args)
    return page is None or page <= 1


def resolve_catalog_total(
    build: CatalogQueryBuild,
    *,
    items: list[Any],
    has_next: bool,
    page: int,
    use_keyset: bool,
    keyset_direction: str,
) -> tuple[int | None, bool]:
    """Return (total, count_deferred). Never COUNTs an unselective filtered set."""
    is_first = (not use_keyset and page == 1) or (use_keyset and keyset_direction == 'first')
    if not has_next and is_first:
        return len(items), False
    if not has_next and not use_keyset:
        return (page - 1) * PAGE_SIZE + len(items), False
    if not build.has_applied_filters:
        return cached_catalog_total(), False
    if build.has_cheap_equality:
        return count_matches(build.count_query), False
    return None, True


def last_page_from_total(total: int, page_size: int = PAGE_SIZE) -> int:
    if total <= 0:
        return 1
    return max(1, math.ceil(total / page_size))
