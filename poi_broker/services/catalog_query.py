"""Build the main-catalog filter/sort query from request args."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Mapping
from urllib.parse import urlencode

from sqlalchemy.orm import Query

from poi_broker import db
from poi_broker.models import Classification, Ztf
from poi_broker.services.search_service import SearchService

COMPLETE_ALERT_ID_RE = re.compile(r'^(?:ztf_candidate|lsst):\d{18,}$')
ALERT_ID_PREFIX_RE = re.compile(r'^(?:ztf|lsst)\D*$')
VALID_PASSBANDS = frozenset({'g', 'R', 'i'})
VALID_PROB_CLASSES = ('cvnova', 'e', 'lpv', 'puls', 'periodic_other', 'quas', 'sn', 'yso')
FILTER_KEYS = (
    'date',
    'date_alert_mjd',
    'alert_id',
    'ztf_object_id',
    'ant_passband',
    'locus_id',
    'locus_ra',
    'locus_dec',
    'magpsf',
    'prob_class',
)
PAGINATION_KEYS = frozenset(
    {
        'page',
        'last',
        'before_mjd',
        'before_alert_id',
        'before_locus_id',
        'after_mjd',
        'after_alert_id',
        'after_locus_id',
    }
)

search_service = SearchService()


@dataclass(frozen=True)
class CatalogQueryBuild:
    list_query: Query
    count_query: Query
    filter_warning: str
    has_applied_filters: bool
    has_cheap_equality: bool
    is_date_only_sort: bool
    date_sort_desc: bool


def _base_query(*, include_classification: bool) -> Query:
    query = db.session.query(Ztf)
    if include_classification:
        query = query.outerjoin(Classification, Ztf.alert_id == Classification.alert_id)
    return query


def _apply_catalog_filters(query: Query, args: Mapping) -> tuple[Query, str, bool, bool, bool]:
    """Apply main-page filters. Returns query, warning, has_filters, cheap_equality, prob_class."""
    filter_warning = ''
    has_applied_filters = False
    has_cheap_equality = False
    has_prob_class_filter = False

    if date_text := args.get('date'):
        parsed_date = search_service.parse_mjd_input(date_text)
        if parsed_date.values:
            query = search_service.apply_mjd_filter(query, Ztf.date_alert_mjd, parsed_date)
            has_applied_filters = True
        else:
            filter_warning += (
                'Date filter cannot be applied - Enter a valid ISO-date or 8-digit integer date '
                'of the form yyyymmdd, e.g. "20201207", or a range, e.g., "20201207 20201209".'
            )

    if mjd_text := args.get('date_alert_mjd'):
        parsed = search_service.parse_float_input(mjd_text)
        if parsed and parsed.values:
            query = search_service.apply_float_filter(query, Ztf.date_alert_mjd, parsed)
            has_applied_filters = True
        else:
            filter_warning += (
                'MJD filter cannot be applied - Enter a valid Modified Julian Date as a number, '
                'e.g. "59190.12", a range, e.g. "59190 59191", or a bound with > or <, e.g. ">59190".'
            )

    if args.get('alert_id'):
        alert_id = args.get('alert_id', '').strip()
        if COMPLETE_ALERT_ID_RE.match(alert_id):
            query = query.filter(Ztf.alert_id == alert_id)
            has_applied_filters = True
            has_cheap_equality = True
        elif ALERT_ID_PREFIX_RE.match(alert_id):
            query = query.filter(Ztf.alert_id.like(f'{alert_id}%'))
            has_applied_filters = True
        elif alert_id != '':
            filter_warning += (
                'Alert ID cannot be filter by partial IDs - Enter a full alert ID, e.g. '
                '"ztf_candidate:335155568501", or "lsst:170094456539709554", or just the catalog '
                'prefix, e.g. "ztf" or "lsst".'
            )

    if object_id_text := args.get('ztf_object_id'):
        query = query.filter(Ztf.ztf_object_id == object_id_text.strip())
        has_applied_filters = True
        has_cheap_equality = True

    if passband_text := args.get('ant_passband'):
        if passband_text in VALID_PASSBANDS:
            query = query.filter(Ztf.ant_passband == passband_text)
            has_applied_filters = True
        else:
            filter_warning += (
                'Passband filter cannot be applied - Enter a valid passband, e.g., "g", "R", or "i".'
            )

    if locus_id_text := args.get('locus_id'):
        query = query.filter(Ztf.locus_id == locus_id_text.strip())
        has_applied_filters = True
        has_cheap_equality = True

    if ra_text := args.get('locus_ra'):
        parsed = search_service.parse_float_input(ra_text, allowed_range=(0.0, 360.0))
        if parsed and parsed.values:
            query = search_service.apply_float_filter(query, Ztf.locus_ra, parsed, decimals=5)
            has_applied_filters = True
        else:
            filter_warning += (
                'Ra filter cannot be applied - Enter a valid number within the range 0° to 360°, '
                'e.g., "118.61421", or range, e.g., "80 90".'
            )

    if dec_text := args.get('locus_dec'):
        parsed = search_service.parse_float_input(dec_text, allowed_range=(-90.0, 90.0))
        if parsed and parsed.values:
            query = search_service.apply_float_filter(query, Ztf.locus_dec, parsed, decimals=5)
            has_applied_filters = True
        else:
            filter_warning += (
                'Dec filter cannot be applied - Enter a valid number within the range -90° to +90°, '
                'e.g., "-20.12345", or range, e.g., "14.5 29".'
            )

    if mag_text := args.get('magpsf'):
        parsed = search_service.parse_float_input(mag_text)
        if parsed and parsed.values:
            query = search_service.apply_float_filter(query, Ztf.ant_mag_corrected, parsed, decimals=3)
            has_applied_filters = True
        else:
            filter_warning += (
                'ant_mag_corrected filter cannot be applied - Enter a valid number, e.g., "18.84", '
                'or range, e.g., "18.8 19.4".'
            )

    if args.get('prob_class'):
        prob_class_value = args.get('prob_class', '').strip()
        if prob_class_value and prob_class_value in VALID_PROB_CLASSES:
            query = query.filter(Classification.prob_class == prob_class_value)
            has_applied_filters = True
            has_prob_class_filter = True
        elif prob_class_value:
            filter_warning += (
                'Classification filter cannot be applied - Enter a valid classification label, '
                f'e.g. "sn". Valid options are: {", ".join(VALID_PROB_CLASSES)}.'
            )

    return query, filter_warning, has_applied_filters, has_cheap_equality, has_prob_class_filter


def _apply_catalog_sorts(query: Query, args: Mapping) -> tuple[Query, bool, bool]:
    """Apply sort keys. Returns query, is_date_only_sort, date_sort_desc."""
    date_sort_desc = True
    if args.get('sort__date'):
        sort_date_order = args.get('sort__date')
        if sort_date_order == 'desc':
            query = query.order_by(Ztf.date_alert_mjd.desc())
            date_sort_desc = True
        if sort_date_order == 'asc':
            query = query.order_by(Ztf.date_alert_mjd.asc())
            date_sort_desc = False
    else:
        query = query.order_by(Ztf.date_alert_mjd.desc())

    is_date_only_sort = True
    sort_alert_order = args.get('sort__alert_id')
    if sort_alert_order:
        is_date_only_sort = False
        if sort_alert_order == 'desc':
            query = query.order_by(Ztf.alert_id.desc())
        if sort_alert_order == 'asc':
            query = query.order_by(Ztf.alert_id.asc())

    sort_object_order = args.get('sort__ztf_object_id')
    if sort_object_order:
        is_date_only_sort = False
        if sort_object_order == 'desc':
            query = query.order_by(Ztf.ztf_object_id.desc())
        if sort_object_order == 'asc':
            query = query.order_by(Ztf.ztf_object_id.asc())

    sort_ra_order = args.get('sort__locus_ra')
    if sort_ra_order:
        is_date_only_sort = False
        if sort_ra_order == 'desc':
            query = query.order_by(Ztf.locus_ra.desc())
        if sort_ra_order == 'asc':
            query = query.order_by(Ztf.locus_ra.asc())

    sort_dec_order = args.get('sort__locus_dec')
    if sort_dec_order:
        is_date_only_sort = False
        if sort_dec_order == 'desc':
            query = query.order_by(Ztf.locus_dec.desc())
        if sort_dec_order == 'asc':
            query = query.order_by(Ztf.locus_dec.asc())

    sort_mag_order = args.get('sort__ant_mag_corrected')
    if sort_mag_order:
        is_date_only_sort = False
        if sort_mag_order == 'desc':
            query = query.order_by(Ztf.ant_mag_corrected.desc())
        if sort_mag_order == 'asc':
            query = query.order_by(Ztf.ant_mag_corrected.asc())

    return query, is_date_only_sort, date_sort_desc


def project_catalog_columns(query: Query) -> Query:
    """Project list columns including classification label."""
    return query.with_entities(
        Ztf.date_alert_mjd,
        Ztf.alert_id,
        Ztf.ztf_object_id,
        Ztf.ant_passband,
        Ztf.locus_id,
        Ztf.locus_ra,
        Ztf.locus_dec,
        Ztf.ant_mag_corrected,
        Classification.prob_class.label('prob_class'),
    )


def build_catalog_query(args: Mapping, *, include_sort: bool = True) -> CatalogQueryBuild:
    """Build list and count queries from the same filter args."""
    list_query, warning, has_filters, cheap_equality, has_prob_class = _apply_catalog_filters(
        _base_query(include_classification=True),
        args,
    )
    count_query, _, _, _, _ = _apply_catalog_filters(
        _base_query(include_classification=has_prob_class),
        args,
    )

    date_sort_desc = True
    is_date_only_sort = True
    if include_sort:
        list_query, is_date_only_sort, date_sort_desc = _apply_catalog_sorts(list_query, args)

    return CatalogQueryBuild(
        list_query=list_query,
        count_query=count_query,
        filter_warning=warning,
        has_applied_filters=has_filters,
        has_cheap_equality=cheap_equality,
        is_date_only_sort=is_date_only_sort,
        date_sort_desc=date_sort_desc,
    )


def _encode_args(args: Mapping, allowed_keys: tuple[str, ...] | None, excluded: frozenset[str]) -> str:
    pairs: list[tuple[str, str]] = []
    getlist = getattr(args, 'getlist', None)
    keys = allowed_keys if allowed_keys is not None else tuple(args.keys())
    seen: set[str] = set()
    for key in keys:
        if key in excluded or key in seen:
            continue
        seen.add(key)
        values = getlist(key) if getlist else [args.get(key)]
        for value in values:
            if value is None or value == '':
                continue
            pairs.append((key, str(value)))
    return urlencode(pairs)


def catalog_query_string(args: Mapping) -> str:
    """Filter + sort query string with pagination/cursor keys stripped."""
    return _encode_args(args, None, PAGINATION_KEYS)


def catalog_filter_query_string(args: Mapping) -> str:
    """Filter-only query string for /api/catalog-count."""
    return _encode_args(args, FILTER_KEYS, PAGINATION_KEYS)


def catalog_href(query_string: str, extra: Mapping[str, str] | None = None) -> str:
    extra_qs = urlencode(list((extra or {}).items()))
    if query_string and extra_qs:
        return f'/?{query_string}&{extra_qs}'
    if extra_qs:
        return f'/?{extra_qs}'
    if query_string:
        return f'/?{query_string}'
    return '/'
