from flask import (render_template, abort, jsonify, request, Response,
                   redirect, url_for, make_response, Blueprint, current_app)
from flask_login import login_required, current_user
import logging
from functools import lru_cache
import csv
import io
from datetime import datetime, timezone
from astropy.time import Time
from astropy.coordinates import EarthLocation
import re
import json

from .helpers import extract_numbers, extract_dates, extract_float_filter, extract_mjd_filter, safe_serialize, result_to_dict, object_as_dict
from . import db, limiter
from .models import Ztf, Crossmatches, User, Favorite, FavoriteGroup, Watchlist, Classification, UserObservatory
from .routes import favorites_bp, filter_bookmarks_bp, visual_query_bp, lightcurve_bp, features_bp, user_observatories_bp
from .constants.features import FEATURE_COLUMNS, default_feature_plot_columns
from .user_settings import user_settings_bp, get_saved_feature_plot_columns, get_saved_last_selected_observatory, UserSettings  # noqa: F401 - imported to register model
from importlib.metadata import version
bokeh_version = version("bokeh")

logger = logging.getLogger(__name__)
main_blueprint = Blueprint('main', __name__)

@lru_cache(maxsize=8192)
def _format_mjd_cached(mjd_value: float) -> str:
    jdate = mjd_value + 2400000.5
    dt = Time(jdate, format='jd', scale='utc').to_datetime(timezone=timezone.utc) # Convert to timezone-aware datetime in UTC
    return  dt.strftime('%Y-%m-%d %H:%M:%S')


@lru_cache(maxsize=1)
def _get_builtin_observatory_options() -> list[dict[str, str]]:
    builtin_names = list(EarthLocation.get_site_names())
    return [{'value': f'builtin:{name}', 'label': name} for name in builtin_names]


def _resolve_selected_observatory(
    selected_meta: dict | None,
    builtin_values: set[str],
    custom_values: set[str]
) -> str | None:
    """Validate and resolve a saved observatory selection using O(1) lookup.
    
    Args:
        selected_meta: Dict from get_saved_last_selected_observatory() with 'source' key
                       ('builtin' or 'custom') and selection identifier (name or id).
        builtin_values: Set of valid builtin option values for efficient lookup.
        custom_values: Set of valid custom option values for efficient lookup.
    
    Returns:
        Valid option value (e.g., 'builtin:Palomar', 'custom:42') or None if selection
        is invalid, deleted, or doesn't match saved metadata.
    """
    if not isinstance(selected_meta, dict):
        return None
    
    source = selected_meta.get('source')
    
    if source == 'custom':
        selected_id = selected_meta.get('id')
        if isinstance(selected_id, int):
            candidate = f'custom:{selected_id}'
            if candidate in custom_values:
                return candidate
    elif source == 'builtin':
        selected_name = selected_meta.get('name')
        if isinstance(selected_name, str):
            candidate = f'builtin:{selected_name}'
            if candidate in builtin_values:
                return candidate
    
    return None


def _build_observatory_context() -> dict[str, list | str | None]:
    """Build observatory context for the main page.
    
    Constructs builtin and custom observatory options, restores the user's last
    selected observatory (if authenticated), and provides a fallback to the first
    available option if the saved selection is invalid or deleted.
    
    Returns:
        dict with keys:
            - 'builtin_options': list[dict] of {'value', 'label'} for astropy sites
            - 'custom_options': list[dict] of {'value', 'label'} for user observatories
                                (empty if not authenticated)
            - 'selected_value': str (e.g., 'builtin:Palomar') or None if no options available
    """
    builtin_options = _get_builtin_observatory_options()
    builtin_values = {opt['value'] for opt in builtin_options}

    custom_options = []
    custom_values = set()
    
    if current_user.is_authenticated:
        custom_rows = (
            UserObservatory.query.filter_by(user_id=current_user.id)
            .order_by(UserObservatory.name.asc())
            .all()
        )
        custom_options = [{'value': f'custom:{row.id}', 'label': row.name} for row in custom_rows]
        custom_values = {opt['value'] for opt in custom_options}

    # Attempt to restore saved selection; gracefully falls back if invalid/deleted
    selected_value = _resolve_selected_observatory(
        get_saved_last_selected_observatory(current_user.id) if current_user.is_authenticated else None,
        builtin_values,
        custom_values
    )

    # Fallback to first available option if no valid saved selection
    if selected_value is None:
        if builtin_options:
            selected_value = builtin_options[0]['value']
        elif custom_options:
            selected_value = custom_options[0]['value']

    return {
        'custom_options': custom_options,
        'builtin_options': builtin_options,
        'selected_value': selected_value,
    }


@main_blueprint.route('/', methods=['GET'])
@limiter.limit(lambda: current_app.config.get('READ_RATE_LIMIT_LAX', '30 per minute'))
def start():
    logger.info('Request with request_args: %s', json.dumps(request.args))
    
    page = request.args.get('page', 1, type=int)
    filter_warning_message = ''

    query = db.session.query(Ztf).outerjoin(Classification, Ztf.alert_id == Classification.alert_id)

    if request.args.get('date'):
        date_input = extract_dates(request.args.get('date'))
        if date_input:
            query = extract_mjd_filter(date_input, Ztf.date_alert_mjd, query)
        else:
            filter_warning_message += 'Date filter cannot be applied - Enter a valid ISO-date or 8-digit integer date of the form yyyymmdd, e.g. "20201207", or a range, e.g., "20201207 20201209".'
    
    if request.args.get('date_alert_mjd'):
        date_input = extract_numbers(request.args.get('date_alert_mjd'))
        if date_input != None:
            query = extract_float_filter(date_input, Ztf.date_alert_mjd, query)
        else:
            filter_warning_message += 'MJD filter cannot be applied - Enter a valid Modified Julian Date as a number, e.g. "59190.12", a range, e.g. "59190 59191", or a bound with > or <, e.g. ">59190".'

    if request.args.get('alert_id'):
        alertId = request.args.get('alert_id', '').strip()
        if re.match(r'^(?:ztf_candidate|lsst):\d{18,}$', alertId): # if it contains 18+ digits, we got a complete alert_id and can query it directly.
            query = query.filter(Ztf.alert_id == alertId)
        elif re.match(r'^(?:ztf|lsst)\D*$', alertId): # otherwise, we allow for partial matching of alert_id to find alerts with a specific prefix, e.g. ztf / lsst
            search = "{}%".format(alertId)
            query = query.filter(Ztf.alert_id.like(search))
        elif alertId != '':
            filter_warning_message += 'Alert ID cannot be filter by partial IDs - Enter a full alert ID, e.g. "ztf_candidate:335155568501", or "lsst:170094456539709554", or just the catalog prefix, e.g. "ztf" or "lsst".'

    if request.args.get('ztf_object_id'):
        query = query.filter(Ztf.ztf_object_id == request.args.get('ztf_object_id'))

    #if request.args.get('filter'):
        #query = query.filter(Ztf.filter == int(request.args.get('filter'))) # 1:g, 2:r, 3:i
    if request.args.get('ant_passband'):
        ant_passband = request.args.get('ant_passband')
        query = query.filter(Ztf.ant_passband == ant_passband) # g, R, i

    if request.args.get('locus_id'):
        query = query.filter(Ztf.locus_id == request.args.get('locus_id'))

    if request.args.get('locus_ra'):
        ra_input = extract_numbers(request.args.get('locus_ra'))
        if ra_input != None:
            query = extract_float_filter(ra_input, Ztf.locus_ra, query, decimals=5)
        else:
            filter_warning_message += 'Ra filter cannot be applied - Enter a valid number, e.g., "118.61421", or range, e.g., "80 90".'

    if request.args.get('locus_dec'):
        dec_input = extract_numbers(request.args.get('locus_dec'))
        if dec_input != None:
            query = extract_float_filter(dec_input, Ztf.locus_dec, query, decimals=5)
        else:
            filter_warning_message += 'Dec filter cannot be applied - Enter a valid number, e.g., "-20.02131", or range, e.g., "18.8 19.4".'

    if request.args.get('magpsf'):
        magpsf_input = extract_numbers(request.args.get('magpsf'))
        if magpsf_input != None:
            query = extract_float_filter(magpsf_input, Ztf.ant_mag_corrected, query, decimals=3)
        else:
            filter_warning_message += 'ant_mag_corrected filter cannot be applied - Enter a valid number, e.g., "18.84", or range, e.g., "18.8 19.4".'

    if request.args.get('prob_class'):
        prob_class_value = request.args.get('prob_class', '').strip()
        valid_prob_classes = ['cvnova', 'e', 'lpv', 'puls', 'periodic_other', 'quas', 'sn', 'yso']
        if prob_class_value and prob_class_value in valid_prob_classes:
            query = query.filter(Classification.prob_class == prob_class_value)
        elif prob_class_value and prob_class_value.strip():  # Warn if non-empty but not in valid options
            filter_warning_message += f'Classification filter cannot be applied - Enter a valid classification label, e.g. "sn". Valid options are: {", ".join(valid_prob_classes)}.'

    #Sort order by date (still sorts by mjd column)
    if request.args.get('sort__date'):
        sort__date_order = request.args.get('sort__date')
        if sort__date_order == 'desc':
            query = query.order_by(Ztf.date_alert_mjd.desc())
        if sort__date_order == 'asc':
            query = query.order_by(Ztf.date_alert_mjd.asc())
    else:
        query = query.order_by(Ztf.date_alert_mjd.desc()) #default sort order

    # Sort order by alert_id
    sort__alert_order = request.args.get('sort__alert_id')
    if sort__alert_order:
        if sort__alert_order == 'desc':
            query = query.order_by(Ztf.alert_id.desc())
        if sort__alert_order == 'asc':
            query = query.order_by(Ztf.alert_id.asc())
    
    # Sort order by ztf_object_id
    sort__object_order = request.args.get('sort__ztf_object_id')
    if sort__object_order:
        if sort__object_order == 'desc':
            query = query.order_by(Ztf.ztf_object_id.desc())
        if sort__object_order == 'asc':
            query = query.order_by(Ztf.ztf_object_id.asc())

    # Sort order by locus_ra
    sort__ra_order = request.args.get('sort__locus_ra')
    if sort__ra_order:
        if sort__ra_order == 'desc':
            query = query.order_by(Ztf.locus_ra.desc())
        if sort__ra_order == 'asc':
            query = query.order_by(Ztf.locus_ra.asc())

    # Sort order by locus_dec
    sort__dec_order = request.args.get('sort__locus_dec')
    if sort__dec_order:
        if sort__dec_order == 'desc':
            query = query.order_by(Ztf.locus_dec.desc())
        if sort__dec_order == 'asc':
            query = query.order_by(Ztf.locus_dec.asc())

    # Sort order by ant_mag_corrected
    sort__mag_order = request.args.get('sort__ant_mag_corrected')
    if sort__mag_order:
        if sort__mag_order == 'desc':
            query = query.order_by(Ztf.ant_mag_corrected.desc())
        if sort__mag_order == 'asc':
            query = query.order_by(Ztf.ant_mag_corrected.asc())

    # Project the classification label directly to avoid per-row lazy loads.
    query = query.with_entities(
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

    #latest = db.session.query(Ztf).order_by(Ztf.date_alert_mjd.desc()).first() # ? IDEA: show latest update date
    #print(f'Latest alert in DB has date_log={latest.date_log} (UTC {Time(latest.date_alert_mjd, format="mjd").iso})') #DEBUG: print latest alert date in MJD and UTC for debugging
    #print(query.statement.compile(compile_kwargs={"literal_binds": True})) #DEBUG: print the resulting SQL query
    paginator = query.paginate(page=page, per_page=100, error_out=True)

    observatory_context = _build_observatory_context()

    return render_template(
        "main.html",
        total_queries=paginator.total,
        table=paginator.items,
        page=paginator.page,
        has_next=paginator.has_next,
        last_page=paginator.pages,
        # ? TODO Pagination query-string re.sub may leave a trailing & in edge cases. TEST
        query_string=re.sub('[&?]?page=\\d+|&$', '', request.query_string.decode('utf-8')), # ? b'' binary string 
        filter_warning = filter_warning_message,
        custom_observatory_options=observatory_context['custom_options'],
        builtin_observatory_options=observatory_context['builtin_options'],
        selected_observatory_value=observatory_context['selected_value'],
        today_utc = datetime.now(timezone.utc).date(),
        bokeh_version = bokeh_version,
        available_feature_columns = FEATURE_COLUMNS,
        default_feature_plot_columns=(get_saved_feature_plot_columns(current_user.id)
                                      if current_user.is_authenticated
                                      else default_feature_plot_columns()),
    )

@main_blueprint.route('/help', methods=['GET'])
def help():
    return render_template(
        "help.html"
    )

@main_blueprint.route('/contact', methods=['GET'])
def contact():
    return render_template(
        "contact.html"
    )

@main_blueprint.route('/profile')
@login_required
def profile():
    """Show user profile and their favorites."""
    try:
        favs = [f.locus_id for f in Favorite.query.filter_by(user_id=current_user.id).order_by(Favorite.created_at.desc()).all()]
    except Exception:
        logger.exception('Failed to load favorites for profile (user_id=%s)', current_user.id)
        favs = []
    return render_template(
        'profile.html', name=current_user.name, favorites=favs
    )

@main_blueprint.route('/download_alerts_csv', methods=['GET'])
@limiter.limit(lambda: current_app.config.get('READ_RATE_LIMIT_MEDIUM', '15 per minute'))
def download_alerts_csv():
    """Download featuretable + classification fields for multiple alert_ids as CSV."""
    alert_ids = [x.strip() for x in request.args.getlist('alert_id') if x and x.strip()]
    if not alert_ids:
        return Response('Missing alert_id', status=400)

    try:
        csv_text, row_count = _build_alerts_csv(alert_ids)
        if row_count == 0:
            return Response('No feature records found for provided alert_ids', status=404)

        response = make_response(csv_text)
        response.headers['Content-Disposition'] = f'attachment; filename="alerts_{row_count}.csv"'
        response.mimetype = 'text/csv'
        return response
    except Exception as e:
        logger.error('Error downloading CSV for alert_ids %s: %s', alert_ids, e, exc_info=True)
        return Response('Error downloading CSV for selected alert ids!', status=500)


def _build_alerts_csv(alert_ids):
    """Build CSV text for alert IDs and return (csv_text, row_count)."""
    unique_alert_ids = list(dict.fromkeys(alert_ids))

    feature_rows = db.session.query(Ztf).filter(Ztf.alert_id.in_(unique_alert_ids)).all()
    if not feature_rows:
        return '', 0

    feature_by_alert_id = {row.alert_id: row for row in feature_rows}
    classification_rows = db.session.query(Classification).filter(Classification.alert_id.in_(unique_alert_ids)).all()
    classification_by_alert_id = {row.alert_id: row for row in classification_rows}

    ordered_feature_rows = [feature_by_alert_id[aid] for aid in unique_alert_ids if aid in feature_by_alert_id]

    first_feature_data = object_as_dict(ordered_feature_rows[0])
    classification_template = object_as_dict(Classification())

    fieldnames = list(first_feature_data.keys()) + [
        col_name for col_name in classification_template.keys() if col_name != 'alert_id'
    ]

    csv_buffer = io.StringIO()
    writer = csv.DictWriter(csv_buffer, fieldnames=fieldnames)
    writer.writeheader()

    for feature_row in ordered_feature_rows:
        merged_row = object_as_dict(feature_row)
        classification_row = classification_by_alert_id.get(feature_row.alert_id)
        if classification_row is not None:
            classification_data = object_as_dict(classification_row)
        else:
            classification_data = classification_template

        for col_name in classification_template.keys():
            if col_name == 'alert_id':
                continue
            merged_row[col_name] = classification_data.get(col_name)

        writer.writerow(merged_row)

    return csv_buffer.getvalue(), len(ordered_feature_rows)


@main_blueprint.route('/query_crossmatches', methods=['GET'])
@limiter.limit(lambda: current_app.config.get('READ_RATE_LIMIT_LAX', '30 per minute'))
def query_crossmatches():
    """Query crossmatches for a given locus id."""
    locusId = request.args.get('locusId') # ex: locusname="ANT2018fywy2"
    if not locusId:
        return Response('Missing locusId', status=400)

    try:
        #query all Crossmatches records from DB where locus id equals given id
        crossmatches_query = db.session.query(Crossmatches)
        crossmatches_query = crossmatches_query.filter(Crossmatches.locus_id == locusId)
        crossmatches_list = result_to_dict(crossmatches_query.all())

        response = current_app.response_class(
            response=safe_serialize(crossmatches_list), #TODO? array[] with single row when using BootstrapTable?
            status=200,
            mimetype='application/json'
        )
        return response
    except Exception as e:
        logger.error('Error querying crossmatches: %s', e, exc_info=True)
        return Response('Error querying crossmatches for selected locus id!', status=500) # Internal Server Error: TODO:Return JSON error message instead of string?

# Register Jinja filters
@main_blueprint.app_template_filter('astro_filter')
def astro_filter(passband):
    if passband == "g":
        return "g"
    elif passband == "R":
        return "R"
    elif passband == "i":
        return "i"
    else:
        return ""

@main_blueprint.app_template_filter('mag_filter')
def mag_filter(num):
    if num: 
        return round(num,3)
    # else:
    #     return ''

@main_blueprint.app_template_filter('format_mjd_readable')
def format_mjd_readable(value):
    if value is None:
        return ''
    
    try:
        mjd_value = float(value)
        return _format_mjd_cached(mjd_value) # Use Astropy for accurate conversion
    except (TypeError, ValueError, OverflowError):
        return ''


def register_blueprints(app):
    """Register all blueprints with the Flask app."""
    app.register_blueprint(main_blueprint)
    app.register_blueprint(favorites_bp)
    app.register_blueprint(filter_bookmarks_bp)
    app.register_blueprint(visual_query_bp)
    app.register_blueprint(lightcurve_bp)
    app.register_blueprint(features_bp)
    app.register_blueprint(user_observatories_bp)
    app.register_blueprint(user_settings_bp)
