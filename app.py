from flask import (Flask, render_template, abort, jsonify, request, Response,
                   redirect, url_for, logging, make_response, Blueprint)
from flask_login import LoginManager, login_required, current_user
import jinja2
# if app.debug is not True:
import logging
logging.basicConfig(handlers=[logging.FileHandler(filename="app.log", 
                                                 encoding='utf-8', mode='a+')],
                    format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                    level=logging.INFO)
from astropy.time import Time
from astropy.coordinates import EarthLocation
import re
import json
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import inspect
#import orm
from . import db
from .models import Ztf, Crossmatches

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.models import Legend

#debug/set_trace breakpoints
import pdb
from pprint import pprint

from .helpers import object_as_dict, result_to_dict

main_blueprint = Blueprint('main', __name__)


@main_blueprint.route('/', methods=['GET', 'POST'])
def start():
    #app.logger.info('Info')
    #app.logger.warning('Warn')
    #logging.error('Exception occurred', exc_info=True) #or: logging.exception()
    logging.info('Request with request_args:' + json.dumps(request.args))
    #pdb.set_trace()
    page = request.args.get('page', 1, type=int)

    data = []
    filter_warning_message = ''
    if request.method == 'GET':
        query = db.session.query(Ztf)
        # Return alerts with a brightness greater than the given value. Ex: ?magpsf=17,18 (range:17-18)
        
        if request.args.get('date_alert_mjd'):
            date_input = extract_numbers(request.args.get('date_alert_mjd'))
            if date_input != None:
                query = extract_float_filter(date_input, Ztf.date_alert_mjd, query)
            else:
                filter_warning_message += 'Date filter cannot be applied - Enter a valid 8-digit integer date of the form yyyymmdd, e.g. "20201207", or range, e.g., "20201207 20201209". You can filter the columns by entering values and then click the "Filter" button.'

        if request.args.get('alert_id'):
            candid_input = extract_numbers(request.args.get('alert_id'))
            if candid_input != None:
                query = query.filter(Ztf.alert_id == int(request.args.get('alert_id')))
            else:
                 filter_warning_message += 'Candid filter cannot be applied - Enter a valid integer, e.g. "1436374650315010006". You can filter the columns by entering values and then click the "Filter" button.'

        if request.args.get('ztf_object_id'):
            query = query.filter(Ztf.ztf_object_id == request.args.get('ztf_object_id'))

        #if request.args.get('jd'):
            #jd_input = extract_numbers(request.args.get('jd'))
            #if jd_input != None:
                #query = extract_float_filter(jd_input, Ztf.jd, query)
            #else:
                #filter_warning_message += 'Jd filter cannot be applied - Enter a valid number, e.g., "2459190.8746528", or range, e.g., "2459190.84 2459190.86". You can filter the columns by entering values and then click the "Filter" button.'

        #if request.args.get('filter'):
            #query = query.filter(Ztf.filter == int(request.args.get('filter'))) # 1:g, 2:r, 3:i
        if request.args.get('ant_passband'):
            query = query.filter(Ztf.ant_passband == request.args.get('ant_passband')) # g, R, i


        if request.args.get('locus_id'):
            query = query.filter(Ztf.locus_id == request.args.get('locus_id'))

        if request.args.get('locus_ra'):
            ra_input = extract_numbers(request.args.get('locus_ra'))
            if ra_input != None:
                query = extract_float_filter(ra_input, Ztf.locus_ra, query)
            else:
                filter_warning_message += 'Ra filter cannot be applied - Enter a valid number, e.g., "118.61421", or range, e.g., "80 90". You can filter the columns by entering values and then click the "Filter" button.'

        if request.args.get('locus_dec'):
            dec_input = extract_numbers(request.args.get('locus_dec'))
            if dec_input != None:
                query = extract_float_filter(dec_input, Ztf.locus_dec, query)
            else:
                filter_warning_message += 'Dec filter cannot be applied - Enter a valid number, e.g., "-20.02131", or range, e.g., "18.8 19.4". You can filter the columns by entering values and then click the "Filter" button.'

        #if request.args.get('magpsf'):
            #magpsf_input = extract_numbers(request.args.get('magpsf'))
            #if magpsf_input != None:
                #query = extract_float_filter(magpsf_input, Ztf.magpsf, query)
            #else:
                #filter_warning_message += 'Magpsf filter cannot be applied - Enter a valid number, e.g., "18.84", or range, e.g., "18.8 19.4". You can filter the columns by entering values and then click the "Filter" button.'

        #Sort order by date (still sorts by mjd column)
        if request.args.get('sort__date'):
            sort__date_order = request.args.get('sort__date')
            if sort__date_order == 'desc':
                query = query.order_by(Ztf.date_alert_mjd.desc())
            if sort__date_order == 'asc':
                query = query.order_by(Ztf.date_alert_mjd.asc())
        else:
            query = query.order_by(Ztf.date_alert_mjd.desc()) #default sort order

        #Sort order by candid
        if request.args.get('sort__candid'):
            sort__candid_order = request.args.get('sort__candid')
            if sort__candid_order == 'desc':
                query = query.order_by(Ztf.alert_id.desc())
            if sort__candid_order == 'asc':
                query = query.order_by(Ztf.alert_id.asc())
        
        #Sort order by objectId
        if request.args.get('sort__objectId'):
            sort__objectId_order = request.args.get('sort__objectId')
            if sort__objectId_order == 'desc':
                query = query.order_by(Ztf.ztf_object_id.desc())
            if sort__objectId_order == 'asc':
                query = query.order_by(Ztf.ztf_object_id.asc())

        #Sort order by jd
        #if request.args.get('sort__jd'):
            #sort__jd_order = request.args.get('sort__jd')
            #if sort__jd_order == 'desc':
                #query = query.order_by(Ztf.jd.desc())
            #if sort__jd_order == 'asc':
                #query = query.order_by(Ztf.jd.asc())

        #Sort order by ra
        if request.args.get('sort__ra'):
            sort__ra_order = request.args.get('sort__ra')
            if sort__ra_order == 'desc':
                query = query.order_by(Ztf.locus_ra.desc())
            if sort__ra_order == 'asc':
                query = query.order_by(Ztf.locus_ra.asc())

        #Sort order by dec
        if request.args.get('sort__dec'):
            sort__dec_order = request.args.get('sort__dec')
            if sort__dec_order == 'desc':
                query = query.order_by(Ztf.locus_dec.desc())
            if sort__dec_order == 'asc':
                query = query.order_by(Ztf.locus_dec.asc())

        ##Sort order by magpsf
        #if request.args.get('sort__magpsf'):
            #sort__magpsf_order = request.args.get('sort__magpsf')
            #if sort__magpsf_order == 'desc':
                #query = query.order_by(Ztf.magpsf.desc())
            #if sort__magpsf_order == 'asc':
                #query = query.order_by(Ztf.magpsf.asc())


        #latest = db.session.query(Ztf).order_by(Ztf.jd.desc()).first() # ? to show latest update date
        #paginator = query.paginate(page, 100, True)
        paginator = query.paginate(page=page, per_page=100, error_out=True)

        #pdb.set_trace()
        # response = {
        #     'has_next': paginator.has_next,
        #     'has_prev': paginator.has_prev,
        #     'results': Alert.serialize_list(paginator.items)
        # }
        #pdb.set_trace()
        print(request.query_string.decode('ascii'))
        print(re.sub('[&?]page=\\d+', '', request.query_string.decode('ascii')))

        site_names = EarthLocation.get_site_names()
        current_date = Time.now().datetime.date()

    return render_template(
        "main.html",
        total_queries=paginator.total,
        table=paginator.items,
        page=paginator.page,
        has_next=paginator.has_next,
        last_page=paginator.pages,
        query_string=re.sub('[&?]?page=\\d+', '', request.query_string.decode('ascii')), # ? b'' binary string
        filter_warning = filter_warning_message,
        observatories = site_names,
        today_utc = current_date
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
    return render_template(
        'profile.html', name=current_user.name
    )

@main_blueprint.route('/query_lightcurve_data', methods=['GET'])
def query_lightcurve_data():
    locusId = request.args.get('locusId')

    #query where locus id equals selected id
    lightcurve_query = db.session.query(Ztf)
    lightcurve_query = lightcurve_query.filter(Ztf.locus_id == locusId) #TODO: load only specific columns
    data = lightcurve_query.all()
    #fields = ['id', 'date_alert_mjd', 'ant_mag_corrected']
    #data = lightcurve_query.options(db.load_only(Ztf.id, Ztf.date_alert_mjd, Ztf.ant_mag_corrected)).all()

    # Creating Plot Figure
    p = figure(height=350, sizing_mode="stretch_width") 
    p.xaxis.axis_label = 'date_alert_mjd'
    p.yaxis.axis_label = 'ant_mag_corrected'
    
    # Defining Plot to be a Scatter Plot
    x_coords = []
    y_coords = []
    #prepare x and y coordinates of rows matched by locus_id of the g band
    for row in (item for item in data if item.ant_passband == 'i'):
        if (row.date_alert_mjd != None and row.ant_mag_corrected != None):
            x_coords += [row.date_alert_mjd]
            y_coords += [row.ant_mag_corrected]
    p.scatter(
        # [i for i in range(10)],
        # [random.randint(1, 50) for j in range(10)],
        x_coords,
        y_coords,
        size=5,
        color="darkkhaki",
        alpha=0.5
    )
    x_coords = []
    y_coords = []
    #prepare x and y coordinates of rows matched by locus_id of the g band
    for row in (item for item in data if item.ant_passband == 'R'):
        if (row.date_alert_mjd != None and row.ant_mag_corrected != None):
            x_coords += [row.date_alert_mjd]
            y_coords += [row.ant_mag_corrected]
    p.scatter(
        # [i for i in range(10)],
        # [random.randint(1, 50) for j in range(10)],
        x_coords,
        y_coords,
        size=5,
        color="indianred",
        alpha=0.5
    )
    x_coords = []
    y_coords = []
    #prepare x and y coordinates of rows matched by locus_id of the g band
    for row in (item for item in data if item.ant_passband == 'g'):
        if (row.date_alert_mjd != None and row.ant_mag_corrected != None):
            x_coords += [row.date_alert_mjd]
            y_coords += [row.ant_mag_corrected]
    p.scatter(
        # [i for i in range(10)],
        # [random.randint(1, 50) for j in range(10)],
        x_coords,
        y_coords,
        size=5,
        color="limegreen",
        alpha=0.5
    )
 
    # Get Chart Components
    script, div = components(p)
 
    # Return the components to the HTML template
    return f'{ div }{ script }'

    return f'''
    <html lang="en">
        <head>
            <script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.3.3.min.js"></script>
            <title>Bokeh Charts</title>
        </head>
        <body>
            <h1>Add Graphs to Flask apps using Python library - Bokeh</h1>
            { div }
            { script }
        </body>
    </html>
    '''

@main_blueprint.route("/locus_plot_csv")
def get_locus_plot():
    locusId = request.args.get('locusId')
    #query where locus id equals selected id
    lightcurve_query = db.session.query(Ztf)
    lightcurve_query = lightcurve_query.filter(Ztf.locus_id == locusId) #TODO: load only specific columns
    data = lightcurve_query.all()
    csv = 'locus_id,date_alert_mjd,ant_mag_corrected\n'
    #prepare x and y coordinates of rows matched by locus_id
    for row in data:
        if (row.date_alert_mjd != None and row.ant_mag_corrected != None):
            csv += f'{row.locus_id},{row.date_alert_mjd},{row.ant_mag_corrected}\n'
    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition":
                 "attachment; filename=myplot.csv"})


#NOTE: Features are filter by objectId & candid
@main_blueprint.route('/query_features', methods=['GET'])
def query_features():
    #objectId = request.args.get('objectId')
    #objectId = request.args.get('ztf_object_id')
    #candid = request.args.get('candid')
    alert_id = request.args.get('alert_id')
    #print('objectId: %s; candid: %s' % ('', alert_id))

    feature_query = db.session.query(Ztf)
    #feature_query = feature_query.filter(Ztf.objectId == objectId)
    feature_query = feature_query.filter(Ztf.alert_id == 'ztf_candidate:'+alert_id)
    data = object_as_dict(feature_query.first())
    #print(data)
    
    response = app.response_class(
        response=json.dumps(data), #TODO? array[] with single row when using BootstrapTable?
        status=200,
        mimetype='application/json'
    )
    return response


@main_blueprint.route('/query_featureplot_data', methods=['GET'])
def query_featureplot_data():
    locusId = request.args.get('locusId')
    selected_features = request.args.get('features')
    #default selected feature list
    feature_list = ['feature_amplitude_magn_r',
        'feature_anderson_darling_normal_magn_r',
        'feature_beyond_1_std_magn_r',
        'feature_beyond_2_std_magn_r',
        'feature_cusum_magn_r']

    if selected_features:
        features_array = selected_features.split(',')
        if features_array and len(features_array) > 0:
            feature_list = features_array
            #print(feature_list)

    #query where locus id equals selected id
    featureplot_query = db.session.query(Ztf)
    featureplot_query = featureplot_query.filter(Ztf.locus_id == locusId) #TODO: load only specific columns
    data = featureplot_query.all()
    #fields = ['id', 'date_alert_mjd', 'ant_mag_corrected']
    #data = lightcurve_query.options(db.load_only(Ztf.id, Ztf.date_alert_mjd, Ztf.ant_mag_corrected)).all()

    # Creating Plot Figure
    p = figure(height=350, sizing_mode="stretch_width") 
    p.xaxis.axis_label = 'date_alert_mjd'
    p.yaxis.axis_label = 'features'
    # Defining Plot to be a Scatter Plot
    x_coords = []
    y_coords = []
    #TODO Features: Retrieve from args.get() - if none: select these as defaults

    #10! colors
    data_colors = [
        "#1f77b4",  # blue
        "#ff7f0e",  # orange
        "#2ca02c",  # green
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#bcbd22",  # yellow-green
        "#17becf"   # cyan
    ]
    #prepare x and y coordinates of rows matched by locus_id
    legend_it = []
    for index, feature in enumerate(feature_list):
        x_coords = []
        y_coords = []
        for row in data:
            if (row.date_alert_mjd != None and row.ant_mag_corrected != None):
                x_coords += [row.date_alert_mjd]
                value = getattr(row, feature)
                #print([getattr(row, feature)])
                y_coords += [value]
        # print(x_coords)
        # print(y_coords)
        pc = p.scatter(
            # [i for i in range(10)],
            # [random.randint(1, 50) for j in range(10)],
            x_coords,
            y_coords,
            size=5,
            color=data_colors[index],
            alpha=0.5
            #,legend_label=feature
        )
        legend_it.append((feature, [pc]))

    legend = Legend(items=legend_it)  #, glyph_height=9, glyph_width=9)
    #legend.click_policy="mute"
    
    # Increasing the glyph height
    # p.legend.glyph_height = 5    
    # # increasing the glyph width
    # p.legend.glyph_width = 5   
    # # Increasing the glyph's label height
    # p.legend.label_height = 5 
    # # Increasing the glyph's label height
    # p.legend.text_font_size = '6px'
    p.add_layout(legend, 'right')
    p.legend.label_text_font_size = '9px'
    p.legend.glyph_width = 12
 
    # Get Chart Components
    script, div = components(p)
    #print(div)
 
    # Return the components to the HTML template
    return f'{ div }{ script }'

@main_blueprint.route('/query_crossmatches', methods=['GET'])
def query_crossmatches():
    """Query crossmatches for a given locus id."""
    locusId = request.args.get('locusId') # locusname="ANT2018fywy2"
    print(locusId)

    try:
        #query all Crossmatches records from DB where locus id equals given id
        crossmatches_query = db.session.query(Crossmatches)
        crossmatches_query = crossmatches_query.filter(Crossmatches.locus_id == locusId)
        dict = result_to_dict(crossmatches_query.all())
        print(dict)

        response = app.response_class(
            response=safe_serialize(dict, locusId), #TODO? array[] with single row when using BootstrapTable?
            status=200,
            mimetype='application/json'
        )
        return response
    except Exception as e:
        logging.error(f"Error querying crossmatches: {e}", exc_info=True)
        return Response(f"{e}", status=500) # Internal Server Error


#if __name__ == '__main__':
#    app.run(debug=True) #uses flask debugger
#    app.run(use_debugger=False, use_reloader=False, passthrough_errors=True) #disable flask debuger and use external debugger instead


