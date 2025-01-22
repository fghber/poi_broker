from flask import Blueprint, render_template, abort, request, make_response

import matplotlib 
matplotlib.use('agg') # OR 'SVG'
import matplotlib.pyplot as plt
from matplotlib import dates
import io
import base64

from bokeh.plotting import figure
from bokeh.models import DatetimeTickFormatter, Range1d, LinearAxis
from bokeh.embed import components
from bokeh.io import export_png
from bokeh.models import Legend

import numpy as np
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
import astropy.units as u
from astropy.wcs import WCS

import matplotlib.pyplot as plt
from astroquery.skyview import SkyView

from timezonefinder import TimezoneFinder
from datetime import datetime
from zoneinfo import ZoneInfo

observing_tool_bp = Blueprint('observing_tool', __name__) # TODO: Maybe add a URL prefix observing_tool/

# @observing_tool_bp.route('/observatories')
# def show():
#     site_names = EarthLocation.get_site_names()

# TODO: Add param for chosen observatory from list!
@observing_tool_bp.route('/query_observing_plot')
def calc_observing_plot():
    
    #http://localhost:5000/query_observing_plot?obs_loc=Rubin%20Observatory&obs_date=2024-11-14&dobs_tz=option_utc&ra=101.28715533&dec=16.71611586
    #/query_observing_data?obs_loc=ALMA&obs_date=2025-01-13&dobs_tz=option_utc&ra=15.7574139&dec=16.215272799999994

    obs_loc = request.args.get('obs_loc') #EarthLocation.of_site('Rubin Observatory')
    obs_date = request.args.get('obs_date')
    year, month, day = obs_date.split("-")
    obs_tz = request.args.get('obs_tz')
    ra = float(request.args.get('ra'))
    dec = float(request.args.get('dec'))

    # Observatory location
    observatory = EarthLocation.of_site(obs_loc) 
    obs_lon, obs_lat = observatory.lon.value, observatory.lat.value
    #0. do not generate any plots if the object is not visible from the observatory
    if (obs_lat - dec >=90):
        return "<p>Object is not visible from your location: declination = {dec} degree, obervatory latitude {lat} degree</p>" 

    tf = TimezoneFinder()
    tz = ZoneInfo(tf.timezone_at(lng=obs_lon, lat=obs_lat))

    #stellar_object = SkyCoord(ra=101.28715533*u.deg, dec=16.71611586*u.deg)
    stellar_object = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

    # Time settings
    midnight_utc = Time(f'{year}-{month}-{day} 00:00:00', scale='utc') #values are already in the correct format no :02d etc needed
    #t_utc = Time(midnight_utc, format='isot', scale='utc') # ? no use
    midnight_zone = midnight_utc.to_datetime(timezone=tz)
    delta_midnight = np.linspace(-12, 12, 1000) * u.hour
    times_range_utc = midnight_utc + delta_midnight
    times_range_zone = midnight_utc + delta_midnight + (tz.utcoffset(midnight_zone).total_seconds() / (60*60)  )*u.hour #potential bug?

    frame = AltAz(obstime=times_range_utc, location=observatory)
    
    object_altazs = stellar_object.transform_to(frame)
    #object_airmass = object_altazs.secz # ? Not used

    moon = get_body("moon", times_range_utc, location=observatory)
    moon_altazs = moon.transform_to(frame)
    moon_alt = moon_altazs.alt.value

    sun = get_body("sun", times_range_utc, location=observatory)
    sun_altazs = sun.transform_to(frame)
    sun_alt = sun_altazs.alt.value

    # (Nested) Observing Plot function (to structure code better)
    def observing_plot():
        plt.style.use(astropy_mpl_style)
        plt.figure(figsize=(8,6.5))
        #plt.margins(0.01, tight=True) #The default margins are rcParams["axes.xmargin"] (default: 0.05) and rcParams["axes.ymargin"] (default: 0.05).
        quantity_support()
        ax = plt.gca()
        if(obs_tz == 'option_utc'):
            timetoplot = times_range_utc
            ax.set_xlabel("Time starting {0} UTC]".format(min(timetoplot).datetime.date()))
        else:
            timetoplot = times_range_zone
            utcoffset = tz.utcoffset(midnight_zone).total_seconds() / (60*60)
            ax.set_xlabel("Time starting {0} [{1}, UTC{2}]".format(min(timetoplot).datetime.date(), tz,f'{utcoffset:+.0f}'))

        # Format the time axis
        xlo, xhi = (timetoplot[0]), (timetoplot[-1])
        ax.set_xlim([xlo.plot_date, xhi.plot_date])
        date_formatter = dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(date_formatter)
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

        plt.fill_between(
            timetoplot.datetime,
            0 * u.deg,
            90 * u.deg,
            sun_altazs.alt < -0 * u.deg,
            color="0.5",
            zorder=0,
        )

        plt.fill_between(
        timetoplot.datetime,
        0 * u.deg,
        90 * u.deg,
        sun_altazs.alt < -18 * u.deg,
        color="k",
        zorder=0,
        )

        plt.scatter(
            timetoplot.datetime,
            moon_altazs.alt.value,s=1,c="lightblue")

        plt.scatter(
            timetoplot.datetime,
            object_altazs.alt.value,s=1,c="orange")

        plt.ylim(0, )
        plt.ylabel("Altitude [deg]")
        ax.set_ylim(0,90)
        airmass_ticks = np.array([1, 2, 3])
        altitude_ticks = 90 - np.degrees(np.arccos(1/airmass_ticks))

        ax2 = ax.twinx()
        ax2.set_yticks(altitude_ticks)
        ax2.set_yticklabels(airmass_ticks)
        ax2.set_ylim(ax.get_ylim())
        ax2.set_ylabel('Airmass')
        plt.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
        print(altitude_ticks)
        ax2.grid(None)

        #plt.show()
        # Create an in-memory buffer
        img_io = io.BytesIO()
        plt.savefig(img_io, format='png')
        img_io.seek(0)

        # Option a: Create a response with the image data
        #response = make_response(img_io.read())
        #response.headers['Content-Type'] = 'image/png'
        # Option b : Encode image to base64
        img_data = base64.b64encode(img_io.getvalue()).decode('utf-8')
        img_obs = f"data:image/png;base64,{img_data}"
        # Close plot
        plt.close()
        return img_obs

    #1: Create an observing plot
    obs_img = observing_plot()

    #2: Create a finder chart
    finder_img = plot_finder_image(stellar_object)

    #3: Create Moon Panel
    moon_panel = ''
    # 1st: Is the moon up at night?
    night_moon_alt = moon_alt[np.where(sun_alt<0)]
    #print(night_moon_alt)
    if(np.max(night_moon_alt) < 0):
        moon_panel = 'Moon down'
    else: 
        #Moon up
        moon_separation = moon.separation(stellar_object, origin_mismatch="ignore")
        moon_panel = get_moon_phase_panel(observatory, midnight_utc, moon_separation)

    #messier1 = FixedTarget.from_name("M1")
    # ra=101.28715533
    # dec=16.71611586
    # target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

    #return response
    return f'<hr><div class="row"><div class="col-md-6"><img src="{obs_img}"></div><div class="col-md-6">&nbsp;&nbsp;&nbsp;<img src="{finder_img}"></div></div>{moon_panel}'

    #Alternative Observing Plot using Bokeh
    '''
    p = figure(x_axis_type="datetime", title="Altitude vs Time", height=350, width=800)
    p.xaxis.axis_label = "Time"
    p.yaxis.axis_label = "Altitude [deg]"

    if obs_tz == 'option_utc':
        timetoplot = times_range_utc
        p.xaxis.axis_label = "Time starting {0} UTC".format(min(timetoplot).datetime.date())
    else:
        timetoplot = times_range_zone
        utcoffset = tz.utcoffset(midnight_zone).total_seconds() / (60*60)
        p.xaxis.axis_label = "Time starting {0} [{1}, UTC{2}]".format(min(timetoplot).datetime.date(), tz, f'{utcoffset:+.0f}')

    # Format the time axis
    p.xaxis.formatter = DatetimeTickFormatter(hours="%H:%M")
    p.xaxis.major_label_orientation = 30

    x = np.append(timetoplot.datetime, timetoplot.datetime[::-1])
    y = np.append([0] * len(timetoplot.datetime), [90] * len(timetoplot.datetime))
    mask = np.append(sun_altazs.alt < -0 * u.deg, sun_altazs.alt[::-1] < -0 * u.deg)
    p.patch(x[mask], y[mask], color="lightgrey", alpha=0.5, legend_label="Sun Altitude < 0 deg")

    mask = np.append(sun_altazs.alt < -18 * u.deg, sun_altazs.alt[::-1] < -18 * u.deg)
    p.patch(x[mask], y[mask], color="black", alpha=0.5, legend_label="Sun Altitude < -18 deg")

    # Scatter plots for moon and object altitudes
    p.circle(timetoplot.datetime, moon_altazs.alt.value, size=1, color="lightblue", legend_label="Moon Altitude")
    p.circle(timetoplot.datetime, object_altazs.alt.value, size=1, color="orange", legend_label="Object Altitude")

    # Set y-axis limits
    p.y_range = Range1d(0, 90)

    # Airmass ticks
    airmass_ticks = np.array([1, 2, 3])
    altitude_ticks = 90 - np.degrees(np.arccos(1/airmass_ticks))

    # Secondary y-axis for airmass
    p.extra_y_ranges = {"airmass": Range1d(start=0, end=90)}
    p.add_layout(LinearAxis(y_range_name="airmass", axis_label="Airmass", major_label_overrides=dict(zip(altitude_ticks, map(str, airmass_ticks)))), 'right')

    # Show grid
    p.grid.grid_line_color = 'grey'
    p.grid.grid_line_dash = 'dashed'
    p.grid.grid_line_width = 0.5

    # Option b: Generate script and div for embedding
    script, div = components(p)
    # Return the components to the HTML template
    #print(f'{ div }{ script }')
    return f'{ div }{ script }'
    '''

def get_moon_phase_panel(observatory, midnight_utc, moon_separation):
    #angle of the tilt of the Moon will be different as seen from different latitudes. 
    #https://astronomy.stackexchange.com/questions/24711/how-does-the-moon-look-like-from-different-latitudes-of-the-earth
    #Calculate lunar orbital phase in radians.

    sun_midnight = get_body("sun", midnight_utc, location=observatory)
    moon_midnight = get_body("moon", midnight_utc, location=observatory)

    elongation = sun_midnight.separation(moon_midnight)
    moon_phase_angle_inc = np.arctan2(sun_midnight.distance*np.sin(elongation),
                moon_midnight.distance - sun_midnight.distance*np.cos(elongation))

    fraction_illuminated = (1 + np.cos(moon_phase_angle_inc))/2.0
    fraction_illuminated_percentage = "{:.0%}".format(fraction_illuminated)
    angle = np.arctan(np.cos(sun_midnight.dec) * np.sin(sun_midnight.ra - moon_midnight.ra), np.sin(sun_midnight.dec) * np.cos(moon_midnight.dec) -
                np.cos(sun_midnight.dec) * np.sin(moon_midnight.dec) * np.cos(sun_midnight.ra - moon_midnight.ra))        
    phase = 0.5 + 0.5 * moon_phase_angle_inc.value * np.sign(angle.value) / np.pi

    phase_name=''
    phase_image=''
    if(phase == 0):
        phase_name='new moon'
        phase_image = 'Moon_new.png'
    if(0 < phase < 0.25):
        phase_name='waxing crescent'
        phase_image = 'Moon_waxingcrescent.png'
    if(phase == 0.25):
        phase_name='first quarter'
        phase_image = 'Moon_firstquarter.png'
    if(0.25 < phase < 0.5):
        phase_name='waxing gibbous'
        phase_image = 'Moon_waxinggibbous.png'
    if(phase == 0.5):
        phase_name='full'
        phase_image = 'Moon_full.png'
    if(0.5 < phase < 0.75):
        phase_name='waning gibbous'
        phase_image = 'Moon_waninggibbous.png'
    if(phase == 0.75):
        phase_name='last quarter'
        phase_image = 'Moon_lastquarter.png'
    if(0.75 < phase < 1):
        phase_name='waning crescent'
        phase_image = 'Moon_waningcrescent.png'
    if(phase == 1):
        phase_name='new moon'
        phase_image = 'Moon_new.png'

    #Rotate the moon picture counterclockwise by (90 - lat_observatory).
    rotation = (90 - observatory.lat.value);

    html = '<hr><div class="row"><div class="col-md-6">'
    html += f'Moon Phase at {midnight_utc} UTC<br>' #TODO:  UTC with 2 digit after ??
    html += f'Phase: {phase_name}<br>'
    html += f'Illumination: {fraction_illuminated_percentage}<br>'
    html += f'separation from moon to object during night: {np.min(moon_separation.arcminute):.3f} to {np.max(moon_separation.arcminute):.3f}' #TODO: round to 3 digits okay?
    html += '</div><div class="col-md-6"><span class="moon-container-square">'
    html += f'<img src="/static/img/{phase_image}" width="96" height="96" style="transform: rotate({rotation}deg);">'
    html += '</span></div></div>'
    #TODO: tilt based on latitude, where I show it on a larger black square
    #this is how it should look like: https://astronomy.stackexchange.com/questions/24711/how-does-the-moon-look-like-from-different-latitudes-of-the-earth
    return html


def plot_finder_image(coord, survey='DSS', fov_radius=10*u.arcmin,
                      log=False, ax=None, grid=False, reticle=False,
                      style_kwargs=None, reticle_style_kwargs=None):
    """
    Plot survey image centered on ``target``.

    Survey images are retrieved from NASA Goddard's SkyView service via
    ``astroquery.skyview.SkyView``.

    If a `~matplotlib.axes.Axes` object already exists, plots the finder image
    on top. Otherwise, creates a new `~matplotlib.axes.Axes`
    object with the finder image.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`
        Coordinates of celestial object

    survey : string
        Name of survey to retrieve image from. For dictionary of
        available surveys, use
        ``from astroquery.skyview import SkyView; SkyView.list_surveys()``.
        Defaults to ``'DSS'``, the Digital Sky Survey.

    fov_radius : `~astropy.units.Quantity`
        Radius of field of view of retrieved image. Defaults to 10 arcmin.

    log : bool, optional
        Take the natural logarithm of the FITS image if `True`.
        False by default.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The `~matplotlib.axes.Axes` object to be drawn on.
        If None, uses the current `~matplotlib.axes.Axes`.

    grid : bool, optional.
        Grid is drawn if `True`. `False` by default.

    reticle : bool, optional
        Draw reticle on the center of the FOV if `True`. Default is `False`.

    style_kwargs : dict or `None`, optional.
        A dictionary of keywords passed into `~matplotlib.pyplot.imshow`
        to set plotting styles.

    reticle_style_kwargs : dict or `None`, optional
        A dictionary of keywords passed into `~matplotlib.pyplot.axvline` and
        `~matplotlib.pyplot.axhline` to set reticle style.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        Matplotlib axes with survey image centered on ``target``

    hdu : `~astropy.io.fits.PrimaryHDU`
        FITS HDU of the retrieved image


    Notes
    -----
    Dependencies:
        In addition to Matplotlib, this function makes use of astroquery.
    """

    #coord = target if not hasattr(target, 'coord') else target.coord
    #print(coord)
    position = coord.icrs
    coordinates = 'icrs'
    #print(position)
    #target_name = None if isinstance(target, SkyCoord) else target.name

    hdu = SkyView.get_images(position=position, coordinates=coordinates,
                             survey=survey, radius=fov_radius)[0][0]
    wcs = WCS(hdu.header)
    
    plt.figure(figsize=(6.5,6.5))
    #plt.margins(10)
    # Set up axes & plot styles if needed.
    if ax is None:
        ax = plt.gcf().add_subplot(projection=wcs) # this makes the coordinates
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('cmap', 'Greys')
    style_kwargs.setdefault('origin', 'lower')

    
    
    lon = ax.coords[0]
    lat = ax.coords[1]

    lon.set_major_formatter('dd:mm:ss.s')
    lat.set_major_formatter('dd:mm')

    
    #other option

    #lon.set_major_formatter('d.d')#('dd:mm:ss.s')
    #lat.set_major_formatter('d.d')#('dd:mm')
    # https://docs.astropy.org/en/stable/visualization/wcsaxes/ticks_labels_grid.html
    
    if log:
        image_data = np.log(hdu.data)
    else:
        image_data = hdu.data
    ax.imshow(image_data, **style_kwargs)

    # Draw reticle
#    if reticle:
#        pixel_width = image_data.shape[0]
#        inner, outer = 0.03, 0.08

#        if reticle_style_kwargs is None:
#            reticle_style_kwargs = {}
#        reticle_style_kwargs.setdefault('linewidth', 2)
#        reticle_style_kwargs.setdefault('color', 'm')

#        ax.axvline(x=0.5*pixel_width, ymin=0.5+inner, ymax=0.5+outer,
#                   **reticle_style_kwargs)
#        ax.axvline(x=0.5*pixel_width, ymin=0.5-inner, ymax=0.5-outer,
#                   **reticle_style_kwargs)
#        ax.axhline(y=0.5*pixel_width, xmin=0.5+inner, xmax=0.5+outer,
#                   **reticle_style_kwargs)
#        ax.axhline(y=0.5*pixel_width, xmin=0.5-inner, xmax=0.5-outer,
#                   **reticle_style_kwargs)

    # Labels, title, grid
    ax.set(xlabel='RA', ylabel='DEC')
    #if target_name is not None:
    #    ax.set_title(target_name)
#    ax.grid(grid)
    
    
    # add marker
    ax.scatter(coord.ra,coord.dec,marker="+",c='r',s=150) ####### for this, use a symbol instead that is a cross that is empty at the center; color red

    # Redraw the figure for interactive sessions.
#    ax.figure.canvas.draw()

    # Create an in-memory buffer
    img_io = io.BytesIO()
    plt.savefig(img_io, format='png')
    img_io.seek(0)

    # Option a: Create a response with the image data
    #response = make_response(img_io.read())
    #response.headers['Content-Type'] = 'image/png'
    # Option b : Encode image to base64
    img_data = base64.b64encode(img_io.getvalue()).decode('utf-8')
    img_finder = f"data:image/png;base64,{img_data}"

    # Close plot
    plt.close()

    return img_finder
