
import numpy as np

import astropy.units as u

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


import matplotlib.pyplot as plt

import astroquery
import matplotlib.pyplot as plt
from astroquery.skyview import SkyView



def plot_finder_image(coord, survey='DSS', fov_radius=10*u.arcmin,
                      log=False, ax=None, grid=False, reticle=False,
                      style_kwargs=None, reticle_style_kwargs=None):
  
    #coord = target if not hasattr(target, 'coord') else target.coord
    print(coord)
    position = coord.icrs
    coordinates = 'icrs'
    print(position)
    #target_name = None if isinstance(target, SkyCoord) else target.name

    print('before')

    hdu = SkyView.get_images(position=position, coordinates=coordinates,
                             survey=survey, radius=fov_radius)[0][0]

    print('after')
    wcs = WCS(hdu.header)
    print('test')

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

    # Labels, title, grid
    ax.set(xlabel='RA', ylabel='DEC')
    #if target_name is not None:
    #    ax.set_title(target_name)
#    ax.grid(grid)
    
    
    # add marker
    ax.scatter(ra,dec,marker="+",c='r') ####### for this, use a symbol instead that is a cross that is empty at the center; color red

    # Redraw the figure for interactive sessions.
#    ax.figure.canvas.draw()
    return ax, hdu

#messier1 = FixedTarget.from_name("M1")

ra=101.28715533
dec=16.71611586

#from astroquery.skyview import SkyView
#SkyView.clear_cache()

#If this function is unavailable, upgrade your version of astroquery. The clear_cache function was introduced in version 0.4.7.dev8479.

target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

ax, hdu = plot_finder_image(target)


plt.show()    
