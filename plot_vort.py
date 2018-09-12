import pandas as pd
import shapely.geometry
from ast import literal_eval
import xarray as xr
import tracpy
import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean.cm as cmo
import matplotlib as mpl
import numpy as np

loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
ds = xr.open_dataset(loc)
grid_file = 'gom03_grd_N050_new.nc'
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid(grid_file, proj=proj)

land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-85.0)
extent = [-98, -80, 18, 31]
f = 6.8388546701932663e-05

# for ind in pts.index[0:1]:
for imodel in range(len(ds['ocean_time'])):

    # imodel=0
    t = pd.Timestamp(ds['ocean_time'][imodel].values)

    # imodel = np.where([ds['ocean_time'] == np.datetime64(ind)])[1][0]
    v = ds['v'][imodel].data  # v grid
    u = ds['u'][imodel].data  # u grid

    vort = (v[:,1:] - v[:,:-1])/(grid.x_v[:,1:] - grid.x_v[:,:-1]) \
            - (u[1:,:] - u[:-1,:])/(grid.y_u[1:,:] - grid.y_u[:-1,:])

    up = tracpy.op.resize(u, 0)  # psi grid
    vp = tracpy.op.resize(v, 1)  # psi grid

    fig = plt.figure(figsize=(8,6))# (9.4, 7.7))
    ax = fig.add_subplot(111, projection=merc)
    # ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=merc)
    ax.set_extent(extent, pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.FixedLocator(np.arange(-100, -70, 2))
    gl.ylocator = mpl.ticker.FixedLocator(np.arange(10, 31, 2))
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    # add background land
    ax.add_feature(land_50m, facecolor='0.8')

    vort[vort<0] = 0

    mappable = ax.pcolormesh(grid.lon_psi, grid.lat_psi, vort/f, transform=pc,
                             cmap=cmo.matter, vmin=0, vmax=0.8)
    cb = fig.colorbar(mappable, shrink=0.8, extend='max')
    cb.set_label('Normalized positive vorticity ($\omega/f$)')

    dd = 15
    ax.quiver(grid.lon_psi[1:-1,1:-1][::dd,::dd],
           grid.lat_psi[1:-1,1:-1][::dd,::dd],
           up[1:-1,1:-1][::dd,::dd],
           vp[1:-1,1:-1][::dd,::dd], transform=pc)

    ax.set_title(t.strftime('%H:00 %b %d, %Y'), fontsize=16)

    fig.savefig('figures/vort/%s.png' % t.isoformat()[:13], bbox_inches='tight', dpi=150)
    plt.close(fig)
