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
loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
ds = xr.open_dataset(loc)

pts = pd.read_csv('points.csv', parse_dates=True, index_col=0, header=0)#, nrows=10)

# remove nan's
pts = pts[~pts['x'].isnull()]
pts = pts[~pts['y'].isnull()]

pts.x = pts.x.apply(literal_eval)
pts.y = pts.y.apply(literal_eval)

# for ind in pts.index:
ind = pts.index[0]

for x, y in zip(pts.loc[ind,'x'],pts.loc[ind,'y']):
    plt.plot(x,y)
# x, y = pts.loc[ind,'x'][0],pts.loc[ind,'y'][0]
# lon, lat, _ = tracpy.tools.interpolate2d(np.asarray(x), np.asarray(y), grid, 'm_ij2ll')
# p = shapely.geometry.Polygon(zip(x,y))



# for how long can I find a polygon that intersects in subsequent time steps?





# Make an animation of quiver arrows on map with polygon identifying cyclones

land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-85.0)
extent = [-98, -80, 18, 31]
colors = ['r', 'orange', 'g', 'b', 'purple']

for ind in pts.index[10:20]:
# ind = pts.index[0]


    imodel = np.where([ds['ocean_time'] == np.datetime64(ind)])[1][0]
    v = ds['v'][imodel].data  # v grid
    u = ds['u'][imodel].data  # u grid
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
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    # add background land
    ax.add_feature(land_50m, facecolor='0.8')

    ax.quiver(grid.lon_psi[1:-1,1:-1][::10,::10],
           grid.lat_psi[1:-1,1:-1][::10,::10],
           up[1:-1,1:-1][::10,::10],
           vp[1:-1,1:-1][::10,::10], transform=pc)


    for x, y, color in zip(pts.loc[ind,'x'],pts.loc[ind,'y'], colors):
    # x, y = pts.loc[ind,'x'][0],pts.loc[ind,'y'][0]
        lon, lat = proj(x, y, inverse=True)

        ax.plot(lon, lat, color=color, lw=2, transform=pc)

    fig.savefig('figures/example/%s.png' % ind.isoformat(), bbox_inches='tight')
    plt.close(fig)
