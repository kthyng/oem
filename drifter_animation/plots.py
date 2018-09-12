import pandas as pd
import shapely.geometry
import xarray as xr
import tracpy
import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean.cm as cmo
import matplotlib as mpl
import numpy as np
import os

# load in gulf ssh
loc = 'output/*.nc'
ds = xr.open_mfdataset(loc)
lats = ds['latitude']
lons = ds['longitude']
# grid_file = '../gom03_grd_N050_new.nc'
# grid = xr.open_dataset(grid_file)
# proj = tracpy.tools.make_proj('nwgom-pyproj')
# grid = tracpy.inout.readgrid(grid_file, proj=proj)

# load in drifter data
dfmix = pd.read_csv('drifter-09_04_18-07_37.csv', usecols=[0,1,4,5], parse_dates=True, index_col=1)
dfmix = dfmix[::-1]
drifters = list(dfmix['DeviceName'].unique())  # drifter names

# interpolate all drifters to standard times
# each df in dfs is a drifter lon/lat
# start when drifter is near eddy, but at hour
dstart = '2018-8-12 18:00'  # dfmix.index[(dfmix['Longitude'] < -90) * (dfmix['Latitude'] < 28)][0].floor('H')
indices = pd.date_range(start=dstart, end=dfmix.index.max(), freq='6H')
dfs = []
for drifter in drifters:
    df = pd.DataFrame(index=indices)
    df = dfmix[dfmix['DeviceName']==drifter].resample('1H', base=0).mean()
    dfs.append(df)

# plot prep
land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-85.0)
extent = [-98, -85, 22, 31]

# loop over drifter times
for index in indices:

    fname = 'figures/%s.png' % (index.isoformat()[:13])

    if os.path.exists(fname):
        continue

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
    ax.add_feature(land_50m, facecolor='0.8', edgecolor='0.2', linewidth=1)

    # use nearest model output
    ssh = ds['zos'].sel(time=index, method='nearest')
    u = ds['uo'].sel(time=index, method='nearest').isel(depth=0)
    v = ds['vo'].sel(time=index, method='nearest').isel(depth=0)

    mappable = ax.pcolormesh(lons, lats, ssh, transform=pc,
                             cmap=cmo.balance, vmin=-0.5, vmax=0.5)
    cb = fig.colorbar(mappable, shrink=0.8, extend='both')
    cb.set_label('Sea surface height [m]')

    # overlay velocities
    dd = 5
    ax.quiver(lons[::dd].values, lats[::dd].values,
           u[::dd,::dd].values, v[::dd,::dd].values, transform=pc, alpha=0.6)


    # overlay drifters
    for i in range(3):
        try:
            lat, lon = dfs[i].loc[index]
            ax.plot(lon, lat, 'ok', transform=pc)
        except:
            pass

    ax.set_title(index.strftime('%b %d %H:%M, %Y'))

    # add labels
    # drifter legend
    ax2 = fig.add_axes([0.125, 0.625, 0.15, 0.2], frameon=False)
    ax2.scatter([], [], c='0.2', s=30, marker='o', label='Drifter')
    ax2.legend(scatterpoints=1, frameon=False, loc='upper left')
    ax2.set_axis_off()
    # label mercator
    ax.text(0.3, 0.96, 'Mercator model output', color='0.2', transform=ax.transAxes)
    # my name
    ax.text(0.8, 0.96, 'K. Thyng', color='0.2', transform=ax.transAxes)


    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
