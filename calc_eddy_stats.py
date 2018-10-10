'''
Calculate eddy statistics from identified eddies from swirl.py.
'''

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
from pyproj import Proj
from tracpy import op
import scipy.ndimage
import ast

aed = cartopy.crs.AzimuthalEquidistant(central_longitude=-90.0, central_latitude=26)
pc = cartopy.crs.PlateCarree()
land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
merc = cartopy.crs.Mercator(central_longitude=-90.0)
extent = [-98, -80, 18, 31]


# read in eddies identified in swirl.py
df = pd.read_csv('calcs/swirl/points.csv', index_col=0, parse_dates=True)

# initialize column to hold polygons
df['p'] = ''
df['p'] = df['p'].apply(list)

# reconstituting polygon shapes and adding them into appropriate time as list
for date in df.index:
    xs = ast.literal_eval(df.loc[date, 'x'])
    ys = ast.literal_eval(df.loc[date, 'y'])

    pstemp = [shapely.geometry.Polygon(np.vstack((x, y)).T) for (x,y) in zip(xs, ys)]

    df.loc[date,'p'] = pstemp

## filter out polygons. make a new column for each identified set of eddies.


# initialize dataframe where each column is an eddy being tracked
dfeddies = pd.DataFrame(index=df.index)

# initialize column to hold continuous polygon
pcount = 0  # counter for polygon columns
colname = 'p%i' % pcount
dfeddies[colname] = ''
dfeddies[colname] = dfeddies[colname].apply(list)
pcount += 1

# loop through dates, then polygons, then dates again, to find continuous eddy
# identifications
for i, date0 in enumerate(dfeddies.index):  # loop over dates
# i=0; date0=df.index[0]
    # import pdb; pdb.set_trace()
    ps0 = df.loc[date0,'p']  # polygons from date0
    for p0 in ps0:  # loop over polygons from date0
    # p0=ps0[0]
        if isinstance(p0, float) and np.isnan(p0):  # this was previously nan'ed out
            continue
        # flag for if I have added p0 to the new polygon column or not
        addp0flag = False
        for j, date1 in enumerate(dfeddies.index[i+1:]):  # loop over subsequent dates
            # date1=df.index[i+1]
            # if i == 1:
            #     import pdb; pdb.set_trace()

            ps1 = df.loc[date1,'p']
            # flag for whether there is an intersection with another eddy for this date
            interflag = False
            # p1=ps1[0]
            for k, p1 in enumerate(ps1):  # loop over polygons from date1
                if isinstance(p1, float) and np.isnan(p1):  # this was previously nan'ed out
                    continue
                if p0.intersects(p1):
                    # add to list, but only if hasn't already been added
                    if not addp0flag:
                        # nan out between 0 and i
                        dfeddies[colname].iloc[0:i] = np.nan
                        dfeddies.loc[date0,colname] = p0
                        addp0flag = True
                        # add in centroid coords and diameter
                        lon, lat = pc.transform_point(p0.centroid.xy[0][0], p0.centroid.xy[1][0], aed)
                        dfeddies.loc[date0,colname + ' lon [centroid]'] = lon
                        dfeddies.loc[date0,colname + ' lat [centroid]'] = lat
                        # calculate equivalent radius for polygon area (km)
                        R = np.sqrt(p0.area/1000**2/np.pi)
                        dfeddies.loc[date0,colname + ' diameter [km]'] = 2*R
                    dfeddies.loc[date1,colname] = p1
                    # add in centroid coords and diameter
                    lon, lat = pc.transform_point(p1.centroid.xy[0][0], p1.centroid.xy[1][0], aed)
                    dfeddies.loc[date1,colname + ' lon [centroid]'] = lon
                    dfeddies.loc[date1,colname + ' lat [centroid]'] = lat
                    # calculate equivalent radius for polygon area (km)
                    R = np.sqrt(p1.area/1000**2/np.pi)
                    dfeddies.loc[date1,colname + ' diameter [km]'] = 2*R
                    interflag = True
                    # replace p0 with p1 since eddy will translate over time
                    p0 = p1
                    # nan out df.loc[date1,'p'][k] (p1) so that it isn't used again on a subsequent loop
                    df.loc[date1,'p'][k] = np.nan

                # should i check for 2 eddies? make multipolygon?

                    # if we have already found a polygon from ps1 to match p0, we
                    # can move on (not checking for multiple eddies currently)
                    break
            # if colname == 'p12':
            #     import pdb; pdb.set_trace()
            # if there was no intersection this time, but was earlier, do these
            # final steps
            if addp0flag and (not interflag or i==len(dfeddies.index)-1):
                # calculate the propagation speed for the times the eddy existed
                lons = dfeddies['%s lon [centroid]' % colname].values
                lats = dfeddies['%s lat [centroid]' % colname].values
                pts = aed.transform_points(pc, lons, lats)
                # distance between subsequent eddy centroids [km]
                dist = np.sqrt((pts[1:,0] - pts[:-1,0])**2 + (pts[1:,1] - pts[:-1,1])**2)/1000
                ind = dfeddies['%s lon [centroid]' % colname].index
                # calculate propagation speed from distance
                dt = (ind[2:] - ind[:-2]).total_seconds()
                # centroid propagation speed in km/day
                speed = abs(dist[1:] - dist[:-1])/dt * (86400)
                dfeddies.loc[2:,colname + ' speed [km/day]'] = speed
                # change all subsequent rows in this column to nans
                dfeddies[colname].loc[date1:] = np.nan
                break
            # if there were no intersections this set of eddies, move onto next polygon p0
            elif not interflag:
                break
        # if this flag is true, there must have been an intersection, so start a new column
        if addp0flag:
            colname = 'p%i' % pcount
            dfeddies[colname] = ''
            dfeddies[colname] = dfeddies[colname].apply(list)
            pcount += 1


# have to remake df because pandas is stupid and won't make an actual copy of
# a dataframe
# read in eddies identified in swirl.py
df = pd.read_csv('calcs/swirl/points.csv', index_col=0, parse_dates=True)

# initialize column to hold polygons
df['p'] = ''
df['p'] = df['p'].apply(list)

# reconstituting polygon shapes and adding them into appropriate time as list
for date in df.index:
    xs = ast.literal_eval(df.loc[date, 'x'])
    ys = ast.literal_eval(df.loc[date, 'y'])

    pstemp = [shapely.geometry.Polygon(np.vstack((x, y)).T) for (x,y) in zip(xs, ys)]

    df.loc[date,'p'] = pstemp


# calculate statistics and make plots based on set of eddies identified

# avg gyre persistence (time), gyre diameter, frequency, timing of eddies,
# propagation speed.
# map and histogram of these numbers

# loop through p_i columns and calculate properties
lastcolnum = int(dfeddies.columns[-1].split(' ')[0][1:])  # int(dfeddies.columns[-1].split('p')[1])
# don't include last p column since should be empty
rownames = ['p%i' % i for i in range(lastcolnum)]
# bool of eddies that are "long enough" to count (more than 7 days currently)
longeddies = dfeddies[rownames].count() > 7
shorteddies = dfeddies[rownames].count() <= 7

# for rowname in
# rowname = 'p0'
metrics = ['persistence [days]', 'mean diameter [km]',
           'mean propagation speed [km/day]']
# index = ['p%i' % i for i in range(int(df.columns[-1].split('p')[1]))]
dfstats = pd.DataFrame(index=rownames, columns=metrics)


# calculate persistance
# don't count eddy times lower than chosen in longeddies
dfstats['persistence [days]'] = dfeddies[rownames].loc[:,longeddies].count()
# print(dfstats['persistence [days]'].mean(), dfstats['persistence [days]'].max())


# calculate mean diameter for each eddy over time
rownames2 = [rowname + ' diameter [km]' for rowname in rownames]  # diameter columns
dfstats['mean diameter [km]'] = dfeddies[rownames2].mean().values
# nan out short time eddies
dfstats['mean diameter [km]'][shorteddies] = np.nan
# print(dfstats['mean diameter [km]'].mean(), dfstats['mean diameter [km]'].max())


# calculate number of eddies per month
# number of days eddy was present for each month of year
numdays = (~dfeddies[rownames].loc[:,longeddies].isnull()).groupby(dfeddies.index.month).sum()
# just counting presence of eddy
freq = (numdays > 0).sum(axis='columns')

# want interannual variability
# number of days eddy present by month and year
s = (~dfeddies[rownames].loc[:,longeddies].isnull()).groupby([(dfeddies.index.year),(dfeddies.index.month)]).sum()
# swap index level
s = s.swaplevel()
# nan out zeros
s[s==0] = np.nan
# standard deviation across years associated with num eddies per month
freqstd = [s.loc[i,:].count(axis=1).std() for i in range(1,13)]


# plot timing of eddies — at what time of year are they present?
# freq.plot()

# calculate mean propagation speed
rownames2 = [rowname + ' speed [km/day]' for rowname in rownames]  # speed columns
dfstats['mean propagation speed [km/day]'] = dfeddies[rownames2].mean().values
# nan out short time eddies
dfstats['mean propagation speed [km/day]'][shorteddies] = np.nan
# print(dfstats['mean propagation speed [km/day]'].mean(), dfstats['mean propagation speed [km/day]'].max())

dfstats.to_csv('calcs/swirl/dfstats.csv')


## print stats/show summary plots ##
sumbase = 'figures/swirl/hycom/withmetrics/summary/'

# timing during the year
nyears = 11  # number of years of output
ax = (freq/nyears).plot(figsize=(10,4), grid=True)
ax.set_xlabel('Months of the year')
ax.set_ylabel('Number of eddies per year\nstandard deviation across 10 years')
ax.fill_between(freq.index, (freq/nyears)-freqstd, (freq/nyears)+freqstd, alpha=0.2, color='k')
ax.set_xlim(1,12)
ax.set_xticks(np.arange(1,13));
ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
plt.savefig(sumbase + 'timeofyear.pdf', bbox_inches='tight')

# diameter
ax = dfstats['mean diameter [km]'].hist(bins=20, grid=False, edgecolor='k')
ax.set_xlabel('Mean eddy diameter across lifetime [km]')
ax.set_ylabel('Number of cyclones longer than 7 days in bins')
# overall mean diameter
dmean = dfstats['mean diameter [km]'].mean()
text0 = 'overall mean diameter: %dkm' % dmean
ax.text(0.4, 0.95, text0, transform=ax.transAxes)
plt.savefig(sumbase + 'diameter.pdf', bbox_inches='tight')

# mean propagation speed
ax = dfstats['mean propagation speed [km/day]'].hist(bins=20, grid=False, edgecolor='k')
ax.set_xlabel('Mean eddy propagation speed across lifetime [km/day]')
ax.set_ylabel('Number of cyclones longer than 7 days in bins')
# overall mean speed
dmean = dfstats['mean propagation speed [km/day]'].mean()
text0 = 'overall mean speed: %dkm/day' % dmean
ax.text(0.4, 0.95, text0, transform=ax.transAxes)
plt.savefig(sumbase + 'speed.pdf', bbox_inches='tight')

# persistence
ax = dfstats['persistence [days]'].hist(bins=20, grid=False, edgecolor='k')
ax.set_xlabel('Persistence of eddy [days]')
ax.set_ylabel('Number of cyclones longer than 7 days in bins')
ylims = ax.get_ylim()
ax.set_ylim(*ylims)
xlims = ax.get_xlim()
ax.set_xticks(np.arange(0,xlims[1], 15))
ax.vlines(30, 0, ylims[1], linewidth=3, alpha=0.3)
ax.vlines(60, 0, ylims[1], linewidth=3, alpha=0.3)
ax.vlines(90, 0, ylims[1], linewidth=3, alpha=0.3)
# number of eddies
neddies0 = (dfstats['persistence [days]']>0).sum()
text0 = '# eddies: %d' % neddies0
ax.text(0.4, 0.95, text0, transform=ax.transAxes)
# over 30 days
neddies30 = (dfstats['persistence [days]']>30).sum()
text30 = '# eddies over 1 month long: %d' % neddies30
ax.text(0.4, 0.9, text30, transform=ax.transAxes)
# over 60 days
neddies60 = (dfstats['persistence [days]']>60).sum()
text60 = '# eddies over 2 months long: %d' % neddies60
ax.text(0.4, 0.85, text60, transform=ax.transAxes)
# over 90 days
neddies90 = (dfstats['persistence [days]']>90).sum()
text90 = '# eddies over 3 months long: %d' % neddies90
ax.text(0.4, 0.8, text90, transform=ax.transAxes)
plt.savefig(sumbase + 'persistance.pdf', bbox_inches='tight')


# eddy diameter vs. persistence (scatter)
####


# map histogram of centroid lon/lat, diameter
rownames2a = [rowname + ' lon [centroid]' for rowname in rownames]
rownames2b = [rowname + ' lat [centroid]' for rowname in rownames]
rownames2c = [rowname + ' diameter [km]' for rowname in rownames]

# convert lon/lat of eddy centroids to projected coordinates
# pts is time by eddy by [x,y,z]
pts = aed.transform_points(pc, dfeddies[rownames2a].values, dfeddies[rownames2b].values)

# remove short eddies
pts = pts[:,np.where(longeddies)[0]]

# save radii out
rads = (dfeddies[rownames2c]*1000/2).values
# remove short eddies
# this time by long eddy
rads = rads[:,np.where(longeddies)[0]]

# reshape to have time/eddy in first axis
pts = pts.reshape((pts[:,:,0].size,pts.shape[2]))[:,:2]
rads = rads.flatten()

# remove nans
inans = ~np.isnan(rads)
pts = pts[inans]
rads = rads[inans]

# make shapely Points of all eddies at any time that are long enough
ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts,rads)]

# plot aggregated eddies on map
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

[ax.add_geometries([ps[i]], crs=aed, facecolor='darkcyan', alpha=0.005, edgecolor='k', linewidth=0.5) for i in range(len(ps))];

fig.savefig('figures/swirl/hycom/map_eddies.png', bbox_inches='tight')
plt.close(fig)




## PLOT ##
# plot SSH in time with eddies identified like before, but color differently
# if used in analysis. Overlay centroid and radius. Compare with means?
# put as much on plot as possible to demonstrate.

# read in model output
loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
ds = xr.open_dataset(loc)
grid_file = 'gom03_grd_N050_new.nc'
grid = xr.open_dataset(grid_file)
LON = grid.lon_rho.values
LAT = grid.lat_rho.values
t = ds['ocean_time']
sshmax = 1.0; sshmin = -sshmax
dt = 8  # use every 8 time outputs (24 hours)
# assign ocean_time to coordinate
ds = ds.assign_coords(ocean_time=ds['ocean_time'])

base = 'figures/swirl/hycom/withmetrics'
os.makedirs(base, exist_ok=True)

for i in range(0,t.size, dt):
# i=24

    date = t[i]
    fname = '%s/%s.png' % (base,str(date.values)[:13])

    if os.path.exists(fname):
        continue

    ssh = ds['zeta'].isel(time_1=i).values
    # ssh resized to psi grid minus 1
    sshpsi = op.resize(op.resize(ssh[1:-1,1:-1], 1), 0)
    u = ds['u'].isel(time_1=i).values  # y by x
    v = ds['v'].isel(time_1=i).values

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


    # sea surface height
    mappable = ax.pcolormesh(LON, LAT, ssh, transform=pc,
                             cmap=cmo.balance, vmin=sshmin, vmax=sshmax)

    dd = 15
    ax.quiver(LON[::dd,::dd], LAT[::dd,::dd], u[::dd,::dd], v[::dd,::dd], transform=pc, alpha=0.6);

    cb = fig.colorbar(mappable, shrink=0.8, extend='both')
    cb.set_label('Sea surface height [m]')

    # add eddy identifications for this time to plot
    datepd = pd.to_datetime(date.values).strftime('%Y-%m-%d')
    polys = df[datepd]['p'][0]
    [ax.add_geometries([poly], crs=aed, facecolor='none',edgecolor='k', linewidth=3) for poly in polys];

    ## Overlay metrics ##

    # show if an eddy identified is being used in statistics
    # find which col has an eddy for this datetime to plot
    # eddies from this datetime that are long enough and aren't nan
    longeddiesnans = longeddies & (~dfeddies[datepd][rownames].isnull()).values[0]
    # dfeddies[datepd][rownames].loc[:,longeddiesnans]

    # iplot = ~dfeddies[datepd].isnull().values[0] & ~dfeddies.columns.str.contains(' ')
    # cols = dfeddies[datepd].iloc[:,iplot].columns
    cols = dfeddies[datepd][rownames].loc[:,longeddiesnans].columns
    if len(cols)>0:  # anything identified this time step
        for col in cols:
            ax.add_geometries([dfeddies[datepd][col].values[0]], crs=aed,
                              facecolor='none', edgecolor='r', linewidth=2, linestyle='--')

    ## centroid and radius ##
    # col = cols[0]
        lonnames = [col + ' lon [centroid]' for col in cols]
        latnames = [col + ' lat [centroid]' for col in cols]
        diamnames = [col + ' diameter [km]' for col in cols]
        pts = aed.transform_points(pc, dfeddies[datepd][lonnames].values, dfeddies[datepd][latnames].values)[0]
        # save radii out
        rads = (dfeddies[datepd][diamnames]*1000/2).values  # radius, in meters
        # make shapely Points of these reconstructed eddies
        ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts[:,:2],rads.T)]
        [ax.add_geometries([ps[i]], crs=aed, facecolor='none', edgecolor='b', linewidth=1) for i in range(len(ps))];

        ## overlay mean diameter size for the plotted eddy, and for all eddies together ##
        # mean per eddy
        rads = (dfstats['mean diameter [km]'][cols]*1000/2).values
        # make shapely Points of these reconstructed eddies
        ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts[:,:2],rads.T)]
        [ax.add_geometries([ps[i]], crs=aed, facecolor='none', edgecolor='b',
                           linewidth=1, linestyle='--') for i in range(len(ps))];
        # overall mean
        rad = dfstats['mean diameter [km]'].mean()*1000/2
        # make shapely Points of these reconstructed eddies
        ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for pt in pts[:,:2]]
        [ax.add_geometries([ps[i]], crs=aed, facecolor='none', edgecolor='b',
                           linewidth=1, linestyle=':') for i in range(len(ps))];

    ## eddy legend ##
    # reconstructed eddy from centroid and mean effective radius for this eddy in time
    lon = -97.75; lat = 30.75; name = 'overall mean'
    pt = aed.transform_point(lon, lat, pc)
    p = shapely.geometry.Point(pt).buffer(20000)
    ax.add_geometries([p], crs=aed, facecolor='none', edgecolor='b', linewidth=1, linestyle=':')
    ax.text(lon+0.25, lat-0.1, name, transform=pc)
    # reconstructed eddy from centroid and mean effective radius for this eddy in time
    lon = -97.75; lat = 30.25; name = 'mean eddy'
    pt = aed.transform_point(lon, lat, pc)
    p = shapely.geometry.Point(pt).buffer(20000)
    ax.add_geometries([p], crs=aed, facecolor='none', edgecolor='b', linewidth=1, linestyle='--')
    ax.text(lon+0.25, lat-0.1, name, transform=pc)
    # reconstructed eddy from centroid and effective radius from area
    lon = -97.75; lat = 29.75; name = 'eddy props'
    pt = aed.transform_point(lon, lat, pc)
    p = shapely.geometry.Point(pt).buffer(20000)
    ax.add_geometries([p], crs=aed, facecolor='none', edgecolor='b', linewidth=1)
    ax.text(lon+0.25, lat-0.1, name, transform=pc)

    ## write persistence and name at centroid ##
    for lonname, latname, col in zip(lonnames, latnames, cols):
        ax.text(dfeddies[datepd][lonname].values, dfeddies[datepd][latname].values,
                col + '-' + str(int(dfstats['persistence [days]'][col])),
                transform=pc, horizontalalignment='center', fontsize=6,
                verticalalignment='center', color='k')



    # label model output
    ax.text(0.3, 0.96, '%s model output' % 'hycom', color='0.2', transform=ax.transAxes)
    # my name
    ax.text(0.8, 0.96, 'K. Thyng', color='0.2', transform=ax.transAxes)

    ax.set_title(pd.to_datetime(date.values).strftime('%b %d %H:%M, %Y'))

    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
