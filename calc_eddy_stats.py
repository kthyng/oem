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

# df HAS TWO ENTRIES PER DAY FOR MOST OF 2009
# remove time, then remove duplicate entries
df.index = df.index.normalize()
idups = df.index.duplicated()
df = pd.DataFrame(index=df.index[~idups], data=df[~idups])

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

# df HAS TWO ENTRIES PER DAY FOR MOST OF 2009
# remove time, then remove duplicate entries
df.index = df.index.normalize()
idups = df.index.duplicated()
df = pd.DataFrame(index=df.index[idups], data=df[idups])

# initialize column to hold polygons
df['p'] = ''
df['p'] = df['p'].apply(list)

# reconstituting polygon shapes and adding them into appropriate time as list
for date in df.index:
    xs = ast.literal_eval(df.loc[date, 'x'])
    ys = ast.literal_eval(df.loc[date, 'y'])

    pstemp = [shapely.geometry.Polygon(np.vstack((x, y)).T) for (x,y) in zip(xs, ys)]

    df.loc[date,'p'] = pstemp

## End read in ##



## Stats ##

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
rownames2 = [rowname + ' speed [km/day]' for rowname in rownames if rowname + ' speed [km/day]' in dfeddies.columns]  # speed columns
rownames3 = [rowname for rowname in rownames if rowname + ' speed [km/day]' in dfeddies.columns]  # speed columns
dfstats.loc[rownames3,'mean propagation speed [km/day]'] = dfeddies[rownames2].mean().values
# nan out short time eddies
dfstats['mean propagation speed [km/day]'][shorteddies] = np.nan
# print(dfstats['mean propagation speed [km/day]'].mean(), dfstats['mean propagation speed [km/day]'].max())

dfstats.to_csv('calcs/swirl/dfstats.csv')

## end calc stats ##



## number of long term (60 days or more) cyclones in time ##



####


## print stats/show summary plots ##
sumbase = 'figures/swirl/hycom/withmetrics/summary/'

# timing during the year
nyears = 10  # number of years of output
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
ax = dfstats['mean propagation speed [km/day]'].plot(kind='hist', bins=20, grid=False, edgecolor='k')
ax.set_xlabel('Mean eddy propagation speed across lifetime [km/day]')
ax.set_ylabel('Number of cyclones longer than 7 days in bins')
# overall mean speed
dmean = dfstats['mean propagation speed [km/day]'].mean()
text0 = 'overall mean speed: %dkm/day' % dmean
ax.text(0.4, 0.95, text0, transform=ax.transAxes)
plt.savefig(sumbase + 'speed.pdf', bbox_inches='tight')

# persistence
ax = dfstats['persistence [days]'].plot(kind='hist', bins=20, grid=False, edgecolor='k')
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


# ## Plot identified eddies on map alone by time step ##
# # loop through p_i columns and calculate properties
# lastcolnum = int(dfeddies.columns[-1].split(' ')[0][1:])  # int(dfeddies.columns[-1].split('p')[1])
# # don't include last p column since should be empty
# rownames = ['p%i' % i for i in range(lastcolnum)]
#
# # bool of eddies that are "long enough" to count (more than 7 days currently)
# longeddies = dfeddies[rownames].count() > 7
#
# base = 'figures/swirl/hycom/eddiesbyregion/'
# os.makedirs(base, exist_ok=True)
#
# for datepd in dfeddies.index:
#
#     # datepd = dfeddies.index[0].strftime('%Y-%m-%d')
#
#     datepd = datepd.strftime('%Y-%m-%d')
#
#     fname = '%s/%s.png' % (base,datepd)
#
#     if os.path.exists(fname):
#         continue
#
#     fig = plt.figure(figsize=(8,6))# (9.4, 7.7))
#     ax = fig.add_subplot(111, projection=merc)
#     # ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=merc)
#     ax.set_extent(extent, pc)
#     gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
#     # the following two make the labels look like lat/lon format
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mpl.ticker.FixedLocator(np.arange(-100, -70, 2))
#     gl.ylocator = mpl.ticker.FixedLocator(np.arange(10, 31, 2))
#     gl.xlabels_top = False  # turn off labels where you don't want them
#     gl.ylabels_right = False
#     # add background land
#     ax.add_feature(land_50m, facecolor='0.8', edgecolor='0.2', linewidth=1)
#
#     # show if an eddy identified is being used in statistics
#     # find which col has an eddy for this datetime to plot
#     # eddies from this datetime that are long enough and aren't nan
#     longeddiesnans = longeddies & (~dfeddies[datepd][rownames].isnull()).values[0]
#
#     cols = dfeddies[datepd][rownames].loc[:,longeddiesnans].columns
#     if len(cols)>0:  # anything identified this time step
#         for col in cols:
#             lonname = col + ' lon [centroid]'
#             latname = col + ' lat [centroid]'
#             diamname = col + ' diameter [km]'
#             # lon and lat of eddy centroid, to determine region of eddy
#             lone, late = dfeddies[lonname].values[0], dfeddies[latname].values[0]
#             # if centroid north of 25.5 more than half of eddy lifetime
#             # and west of 90: slope eddy (NW)
#             isnorth = ((dfeddies[latname] > 25.5).sum()/dfeddies[latname].count()) > 0.5
#             iswest = ((dfeddies[lonname] < -90).sum()/dfeddies[lonname].count()) > 0.5
#             # if centroid east of 85 more than half of eddy lifetime: tortuga eddy
#             iseast = ((dfeddies[lonname] > -85).sum()/dfeddies[lonname].count()) > 0.5
#             # if centroid south of 22 more than half of eddy lifetime: campeche eddy
#             issouth = ((dfeddies[latname] < 22).sum()/dfeddies[latname].count()) > 0.5
#             if isnorth and iswest:  # slope
#                 ceddy = 'r'
#             elif iseast:  # tortuga
#                 ceddy = 'orange'
#             elif issouth:  # campeche
#                 ceddy = 'g'
#             else:  # parasitic
#                 ceddy = 'b'
#
#             # pts = aed.transform_points(pc, dfeddies[datepd][lonnames].values, dfeddies[datepd][latnames].values)[0]
#             # # save radii out
#             # rads = (dfeddies[datepd][diamnames]*1000/2).values  # radius, in meters
#             # # make shapely Points of these reconstructed eddies
#             # ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts[:,:2],rads.T)]
#             # [ax.add_geometries([ps[i]], crs=aed, facecolor='none', edgecolor='b', linewidth=1) for i in range(len(ps))];
#
#             ax.add_geometries([dfeddies[datepd][col].values[0]], crs=aed, alpha=0.5,
#                               facecolor=ceddy, edgecolor='none', linewidth=2, linestyle='-')
#
#     # label model output
#     ax.text(0.3, 0.96, '%s model output' % 'hycom', color='0.2', transform=ax.transAxes)
#     # my name
#     ax.text(0.8, 0.96, 'K. Thyng', color='0.2', transform=ax.transAxes)
#
#     ax.set_title(datepd)
#
#     fig.savefig(fname, bbox_inches='tight')
#     plt.close(fig)

####

##### cyclone stats by region #####

## Plot in time of cyclones detected by region and size ##

# rownames of only long eddies
rownameslong = ['p%i' % i for i in range(lastcolnum) if longeddies[i]]

lonnames = [col + ' lon [centroid]' for col in rownameslong]
latnames = [col + ' lat [centroid]' for col in rownameslong]
diamnames = [col + ' diameter [km]' for col in rownameslong]
# lon and lat of eddy centroid, to determine region of eddy
lone, late = dfeddies[lonnames], dfeddies[latnames]
# lone, late = dfeddies[datepd][lonnames].values[0], dfeddies[datepd][latnames].values[0]
# if centroid north of 25.5 more than half of eddy lifetime
# and west of 90: slope eddy (NW)
isnorth = ((late > 25.5).sum(axis=0)/late.count(axis=0)) > 0.5
iswest = ((lone < -90).sum(axis=0)/lone.count(axis=0)) > 0.5
# if centroid east of 85 more than half of eddy lifetime: tortuga eddy
iseast = ((lone > -85).sum()/lone.count()) > 0.5
# if centroid south of 22 more than half of eddy lifetime: campeche eddy
issouth = ((late < 22).sum()/late.count()) > 0.5

islope = np.where(isnorth.values * iswest.values)[0]
itort = np.where(iseast.values)[0]
icamp = np.where(issouth.values)[0]
iother = np.where(~(isnorth.values * iswest.values) * ~iseast.values
                  * ~issouth.values)[0]

# dataframe in time to save timing
dfepresent = ~dfeddies[rownameslong].isnull()  # when eddies are present




## number of eddies in time present by region ##
dfepresentnan = ~dfeddies[rownameslong].isnull()  # when eddies are present
dfepresentnan[~dfepresentnan] = np.nan


## histograms ##

# lifetime of eddies by region in histogram
fig, axes = plt.subplots(4,1, sharex=True, sharey=True, figsize=(10,6))
[ax.set_yticks(np.arange(0,175,25)) for ax in axes]
[ax.set_xticks(np.arange(0,150,15)) for ax in axes]
dfeddies[rownameslong].iloc[:,islope].count().plot(kind='hist', color='r', alpha=0.3, bins=20, range=(0,120), ax=axes[0], grid=True, logy=True, edgecolor='k')
dfeddies[rownameslong].iloc[:,itort].count().plot(kind='hist', color='orange', alpha=0.3, bins=20, range=(0,120), ax=axes[1], grid=True, logy=True, edgecolor='k')
dfeddies[rownameslong].iloc[:,icamp].count().plot(kind='hist', color='g', alpha=0.3, bins=20, range=(0,120), ax=axes[2], grid=True, logy=True, edgecolor='k')
dfeddies[rownameslong].iloc[:,iother].count().plot(kind='hist', color='b', alpha=0.3, bins=20, range=(0,120), ax=axes[3], grid=True, logy=True, edgecolor='k')
axes[0].set_title('number of eddies over 10 year period')
axes[-1].set_xlabel('Cyclone lifetime [days]')
labels = ['Slope cyclones', 'Tortuga', 'Campeche', 'Parasitic']
colors = ['r', 'orange', 'g', 'b']
[ax.text(90, 50, label, color=color, fontsize=12) for ax, color, label in zip(axes, colors, labels)]
fname = 'figures/swirl/hycom/eddiesbyregion/stats_hist_lifetime'
fig.savefig('%s.png' % fname, bbox_inches='tight')

# mean cyclone diamater by region in histogram
rownames2 = [rowname + ' diameter [km]' for rowname in rownameslong]  # diameter columns
fig, axes = plt.subplots(4,1, sharex=True, sharey=True, figsize=(10,6))
[ax.set_yticks(np.arange(0,300,25)) for ax in axes]
[ax.set_xticks(np.arange(0,300,15)) for ax in axes]
dfeddies[rownames2].iloc[:,islope].mean().plot(kind='hist', color='r', alpha=0.3, bins=20, range=(75,225), ax=axes[0], grid=True, logy=True, edgecolor='k')
dfeddies[rownames2].iloc[:,itort].mean().plot(kind='hist', color='orange', alpha=0.3, bins=20, range=(75,225), ax=axes[1], grid=True, logy=True, edgecolor='k')
dfeddies[rownames2].iloc[:,icamp].mean().plot(kind='hist', color='g', alpha=0.3, bins=20, range=(75,225), ax=axes[2], grid=True, logy=True, edgecolor='k')
dfeddies[rownames2].iloc[:,iother].mean().plot(kind='hist', color='b', alpha=0.3, bins=20, range=(75,225), ax=axes[3], grid=True, logy=True, edgecolor='k')
axes[0].set_title('number of eddies over 10 year period')
axes[-1].set_xlabel('Cyclone mean diameter [km]')
labels = ['Slope cyclones', 'Tortuga', 'Campeche', 'Parasitic']
colors = ['r', 'orange', 'g', 'b']
[ax.text(195, 30, label, color=color, fontsize=12) for ax, color, label in zip(axes, colors, labels)]
fname = 'figures/swirl/hycom/eddiesbyregion/stats_hist_diam'
fig.savefig('%s.png' % fname, bbox_inches='tight')

# # mean propagation speed by region in histogram
# # not all cyclones are being included, so not doing this for now
# fig, axes = plt.subplots(4,1, sharex=True, sharey=True, figsize=(10,6))
# [ax.set_yticks(np.arange(0,300,10)) for ax in axes]
# [ax.set_xticks(np.arange(0,300,5)) for ax in axes]
# key = 'mean propagation speed [km/day]'
# dfstats[key].iloc[islope].plot(kind='hist', color='r', alpha=0.3, bins=20, range=(0,35), ax=axes[0], grid=True, logy=True, edgecolor='k')
# dfstats[key].iloc[itort].plot(kind='hist', color='orange', alpha=0.3, bins=20, range=(0,35), ax=axes[1], grid=True, logy=True, edgecolor='k')
# dfstats[key].iloc[icamp].plot(kind='hist', color='g', alpha=0.3, bins=20, range=(0,35), ax=axes[2], grid=True, logy=True, edgecolor='k')
# dfstats[key].iloc[iother].plot(kind='hist', color='b', alpha=0.3, bins=20, range=(0,35), ax=axes[3], grid=True, logy=True, edgecolor='k')
# axes[0].set_title('number of eddies over 11 year period')
# axes[-1].set_xlabel('Cyclone mean propagation speed [km/day]')
# labels = ['Slope cyclones', 'Tortuga', 'Campeche', 'Parasitic']
# colors = ['r', 'orange', 'g', 'b']
# [ax.text(25, 5, label, color=color, fontsize=12) for ax, color, label in zip(axes, colors, labels)]
# fname = 'figures/swirl/hycom/eddiesbyregion/stats_hist_speed'
# fig.savefig('%s.png' % fname, bbox_inches='tight')

# cyclones per month of the year across years, by region

# calculate number of eddies per month
# number of days eddy was present for each month of year
numdays = (~dfeddies[rownameslong].isnull()).groupby(dfeddies.index.month).sum()
# just counting presence of eddy
freq_islope = (numdays.iloc[:,islope] > 0).sum(axis='columns')
freq_itort = (numdays.iloc[:,itort] > 0).sum(axis='columns')
freq_icamp = (numdays.iloc[:,icamp] > 0).sum(axis='columns')
freq_iother = (numdays.iloc[:,iother] > 0).sum(axis='columns')
freqstd_islope = [s.loc[i,:].iloc[:,islope].count(axis=1).std() for i in range(1,13)]
freqstd_itort = [s.loc[i,:].iloc[:,itort].count(axis=1).std() for i in range(1,13)]
freqstd_icamp = [s.loc[i,:].iloc[:,icamp].count(axis=1).std() for i in range(1,13)]
freqstd_iother = [s.loc[i,:].iloc[:,iother].count(axis=1).std() for i in range(1,13)]

nyears = 11  # number of years of output
ax = (freq_islope/nyears).plot(figsize=(10,4), grid=True, color='r', lw=2)
(freq_itort/nyears).plot(figsize=(10,4), grid=True, color='orange', ax=ax, lw=2)
(freq_icamp/nyears).plot(figsize=(10,4), grid=True, color='g', ax=ax, lw=2)
(freq_iother/nyears).plot(figsize=(10,4), grid=True, color='b', ax=ax, lw=2)
ax.set_xlabel('Months of the year')
ax.set_ylabel('Number of eddies per year\nstandard deviation across 11 years')
ax.fill_between(freq.index, (freq_islope/nyears)-freqstd_islope, (freq_islope/nyears)+freqstd_islope, alpha=0.2, color='r')
ax.fill_between(freq.index, (freq_itort/nyears)-freqstd_itort, (freq_itort/nyears)+freqstd_itort, alpha=0.2, color='orange')
ax.fill_between(freq.index, (freq_icamp/nyears)-freqstd_icamp, (freq_icamp/nyears)+freqstd_icamp, alpha=0.2, color='g')
ax.fill_between(freq.index, (freq_iother/nyears)-freqstd_iother, (freq_iother/nyears)+freqstd_iother, alpha=0.2, color='b')
ax.set_xlim(1,12)
ax.set_xticks(np.arange(1,13));
ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
fname = 'figures/swirl/hycom/eddiesbyregion/stats_eddy_freq'
plt.savefig('%s.png' % fname, bbox_inches='tight')

#########

### PLOT ALL EDDYS but colored by region ###
# map histogram of centroid lon/lat, diameter
rislopelon = [rownameslong[i] + ' lon [centroid]' for i in islope]
rislopelat = [rownameslong[i] + ' lat [centroid]' for i in islope]
rislopedia = [rownameslong[i] + ' diameter [km]' for i in islope]

ritortlon = [rownameslong[i] + ' lon [centroid]' for i in itort]
ritortlat = [rownameslong[i] + ' lat [centroid]' for i in itort]
ritortdia = [rownameslong[i] + ' diameter [km]' for i in itort]

ricamplon = [rownameslong[i] + ' lon [centroid]' for i in icamp]
ricamplat = [rownameslong[i] + ' lat [centroid]' for i in icamp]
ricampdia = [rownameslong[i] + ' diameter [km]' for i in icamp]

riotherlon = [rownameslong[i] + ' lon [centroid]' for i in iother]
riotherlat = [rownameslong[i] + ' lat [centroid]' for i in iother]
riotherdia = [rownameslong[i] + ' diameter [km]' for i in iother]


# convert lon/lat of eddy centroids to projected coordinates
# pts is time by eddy by [x,y,z]
pts_islope = aed.transform_points(pc, dfeddies[rislopelon].values, dfeddies[rislopelat].values)
pts_itort = aed.transform_points(pc, dfeddies[ritortlon].values, dfeddies[ritortlat].values)
pts_icamp = aed.transform_points(pc, dfeddies[ricamplon].values, dfeddies[ricamplat].values)
pts_iother = aed.transform_points(pc, dfeddies[riotherlon].values, dfeddies[riotherlat].values)
# pts = aed.transform_points(pc, dfeddies[rownames2a].values, dfeddies[rownames2b].values)

# # remove short eddies
# pts = pts[:,np.where(longeddies)[0]]

# save radii out
rads_islope = (dfeddies[rislopedia]*1000/2).values
rads_itort = (dfeddies[ritortdia]*1000/2).values
rads_icamp = (dfeddies[ricampdia]*1000/2).values
rads_iother = (dfeddies[riotherdia]*1000/2).values
# rads = (dfeddies[rownames2c]*1000/2).values
# # remove short eddies
# # this time by long eddy
# rads = rads[:,np.where(longeddies)[0]]

# reshape to have time/eddy in first axis
pts_islope = pts_islope.reshape((pts_islope[:,:,0].size,pts_islope.shape[2]))[:,:2]
pts_itort = pts_itort.reshape((pts_itort[:,:,0].size,pts_itort.shape[2]))[:,:2]
pts_icamp = pts_icamp.reshape((pts_icamp[:,:,0].size,pts_icamp.shape[2]))[:,:2]
pts_iother = pts_iother.reshape((pts_iother[:,:,0].size,pts_iother.shape[2]))[:,:2]
rads_islope = rads_islope.flatten()
rads_itort = rads_itort.flatten()
rads_icamp = rads_icamp.flatten()
rads_iother = rads_iother.flatten()
# pts = pts.reshape((pts[:,:,0].size,pts.shape[2]))[:,:2]
# rads = rads.flatten()

# remove nans
inans = ~np.isnan(rads_islope)
pts_islope = pts_islope[inans]
rads_islope = rads_islope[inans]

inans = ~np.isnan(rads_itort)
pts_itort = pts_itort[inans]
rads_itort = rads_itort[inans]

inans = ~np.isnan(rads_icamp)
pts_icamp = pts_icamp[inans]
rads_icamp = rads_icamp[inans]

inans = ~np.isnan(rads_iother)
pts_iother = pts_iother[inans]
rads_iother = rads_iother[inans]
# inans = ~np.isnan(rads)
# pts = pts[inans]
# rads = rads[inans]

# make shapely Points of all eddies at any time that are long enough
ps_islope = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_islope,rads_islope)]
ps_itort = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_itort,rads_itort)]
ps_icamp = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_icamp,rads_icamp)]
ps_iother = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_iother,rads_iother)]
# ps = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts,rads)]

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

[ax.add_geometries([ps_islope[i]], crs=aed, facecolor='r', alpha=0.005) for i in range(len(ps_islope))];
[ax.add_geometries([ps_itort[i]], crs=aed, facecolor='orange', alpha=0.005) for i in range(len(ps_itort))];
[ax.add_geometries([ps_icamp[i]], crs=aed, facecolor='g', alpha=0.005) for i in range(len(ps_icamp))];
[ax.add_geometries([ps_iother[i]], crs=aed, facecolor='b', alpha=0.005) for i in range(len(ps_iother))];
# [ax.add_geometries([ps[i]], crs=aed, facecolor='r', alpha=0.5) for i in range(len(ps[:5]))];

fig.savefig('figures/swirl/hycom/eddiesbyregion/map_eddies.png', bbox_inches='tight')
plt.close(fig)

######


## how many of which eddies that last how long at once ##
# remove timing from days
cols_islope = [rownameslong[i] for i in islope]
cols_itort = [rownameslong[i] for i in itort]
cols_icamp = [rownameslong[i] for i in icamp]
cols_iother = [rownameslong[i] for i in iother]
years = list(np.arange(2002,2013))
years.pop(-2)
# fig, axes = plt.subplots(1,1, figsize=(8.5,1))
fig, axes = plt.subplots(10,1, figsize=(8.5,11))
for i, year in enumerate(years):
# year = 2009
    year = str(year); ax=axes[i]
    dy = 0; color = 'r'; Y = 0
    for col in cols_islope:

        # skip loop if cyclone is not present in this year
        if dfepresent[year].iloc[:,islope][col].sum() == 0:
            continue
        dy += 0.05
        # if first time of col's cyclone is when there are no other cyclones of this
        # type, put back to initial level
        irow = np.where(dfepresent[year].iloc[:,islope][col])[0][0]
        if dfepresentnan[year].iloc[:,islope].count(axis=1).iloc[irow-5] == 0 \
            and dfepresentnan[year].iloc[:,islope].count(axis=1).iloc[irow] == 1:
            dy = 0
        (dfepresentnan[year].iloc[:,islope][col] - Y - dy).plot(color=color, lw=1.5, alpha=0.5, ax=ax)

    dy = 0; color = 'orange'; Y = 0.3
    for col in cols_itort:

        # skip loop if cyclone is not present in this year
        if dfepresent[year].iloc[:,itort][col].sum() == 0:
            continue
        dy += 0.05
        # if first time of col's cyclone is when there are no other cyclones of this
        # type, put back to initial level
        irow = np.where(dfepresent[year].iloc[:,itort][col])[0][0]
        if dfepresentnan[year].iloc[:,itort].count(axis=1).iloc[irow-5] == 0 \
            and dfepresentnan[year].iloc[:,itort].count(axis=1).iloc[irow] == 1:
            dy = 0
        (dfepresentnan[year].iloc[:,itort][col] - Y - dy).plot(color=color, lw=1.5, alpha=0.5, ax=ax)

    dy = 0; color = 'g'; Y = 0.6
    for col in cols_icamp:

        # skip loop if cyclone is not present in this year
        if dfepresent[year].iloc[:,icamp][col].sum() == 0:
            continue
        dy += 0.05
        # if first time of col's cyclone is when there are no other cyclones of this
        # type, put back to initial level
        irow = np.where(dfepresent[year].iloc[:,icamp][col])[0][0]
        if dfepresentnan[year].iloc[:,icamp].count(axis=1).iloc[irow-5] == 0 \
            and dfepresentnan[year].iloc[:,icamp].count(axis=1).iloc[irow] == 1:
            dy = 0
        (dfepresentnan[year].iloc[:,icamp][col] - Y - dy).plot(color=color, lw=1.5, alpha=0.5, ax=ax)

    dy = 0; color = 'b'; Y = 0.9
    for col in cols_iother:

        # skip loop if cyclone is not present in this year
        if dfepresent[year].iloc[:,iother][col].sum() == 0:
            continue
        dy += 0.05
        # if first time of col's cyclone is when there are no other cyclones of this
        # type, put back to initial level
        irow = np.where(dfepresent[year].iloc[:,iother][col])[0][0]
        if dfepresentnan[year].iloc[:,iother].count(axis=1).iloc[irow-5] == 0 \
            and dfepresentnan[year].iloc[:,iother].count(axis=1).iloc[irow] == 1:
            dy = 0
        (dfepresentnan[year].iloc[:,iother][col] - Y - dy).plot(color=color, lw=1.5, alpha=0.5, ax=ax)

    ax.set_yticklabels(''); ax.set_yticks([])
    # minor = mpl.dates.DayLocator(bymonthday=14)
    # ax.xaxis.set_minor_locator(minor)
    # weekly major ticks
    major = mpl.dates.DayLocator(bymonthday=1)
    ax.xaxis.set_major_locator(major)
    ax.set_xlim(year + '-1-1', year + '-12-31')
    ax.grid(which='major', lw=1.5, color='k', alpha=0.25)
    # ax.grid(which='minor', lw=1.0, color='k', alpha=0.5)
    ax.set_ylabel(year, fontsize=14)

    if year != '2012':
        ax.set_xticklabels('')#; ax.set_xticks([])
    else:
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b'))
fig.subplots_adjust(bottom=0.03, top=0.99, hspace=0.05, left=0.1, right=0.92)
fig.savefig('figures/swirl/hycom/eddiesbyregion/eddiesintime.png', bbox_inches='tight')
####


### Eddy size vs. lifetime by region ###

# rownames of only long eddies
rownameslong = ['p%i' % i for i in range(lastcolnum) if longeddies[i]]

lonnames = [col + ' lon [centroid]' for col in rownameslong]
latnames = [col + ' lat [centroid]' for col in rownameslong]
diamnames = [col + ' diameter [km]' for col in rownameslong]
# lon and lat of eddy centroid, to determine region of eddy
lone, late = dfeddies[lonnames], dfeddies[latnames]
# lone, late = dfeddies[datepd][lonnames].values[0], dfeddies[datepd][latnames].values[0]
# if centroid north of 25.5 more than half of eddy lifetime
# and west of 90: slope eddy (NW)
isnorth = ((late > 25.5).sum(axis=0)/late.count(axis=0)) > 0.5
iswest = ((lone < -90).sum(axis=0)/lone.count(axis=0)) > 0.5
# if centroid east of 85 more than half of eddy lifetime: tortuga eddy
iseast = ((lone > -85).sum()/lone.count()) > 0.5
# if centroid south of 22 more than half of eddy lifetime: campeche eddy
issouth = ((late < 22).sum()/late.count()) > 0.5

islope = np.where(isnorth.values * iswest.values)[0]
itort = np.where(iseast.values)[0]
icamp = np.where(issouth.values)[0]
iother = np.where(~(isnorth.values * iswest.values) * ~iseast.values
                  * ~issouth.values)[0]

# correlation for slope eddies
# dfeddies[rownameslong].iloc[:,islope].count().plot(kind='hist', color='r', alpha=0.3, bins=20, range=(0,120), ax=axes[0], grid=True, logy=True, edgecolor='k')
ax = dfstats.loc[rownameslong,:].iloc[islope,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='r', alpha=0.5)
dfstats.loc[rownameslong,:].iloc[itort,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='orange', ax=ax, alpha=0.5)
dfstats.loc[rownameslong,:].iloc[icamp,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='g', ax=ax, alpha=0.5)
dfstats.loc[rownameslong,:].iloc[iother,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='b', ax=ax, alpha=0.5)

fig, axes = plt.subplots(2,2,figsize=(8,7), sharex=True, sharey=True)
dfstats.loc[rownameslong,:].iloc[islope,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='r', alpha=0.5, ax=axes[0][0])
dfstats.loc[rownameslong,:].iloc[itort,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='orange', ax=axes[0][1], alpha=0.5)
dfstats.loc[rownameslong,:].iloc[icamp,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='g', ax=axes[1][0], alpha=0.5)
dfstats.loc[rownameslong,:].iloc[iother,:].plot(kind='scatter', x='persistence [days]', y='mean diameter [km]', color='b', ax=axes[1][1], alpha=0.5)
fig.suptitle('Identified cyclones over 10 years', y=0.91)
names = ['Slope eddies', 'Tortuga', 'Campeche', 'Other']
colors = ['r', 'orange', 'g', 'b']
for ax, name, color in zip(axes.flatten(),names, colors):
    ax.text(0.6, 0.9, name, color=color, transform=ax.transAxes, fontsize=12)
fig.savefig('figures/swirl/hycom/eddiesbyregion/diamvlife.png', bbox_inches='tight')
######


## plot aggregated eddies on map but ONLY LONG ##
lastcolnum = int(dfeddies.columns[-1].split(' ')[0][1:])  # int(dfeddies.columns[-1].split('p')[1])
# don't include last p column since should be empty
rownames = ['p%i' % i for i in range(lastcolnum)]
longeddies30 = dfeddies[rownames].count() > 30
rownameslong30 = ['p%i' % i for i in range(lastcolnum) if longeddies30[i]]

lonnames = [col + ' lon [centroid]' for col in rownameslong30]
latnames = [col + ' lat [centroid]' for col in rownameslong30]
diamnames = [col + ' diameter [km]' for col in rownameslong30]
# lon and lat of eddy centroid, to determine region of eddy
lone, late = dfeddies[lonnames], dfeddies[latnames]
# lone, late = dfeddies[datepd][lonnames].values[0], dfeddies[datepd][latnames].values[0]
# if centroid north of 25.5 more than half of eddy lifetime
# and west of 90: slope eddy (NW)
isnorth = ((late > 25.5).sum(axis=0)/late.count(axis=0)) > 0.5
iswest = ((lone < -90).sum(axis=0)/lone.count(axis=0)) > 0.5
# if centroid east of 85 more than half of eddy lifetime: tortuga eddy
iseast = ((lone > -85).sum()/lone.count()) > 0.5
# if centroid south of 22 more than half of eddy lifetime: campeche eddy
issouth = ((late < 22).sum()/late.count()) > 0.5

islope30 = np.where(isnorth.values * iswest.values)[0]
itort30 = np.where(iseast.values)[0]
icamp30 = np.where(issouth.values)[0]
iother30 = np.where(~(isnorth.values * iswest.values) * ~iseast.values
                  * ~issouth.values)[0]

# map histogram of centroid lon/lat, diameter
rislopelon30 = [rownameslong30[i] + ' lon [centroid]' for i in islope30]
rislopelat30 = [rownameslong30[i] + ' lat [centroid]' for i in islope30]
rislopedia30 = [rownameslong30[i] + ' diameter [km]' for i in islope30]

ritortlon30 = [rownameslong30[i] + ' lon [centroid]' for i in itort30]
ritortlat30 = [rownameslong30[i] + ' lat [centroid]' for i in itort30]
ritortdia30 = [rownameslong30[i] + ' diameter [km]' for i in itort30]

ricamplon30 = [rownameslong30[i] + ' lon [centroid]' for i in icamp30]
ricamplat30 = [rownameslong30[i] + ' lat [centroid]' for i in icamp30]
ricampdia30 = [rownameslong30[i] + ' diameter [km]' for i in icamp30]

riotherlon30 = [rownameslong30[i] + ' lon [centroid]' for i in iother30]
riotherlat30 = [rownameslong30[i] + ' lat [centroid]' for i in iother30]
riotherdia30 = [rownameslong30[i] + ' diameter [km]' for i in iother30]


# convert lon/lat of eddy centroids to projected coordinates
# pts is time by eddy by [x,y,z]
pts_islope30 = aed.transform_points(pc, dfeddies[rislopelon30].values, dfeddies[rislopelat30].values)
pts_itort30 = aed.transform_points(pc, dfeddies[ritortlon30].values, dfeddies[ritortlat30].values)
pts_icamp30 = aed.transform_points(pc, dfeddies[ricamplon30].values, dfeddies[ricamplat30].values)
pts_iother30 = aed.transform_points(pc, dfeddies[riotherlon30].values, dfeddies[riotherlat30].values)

# save radii out
rads_islope30 = (dfeddies[rislopedia30]*1000/2).values
rads_itort30 = (dfeddies[ritortdia30]*1000/2).values
rads_icamp30 = (dfeddies[ricampdia30]*1000/2).values
rads_iother30 = (dfeddies[riotherdia30]*1000/2).values

# reshape to have time/eddy in first axis
pts_islope30 = pts_islope30.reshape((pts_islope30[:,:,0].size,pts_islope30.shape[2]))[:,:2]
pts_itort30 = pts_itort30.reshape((pts_itort30[:,:,0].size,pts_itort30.shape[2]))[:,:2]
pts_icamp30 = pts_icamp30.reshape((pts_icamp30[:,:,0].size,pts_icamp30.shape[2]))[:,:2]
pts_iother30 = pts_iother30.reshape((pts_iother30[:,:,0].size,pts_iother30.shape[2]))[:,:2]
rads_islope30 = rads_islope30.flatten()
rads_itort30 = rads_itort30.flatten()
rads_icamp30 = rads_icamp30.flatten()
rads_iother30 = rads_iother30.flatten()

# remove nans
inans = ~np.isnan(rads_islope30)
pts_islope30 = pts_islope30[inans]
rads_islope30 = rads_islope30[inans]

inans = ~np.isnan(rads_itort30)
pts_itort30 = pts_itort30[inans]
rads_itort30 = rads_itort30[inans]

inans = ~np.isnan(rads_icamp30)
pts_icamp30 = pts_icamp30[inans]
rads_icamp30 = rads_icamp30[inans]

inans = ~np.isnan(rads_iother30)
pts_iother30 = pts_iother30[inans]
rads_iother30 = rads_iother30[inans]

ps_islope30 = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_islope30,rads_islope30)]
ps_itort30 = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_itort30,rads_itort30)]
ps_icamp30 = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_icamp30,rads_icamp30)]
ps_iother30 = [shapely.geometry.Point(pt[0], pt[1]).buffer(rad) for (pt, rad) in zip(pts_iother30,rads_iother30)]

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

[ax.add_geometries([ps_islope30[i]], crs=aed, facecolor='r', alpha=0.005) for i in range(len(ps_islope30))];
[ax.add_geometries([ps_itort30[i]], crs=aed, facecolor='orange', alpha=0.005) for i in range(len(ps_itort30))];
[ax.add_geometries([ps_icamp30[i]], crs=aed, facecolor='g', alpha=0.005) for i in range(len(ps_icamp30))];
[ax.add_geometries([ps_iother30[i]], crs=aed, facecolor='b', alpha=0.005) for i in range(len(ps_iother30))];
# [ax.add_geometries([ps[i]], crs=aed, facecolor='r', alpha=0.5) for i in range(len(ps[:5]))];
ax.set_title('Cyclones over 30 days')
fig.savefig('figures/swirl/hycom/eddiesbyregion/map_eddies_30.png', bbox_inches='tight')
plt.close(fig)

####

## Make excel file of all centroid lon/lat ##

# drop all rows that aren't centroid coords
lonnames = [col + ' lon [centroid]' for col in rownameslong]
latnames = [col + ' lat [centroid]' for col in rownameslong]
dfeddies[[lonnames+latnames]]
dfeddies.loc[:,lonnames+latnames].to_csv('centroids.csv')
# rlons = [rownameslong[i] + ' lon [centroid]' for i in range(len(rownameslong))]
# rlats = [rownameslong[i] + ' lat [centroid]' for i in range(len(rownameslong))]
#
# # inds = np.where(~dfeddies[rlons].isnull())
#
# for rlon, rlat in zip(rlons, rlats):
#     inds = np.where(~dfeddies[rlon].isnull())[0]
#     dfeddies[[rlon,rlat]].iloc[inds]
####

# isnorth = ((dfeddies[latnames] > 25.5).sum()/dfeddies[latnames].count()) > 0.5
# iswest = ((dfeddies[lonnames] < -90).sum()/dfeddies[lonnames].count()) > 0.5
# # if centroid east of 85 more than half of eddy lifetime: tortuga eddy
# iseast = ((dfeddies[lonnames] > -85).sum()/dfeddies[lonnames].count()) > 0.5
# # if centroid south of 22 more than half of eddy lifetime: campeche eddy
# issouth = ((dfeddies[latnames] < 22).sum()/dfeddies[latnames].count()) > 0.5
# islope = np.where(isnorth.values * iswest.values)[0]
# itort = np.where(iseast.values)[0]
# icamp = np.where(issouth.values)[0]
# # other
# otherinds = np.concatenate((islope,itort,icamp))
# # iother = list(np.arange(0,len(neweddies)))
# # for ind in otherinds:
# #     iother.pop(ind)
# # iother = np.asarray(iother)
# othereddies = list(set(cols) & set(neweddies))
# iother = []
# for othereddy in othereddies:
#     iother.append(cols.get_loc(othereddy))
# iother = np.asarray(iother)

# Can pre-process and create dataframe with color and level info prescribed,
# then make the plot itself from the dataframe

# fig, ax = plt.subplots(1,1, figsize=(12,6))
#
# for i, datepd in enumerate(dfeddies.index[:30]):
#
#     datepd = datepd.strftime('%Y-%m-%d')
#
#     # previous date, to get previous eddies
#     # if i>0:
#     datepdold = dfeddies.index[i-1].strftime('%Y-%m-%d')
#
# # datepd = dfeddies.index[0].strftime('%Y-%m-%d')
#
#     # show if an eddy identified is being used in statistics
#     # find which col has an eddy for this datetime to plot
#     # eddies from this datetime that are long enough and aren't nan
#     longeddiesnans = longeddies & (~dfeddies[datepd][rownames].isnull()).values[0]
#     longeddiesnansold = longeddies & (~dfeddies[datepdold][rownames].isnull()).values[0]
#
#     cols = dfeddies[datepd][rownames].loc[:,longeddiesnans].columns
#     colsold = dfeddies[datepdold][rownames].loc[:,longeddiesnansold].columns
#
#     if len(cols)>0:  # anything identified this time step
#
#         ceddy = ['b'] * len(cols)  # blue is default color, parasitic eddy
#         ylevel = [1] * len(cols)
#
#         # first see if any of the eddies are the same as last time step
#         repeateddies = set(cols.values) & set(colsold.values)
#         if len(repeateddies) > 0:
#             for colh in repeateddies:
#                 iedold = colsold.get_loc(colh)  # index of this eddy from previous time step
#                 ied = cols.get_loc(colh)  # index of this eddy from this time step
#                 ylevel[ied] = ylevelold[iedold]  # use same ylevel as for eddy from previous time step since the same
#                 ceddy[ied] = ceddyold[iedold]  # use same ylevel as for eddy from previous time step since the same
#
#         # now deal with new eddies from this time step
#         neweddies = list(set(cols.values) - repeateddies)
#         if len(neweddies) > 0:
#             # for colh in neweddies:
#
#             lonnames = [col + ' lon [centroid]' for col in neweddies]
#             latnames = [col + ' lat [centroid]' for col in neweddies]
#             diamnames = [col + ' diameter [km]' for col in neweddies]
#             # lon and lat of eddy centroid, to determine region of eddy
#             lone, late = dfeddies[datepd][lonnames].values[0], dfeddies[datepd][latnames].values[0]
#             # if centroid north of 25.5 more than half of eddy lifetime
#             # and west of 90: slope eddy (NW)
#             isnorth = ((dfeddies[latnames] > 25.5).sum()/dfeddies[latnames].count()) > 0.5
#             iswest = ((dfeddies[lonnames] < -90).sum()/dfeddies[lonnames].count()) > 0.5
#             # if centroid east of 85 more than half of eddy lifetime: tortuga eddy
#             iseast = ((dfeddies[lonnames] > -85).sum()/dfeddies[lonnames].count()) > 0.5
#             # if centroid south of 22 more than half of eddy lifetime: campeche eddy
#             issouth = ((dfeddies[latnames] < 22).sum()/dfeddies[latnames].count()) > 0.5
#             islope = np.where(isnorth.values * iswest.values)[0]
#             itort = np.where(iseast.values)[0]
#             icamp = np.where(issouth.values)[0]
#             # other
#             otherinds = np.concatenate((islope,itort,icamp))
#             # iother = list(np.arange(0,len(neweddies)))
#             # for ind in otherinds:
#             #     iother.pop(ind)
#             # iother = np.asarray(iother)
#             othereddies = list(set(cols) & set(neweddies))
#             iother = []
#             for othereddy in othereddies:
#                 iother.append(cols.get_loc(othereddy))
#             iother = np.asarray(iother)
#             # if otherinds.size > 0:
#             #     iother = np.where(list(set(neweddies) - set(neweddies[islope]) - set(neweddies[itort]) - set(neweddies[icamp])) == neweddies)[0]
#             # else:
#             if islope.size > 0:
#                 ceddy[islope[0]] = 'r'
#                 for ind in islope:  # loop over these eddies
#                     ylev = 4#; j = 0
#                     while ylev in ylevel:
#                         ylev -= 0.25
#                     ylevel[ind] = ylev
#                 # inds = ylevel == 4
#                 # if len(inds)>0:  # if already have this value
#                 #     for j, ind in enumerate(islope[0]):  # loop over indices
#                 #         ylevel[ind] = 4 - j*0.25
#             if itort.size > 0:
#                 ceddy[itort[0]] = 'orange'
#                 for ind in itort:  # loop over these eddies
#                     ylev = 3#; j = 0
#                     while ylev in ylevel:
#                         ylev -= 0.25
#                     ylevel[ind] = ylev
#                 # inds = ylevel == 3
#                 # if len(inds)>0:  # if already have this value
#                 #     for j, ind in enumerate(itort[0]):  # loop over indices
#                 #         ylevel[ind] = 3 - j*0.25
#             if icamp.size > 0:
#                 ceddy[icamp[0]] = 'g'
#                 for ind in icamp:  # loop over these eddies
#                     ylev = 2#; j = 0
#                     while ylev in ylevel:
#                         ylev -= 0.25
#                     ylevel[ind] = ylev
#                 # inds = ylevel == 2
#                 # if len(inds)>0:  # if already have this value
#                 #     for j, ind in enumerate(icamp[0]):  # loop over indices
#                 #         ylevel[ind] = 2 - j*0.25
#             if iother.size > 0:
#                 ceddy[iother[0]] = 'b'
#                 for ind in iother:  # loop over these eddies
#                     ylev = 1#; j = 0
#                     while ylev in ylevel:
#                         ylev -= 0.25
#                     ylevel[ind] = ylev
#
#                 # inds = ylevel == 1
#                 # if len(inds)>0:  # if already have this value
#                 #     j = 0
#                 #     for ind in iother[0]:  # loop over indices
#                 #         ylevel[ind] = 1 - j*0.25
#
#             # check if any of the new values are already used
#
#
#
#         # save ylevels for next loop in case same eddy present
#         ylevelold = ylevel.copy()
#         ceddyold = ceddy.copy()
#
#         diamnames = [col + ' diameter [km]' for col in cols]
#         rads = (dfeddies[datepd][diamnames]*1000/2).values
#         ax.scatter([datepd] * len(ylevel), ylevel, c=ceddy, s=2*np.pi*rads/1000/5, alpha=0.25)

####





## PLOT ##
# plot SSH in time with eddies identified like before, but color differently
# if used in analysis. Overlay centroid and radius. Compare with means?
# put as much on plot as possible to demonstrate.

withmetrics = False
longeddies = True  # color long-term eddies distinct from short

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

base = 'figures/swirl/hycom'
if withmetrics:
    base += '/withmetrics/'
if longeddies:
    base += '/longeddies/'
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
    if withmetrics:
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
