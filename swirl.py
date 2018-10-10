'''Calculate swirl strength.'''

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

whichmodel = 'hycom'  # 'mercator'  # 'hycom'

plotdrifters = False
spacefilter = False
timeaverage = False

# year = 2002

# make projection
inputs = {'proj': 'aeqd', 'ellps': 'clrk66', 'datum': 'NAD27',
                  'lat_0': 26, 'lon_0': -90,
                  'x_0': 0, 'y_0': 0}
proj = Proj(**inputs)

# read in model output
if whichmodel == 'mercator':
    loc = 'drifter_animation/output/*.nc'
    ds = xr.open_mfdataset(loc)
    lats = ds['latitude']
    lons = ds['longitude']
    LON, LAT = np.meshgrid(lons, lats)  # get grid of lon/lats
    xs, ys = proj(LON, LAT)  # project
    dxs = xs[:,1:] - xs[:,:-1]
    dys = ys[1:,:] - ys[:-1,:]
    dA = op.resize(dxs, 0) * op.resize(dys, 1)  # area on vort grid, meters
    xeigs = op.resize(op.resize(xs, -1), -1)
    yeigs = op.resize(op.resize(ys, -1), -1)
    loneigs, lateigs = proj(xeigs, yeigs, inverse=True)
    t = ds['time']
    sshmax = 0.5; sshmin = -sshmax
    dt = 1  # use every time output
    # for setting up final dataframe and selecting times
    # dates = pd.to_datetime(ds['time'].values)
    # dfdates = pd.DataFrame(index=dates)
elif whichmodel == 'hycom':
    loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
    ds = xr.open_dataset(loc)
    grid_file = 'gom03_grd_N050_new.nc'
    grid = xr.open_dataset(grid_file)
    LON = grid.lon_rho.values
    LAT = grid.lat_rho.values
    xs, ys = proj(LON, LAT)  # project
    xsm1 = xs[1:-1,1:-1]
    ysm1 = ys[1:-1,1:-1]
    dxs = xs[:,1:] - xs[:,:-1]
    dys = ys[1:,:] - ys[:-1,:]
    # dA = op.resize(dxs, 0) * op.resize(dys, 1)  # area on vort grid, meters
    # dA = (xsm1[:,1:] - xsm1[:,:-1])*(ysm1[1:,:] - ysm1[:-1,:])
    dx = 1/grid.pm[1:-1,1:-1].values
    dx = op.resize(op.resize(dx, 1), 0)
    dy = 1/grid.pn[1:-1,1:-1].values
    dy = op.resize(op.resize(dy, 1), 0)
    dA = dx*dy
    # xeigs = op.resize(op.resize(xs, -1), -1)
    # yeigs = op.resize(op.resize(ys, -1), -1)
    # loneigs, lateigs = proj(xeigs, yeigs, inverse=True)
    loneigs = grid.lon_psi[1:-1,1:-1].values
    lateigs = grid.lat_psi[1:-1,1:-1].values
    xeigs, yeigs = proj(loneigs, lateigs)
    xu, yu = proj(grid.lon_u.values, grid.lat_u.values)
    xv, yv = proj(grid.lon_v.values, grid.lat_v.values)
    t = ds['ocean_time']
    sshmax = 1.0; sshmin = -sshmax
    dt = 8  # use every 8 time outputs (24 hours)
    dtname = '24H'
    filtersize = 9
    # assign ocean_time to coordinate
    ds = ds.assign_coords(ocean_time=ds['ocean_time'])
    # for setting up final dataframe and selecting times
    # dates = pd.to_datetime(ds['ocean_time'].values)
    # dfdates = pd.DataFrame(index=dates)

def dofilter(up, vp, size=3):
    # filter velocities a little
    up0 = up.copy()
    up0[np.isnan(up0)] = 0  # fill nan's with 0s to use in filter
    upf = scipy.ndimage.uniform_filter(up0, size=size)
    vp0 = vp.copy()
    vp0[np.isnan(vp0)] = 0  # fill nan's with 0s to use in filter
    vpf = scipy.ndimage.uniform_filter(vp0, size=size)
    # reinstate nan's
    up = upf.copy()
    up[up == 0] = np.nan
    vp = vpf.copy()
    vp[vp == 0] = np.nan

    return up, vp

# use pts to calculate circulation later
# have a Point for each grid node in vortex calculation
pts = [shapely.geometry.Point(x,y) for (x,y) in zip(xeigs.flatten(), yeigs.flatten())]

if plotdrifters:
    # load in drifter data
    dfmix = pd.read_csv('drifter_animation/drifter-09_04_18-07_37.csv', usecols=[0,1,4,5], parse_dates=True, index_col=1)
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

# east boundary
xeast, yeast = proj(-82,26)
# south boundary
xsouth, ysouth = proj(-88,19)

# plot prep
land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-90.0)
# extent = [-98, -85, 22, 31]
extent = [-98, -80, 18, 31]
aed = cartopy.crs.AzimuthalEquidistant(central_longitude=-90.0, central_latitude=26)

base = 'figures/swirl/%s' % whichmodel
if plotdrifters:
    base += '/drifters'
if spacefilter:
    base += '/spacefilter_size%i' % filtersize
if timeaverage:
    base += '/timeaverage_dt%i' % dt
os.makedirs(base, exist_ok=True)

ptsname = 'calcs/swirl/'
if plotdrifters:
    ptsname += '/drifters'
if spacefilter:
    ptsname += '/spacefilter_size%i' % filtersize
if timeaverage:
    ptsname += '/timeaverage_dt%i' % dt
os.makedirs(ptsname, exist_ok=True)
# ptsname += '/points%i.csv' % year
ptsname += '/points.csv'

# set up dataframe to store polygon points

# index = pd.date_range(start=dfdates[str(year)].index[0],
#                       end=dfdates[str(year)].index[-1], freq=dtname)
df = pd.DataFrame(dtype=object, columns=['x','y'])
# df = pd.DataFrame(dtype=object, columns=['p'])


# plot
pathssave = []  # save paths over multiple days
for i in range(0,t.size, dt):
# i = 7

    date = t[i]
    fname = '%s/%s.png' % (base,str(date.values)[:13])

    if os.path.exists(fname):
        continue

    if whichmodel == 'mercator':
        # use nearest model output
        ssh = ds['zos'].isel(time=i).values
        # ssh resized in a grid
        sshpsi = op.resize(op.resize(ssh, 1), 0)
        u = ds['uo'].isel(depth=0, time=i).values  # y by x
        v = ds['vo'].isel(depth=0, time=i).values

        if spacefilter:
            u, v = dofilter(u, v, filtersize)

        ## Swirl strength ##
        # calculate velocity gradient
        dudx = op.resize((u[:,1:] - u[:,:-1])/dxs, 0)
        dudy = op.resize((u[1:,:] - u[:-1,:])/dys, 1)
        dvdx = op.resize((v[:,1:] - v[:,:-1])/dxs, 0)
        dvdy = op.resize((v[1:,:] - v[:-1,:])/dys, 1)

        # calculate velocity gradient matrix G
        G = np.array([[dudx, dudy], [dvdx, dvdy]])

        # want the positive imaginary eigenvalue of matrix G
        Gtemp = G.T.copy()
        inds = np.isnan(Gtemp)  # save nans to remove values later
        Gtemp[inds] = 0  # can't have nans in matrix
        eigs = np.imag(np.linalg.eig(Gtemp)[0])
        # keep first eigenvalue, which is the first one
        eigs = eigs[:,:,0]
        # if the sum of the inds over a 2x2 matrix is > 0, there was a nan in it
        inans = inds.sum(axis=-1).sum(axis=-1) > 0
        eigs[inans] = np.nan
        # transpose back
        eigs = eigs.T

        ## vertical vorticity ##
        vort = dvdx - dudy

        ## swirl with direction ##
        swirl = eigs*np.sign(vort)

    elif whichmodel == 'hycom':

        if timeaverage:
            ssh = ds['zeta'].isel(time_1=slice(i,i+dt)).values
            ssh = ssh.mean(axis=0)  # time mean
            # ssh resized to psi grid minus 1
            sshpsi = op.resize(op.resize(ssh[1:-1,1:-1], 1), 0)
            u = ds['u'].isel(time_1=slice(i,i+dt)).values  # y by x
            v = ds['v'].isel(time_1=slice(i,i+dt)).values
            u = u.mean(axis=0)  # time mean
            v = v.mean(axis=0)  # time mean
        else:
            ssh = ds['zeta'].isel(time_1=i).values
            # ssh resized to psi grid minus 1
            sshpsi = op.resize(op.resize(ssh[1:-1,1:-1], 1), 0)
            u = ds['u'].isel(time_1=i).values  # y by x
            v = ds['v'].isel(time_1=i).values

        if spacefilter:
            u, v = dofilter(u, v, filtersize)

        ## Swirl strength ##
        # calculate velocity gradient
        # get onto inner psi grid
        dudx = (u[:,1:] - u[:,:-1])/(xu[:,1:] - xu[:,:-1])
        dudx = op.resize(op.resize(dudx, 1), 0)
        dudx = dudx[1:-1,:]
        dudy = (u[1:,:] - u[:-1,:])/(yu[1:,:] - yu[:-1,:])
        dudy = dudy[1:-1,1:-1]

        dvdx = (v[:,1:] - v[:,:-1])/(xv[:,1:] - xv[:,:-1])
        dvdx = dvdx[1:-1,1:-1]
        dvdy = (v[1:,:] - v[:-1,:])/(yv[1:,:] - yv[:-1,:])
        dvdy = op.resize(op.resize(dvdy, 0), 1)
        dvdy = dvdy[:,1:-1]

        # calculate velocity gradient matrix G
        G = np.array([[dudx, dudy], [dvdx, dvdy]])

        # want the positive imaginary eigenvalue of matrix G
        Gtemp = G.T.copy()
        inds = np.isnan(Gtemp)  # save nans to remove values later
        Gtemp[inds] = 0  # can't have nans in matrix
        eigs = np.imag(np.linalg.eig(Gtemp)[0])
        # keep first eigenvalue, which is the first one
        eigs = eigs[:,:,0]
        # if the sum of the inds over a 2x2 matrix is > 0, there was a nan in it
        inans = inds.sum(axis=-1).sum(axis=-1) > 0
        eigs[inans] = np.nan
        # transpose back
        eigs = eigs.T

        ## vertical vorticity, inner psi grid ##
        vort = dvdx - dudy

        ## swirl with direction ##
        swirl = eigs*np.sign(vort)



    # contour at a minimum value of swirl to find isosurfaces
    cs = plt.contour(xeigs, yeigs, swirl, [0.5e-6])
    # cs = plt.contour(xeigs, yeigs, eigs[i], [3e-6], transform=aed)
    plt.close(plt.gcf())


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


    # # eigs is swirl without direction
    # mappable = ax.pcolormesh(loneigs, lateigs, eigs, transform=pc,
    #                          cmap=cmo.speed, vmin=0, vmax=3e-5)
    # # swirl has direction according to vorticity
    # mappable = ax.pcolormesh(loneigs, lateigs, swirl, transform=pc,
    #                          cmap=cmo.curl, vmin=-3e-5, vmax=3e-5)
    if whichmodel == 'mercator':
        # sea surface height
        mappable = ax.pcolormesh(lons, lats, ssh, transform=pc,
                                 cmap=cmo.balance, vmin=sshmin, vmax=sshmax)
    elif whichmodel == 'hycom':
        # sea surface height
        mappable = ax.pcolormesh(LON, LAT, ssh, transform=pc,
                                 cmap=cmo.balance, vmin=sshmin, vmax=sshmax)
        # # eigs is swirl without direction
        # mappable = ax.pcolormesh(loneigs, lateigs, eigs, transform=pc,
        #                          cmap=cmo.speed, vmin=0, vmax=3e-5)
        # # swirl has direction according to vorticity
        # mappable = ax.pcolormesh(loneigs, lateigs, swirl, transform=pc,
        #                          cmap=cmo.curl, vmin=-3e-5, vmax=3e-5)

    # overlay velocities
    if whichmodel == 'mercator':
        dd = 5
        ax.quiver(lons[::dd].values, lats[::dd].values, u[::dd,::dd], v[::dd,::dd], transform=pc, alpha=0.6);
    elif whichmodel == 'hycom':
        dd = 15
        ax.quiver(LON[::dd,::dd], LAT[::dd,::dd], u[::dd,::dd], v[::dd,::dd], transform=pc, alpha=0.6);

    cb = fig.colorbar(mappable, shrink=0.8, extend='both')
    cb.set_label('Sea surface height [m]')

    # loop through identified eddies from contours to calculate area and circulation
    # and eliminate bad choices based on these criteria
    paths = cs.collections[0].get_paths()

    # has to have enough points to be big enough
    paths2 = [path for path in paths if path.vertices.size>90]

    xsave = []  # to store accepted polygon coordinates
    ysave = []
    # psave = []
    paths3 = []
    for path in paths2:
    # path = paths[0]

        # has to have an area to count
        try:
            p = shapely.geometry.Polygon(path.vertices)
            assert p.area
        except:
            continue


        # eliminate if east of east boundary or south of south boundary
        if (p.centroid.xy[0][0] > xeast) or (p.centroid.xy[1][0] < ysouth):
            continue

        # calculate equivalent radius for polygon area (km)
        R = np.sqrt(p.area/1000**2/np.pi)

        # radius should be right size to be an eddy (30-140km?)
        # if not ((R>30) and (R<140)):
        if R < 40:# or R > 100:
            continue

    # # test plot
    # # plt.contour(xeigs, yeigs, swirl[i], [3e-6], transform=aed)


        ## circulation strength ##
        # boolean array of whether vortex x,y location is inside polygon p
        contains = [p.contains(pt) for pt in pts]
        # reshape to match vort array
        contains = np.asarray(contains).reshape(xeigs.shape)
        # add up vorticity and multiply by area for grid nodes for which contains==True
        circ = (vort[contains]*dA[contains]).sum()  # m^2/s

        # need strong enough circulation
        if circ < 0.7e5:
            continue

        # # need right size aspect ratio of bounds
        # W = abs(p.bounds[0] - p.bounds[2])
        # H = abs(p.bounds[1] - p.bounds[3])
        # ratio = W/H
        # if (ratio < 0.65) or (ratio > 1.5):
        #     continue

        # skip things that are above MSL — not sure what to think about them
        if sshpsi[contains].sum() > -4:
            continue

        paths3.append(path)

        # ax.add_geometries([p.convex_hull], crs=aed, facecolor='none', edgecolor='g')

    # compare neighboring eddies to see if they should be combined
    k = 0
    while k < len(paths3)-1:
        p0 = shapely.geometry.Polygon(paths3[k].vertices)
        p1 = shapely.geometry.Polygon(paths3[k+1].vertices)

        # if they overlap, combine and remove previous path
        buffer = 20000
        if p0.buffer(buffer).intersects(p1):
            paths3[k] = p0.buffer(buffer).union(p1).convex_hull
            paths3.pop(k+1)
        k += 1

    # save pointsmaking up polygons
    for path in paths3:

        if not isinstance(path, shapely.geometry.polygon.Polygon):
            p = shapely.geometry.Polygon(path.vertices)
        else:  # is already a polygon
            p = path
        # p = shapely.geometry.Polygon(path.vertices)

        # save convex hull since don't want swirl contours
        p = p.convex_hull

        # reduce number of points used to define polygon
        p = p.simplify(10000)
        # save whatever polygon points make it to this point
        xsave.append(list(p.exterior.xy[0]))
        ysave.append(list(p.exterior.xy[1]))
        # psave.append(p)

    # put into dataframe at corresponding time
    # df.loc[pd.Timestamp(ds['ocean_time'][i].values),'p'] = psave
    df.loc[pd.Timestamp(ds['ocean_time'][i].values),'x'] = xsave
    df.loc[pd.Timestamp(ds['ocean_time'][i].values),'y'] = ysave
    df.to_csv(ptsname)


    # add paths3 to pathssave now
    pathssave.insert(0, paths3)
    # if we have more than 7 days' worth, remove oldest one
    if len(pathssave) > 7:
        pathssave.pop(-1)
    alphas = np.linspace(1, 0.2, 7)

    # add eddy identifications to plot
    for paths3, alpha in zip(pathssave, alphas):  # loop over days of paths
        for path in paths3:  # loop over paths for a single day
            if not isinstance(path, shapely.geometry.polygon.Polygon):
                p = shapely.geometry.Polygon(path.vertices)
            else:
                p = path
            ax.add_geometries([p.convex_hull], crs=aed, facecolor='none',
                              edgecolor='k', linewidth=3, alpha=alpha)

    # overlay drifters
    if plotdrifters:
        for j in range(3):
            try:
                lat, lon = dfs[j].loc[pd.to_datetime(t[i].values)]
                ax.plot(lon, lat, 'ok', transform=pc)
            except:
                pass
        # add labels
        # drifter legend
        ax2 = fig.add_axes([0.125, 0.625, 0.15, 0.2], frameon=False)
        ax2.scatter([], [], c='0.2', s=30, marker='o', label='Drifter')
        ax2.legend(scatterpoints=1, frameon=False, loc='upper left')
        ax2.set_axis_off()

    # label model output
    ax.text(0.3, 0.96, '%s model output' % whichmodel, color='0.2', transform=ax.transAxes)
    # my name
    ax.text(0.8, 0.96, 'K. Thyng', color='0.2', transform=ax.transAxes)

    ax.set_title(pd.to_datetime(date.values).strftime('%b %d %H:%M, %Y'))

    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
