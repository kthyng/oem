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



# read in model output
loc = 'drifter_animation/output/*.nc'
ds = xr.open_mfdataset(loc)
lats = ds['latitude']
lons = ds['longitude']

# make projection
inputs = {'proj': 'aeqd', 'ellps': 'clrk66', 'datum': 'NAD27',
                  'lat_0': 26, 'lon_0': -90,
                  'x_0': 0, 'y_0': 0}
proj = Proj(**inputs)
# proj = tracpy.tools.make_proj('nwgom-pyproj')#, **{'llcrnrlat': 18, 'urcrnrlon': -80})

LON, LAT = np.meshgrid(lons, lats)  # get grid of lon/lats
xs, ys = proj(LON, LAT)  # project
dxs = xs[:,1:] - xs[:,:-1]
dys = ys[1:,:] - ys[:-1,:]
dA = op.resize(dxs, 0) * op.resize(dys, 1)  # area on vort grid, meters
xeigs = op.resize(op.resize(xs, -1), -1)
yeigs = op.resize(op.resize(ys, -1), -1)
loneigs, lateigs = proj(xeigs, yeigs, inverse=True)

# use pts to calculate circulation later
# have a Point for each grid node in vortex calculation
pts = [shapely.geometry.Point(x,y) for (x,y) in zip(xeigs.flatten(), yeigs.flatten())]


# east boundary
xeast, yeast = proj(-82,26)
# south boundary
xsouth, ysouth = proj(-88,19)

# use nearest model output
ssh = ds['zos'].values#.sel(time=index, method='nearest')
# ssh resized in a grid
sshpsi = op.resize(op.resize(ssh, 2), 1)
u = ds['uo'].isel(depth=0).values  # time by y by x
v = ds['vo'].isel(depth=0).values

## Swirl strength ##
# calculate velocity gradient
# dudx = u.diff('longitude')/dxs
# dudy = u.diff('latitude')/dys
# dvdx = op.resize((v[:,:,1:].values - v[:,:,:-1].values)/(xs[:,1:] - xs[:,:-1]), 1)
# dvdy = op.resize((v[:,1:,:].values - v[:,:-1,:].values)/(ys[1:,:] - ys[:-1,:]), 2)
dudx = op.resize((u[:,:,1:] - u[:,:,:-1])/dxs, 1)
dudy = op.resize((u[:,1:,:] - u[:,:-1,:])/dys, 2)
dvdx = op.resize((v[:,:,1:] - v[:,:,:-1])/dxs, 1)
dvdy = op.resize((v[:,1:,:] - v[:,:-1,:])/dys, 2)

# calculate velocity gradient matrix G
G = np.array([[dudx, dudy], [dvdx, dvdy]])

# want the positive imaginary eigenvalue of matrix G
Gtemp = G.T.copy()
inds = np.isnan(Gtemp)  # save nans to remove values later
Gtemp[inds] = 0  # can't have nans in matrix
eigs = np.imag(np.linalg.eig(Gtemp)[0])
# keep first eigenvalue, which is the first one
eigs = eigs[:,:,:,0]
# if the sum of the inds over a 2x2 matrix is > 0, there was a nan in it
inans = inds.sum(axis=-1).sum(axis=-1) > 0
eigs[inans] = np.nan
# transpose back
eigs = eigs.T

## vertical vorticity ##
vort = dvdx - dudy

## swirl with direction ##
swirl = eigs*np.sign(vort)



# plot prep
land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m')
pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-90.0)
# extent = [-98, -85, 22, 31]
extent = [-98, -80, 18, 31]
aed = cartopy.crs.AzimuthalEquidistant(central_longitude=-90.0, central_latitude=26)


# plot
pathssave = []  # save paths over multiple days
for i in range(ds['time'].size):
# i = 7

    date = ds['time'][i]

    fname = 'figures/swirl/%s.png' % (str(date.values)[:10])

    if os.path.exists(fname):
        continue

    # contour at a minimum value of swirl to find isosurfaces
    cs = plt.contour(xeigs, yeigs, swirl[i], [0.5e-6])
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
    # mappable = ax.pcolormesh(loneigs, lateigs, eigs[i], transform=pc,
    #                          cmap=cmo.speed, vmin=0, vmax=3e-5)
    # # swirl has direction according to vorticity
    # mappable = ax.pcolormesh(loneigs, lateigs, swirl[i], transform=pc,
    #                          cmap=cmo.curl, vmin=-3e-5, vmax=3e-5)
    # sea surface height
    mappable = ax.pcolormesh(lons, lats, ssh[i], transform=pc,
                             cmap=cmo.balance, vmin=-0.5, vmax=0.5)
    # overlay velocities
    dd = 5
    ax.quiver(lons[::dd].values, lats[::dd].values, u[i,::dd,::dd], v[i,::dd,::dd], transform=pc, alpha=0.6);

    cb = fig.colorbar(mappable, shrink=0.8, extend='both')
    cb.set_label('Sea surface height [m]')



    # loop through identified eddies from contours to calculate area and circulation
    # and eliminate bad choices based on these criteria
    paths = cs.collections[0].get_paths()

    # has to have enough points to be big enough
    paths2 = [path for path in paths if path.vertices.size>90]

    xsave = []  # to store accepted polygon coordinates
    ysave = []
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
        circ = (vort[i,contains]*dA[contains]).sum()  # m^2/s

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
        if sshpsi[i,contains].sum() > -4:
            continue

        paths3.append(path)

        # ax.add_geometries([p.convex_hull], crs=aed, facecolor='none', edgecolor='g')

        # # save convex hull since don't want swirl contours
        # p = p.convex_hull
        #
        # # reduce number of points used to define polygon
        # p = p.simplify(10000)
        # # save whatever polygon points make it to this point
        # xsave.append(list(p.exterior.xy[0]))
        # ysave.append(list(p.exterior.xy[1]))

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

    # add paths3 to pathssave now
    pathssave.append(paths3)
    # if we have more than 7 days' worth, remove oldest one
    if len(pathssave) > 7:
        pathssave.pop(0)
    alphas = np.linspace(1,0.2, 7)

    # add eddy identifications to plot
    for paths3, alpha in zip(pathssave, alphas):  # loop over days of paths
        for path in paths3:  # loop over paths for a single day
            if not isinstance(path, shapely.geometry.polygon.Polygon):
                p = shapely.geometry.Polygon(path.vertices)
            else:
                p = path
            ax.add_geometries([p.convex_hull], crs=aed, facecolor='none',
                              edgecolor='k', linewidth=3, alpha=alpha)


    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
