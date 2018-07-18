import shapely.geometry
import matplotlib.pyplot as plt
import numpy as np
import tracpy
import xarray as xr
import cmocean.cm as cmo
import pandas as pd

loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
ds = xr.open_dataset(loc)
grid_file = '../experiments/gom03_grd_N050_new.nc'
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid(grid_file, proj=proj)

# initialize vec (inner psi grid)
vec = np.ones((5,*grid.x_psi[1:-1,1:-1].shape))

# set up dataframe
index = pd.date_range(start=ds['ocean_time'][0].values, end=ds['ocean_time'][-1].values, freq='3H')
df = pd.DataFrame(index=index, dtype=object, columns=['Points'])

for i in range(6):
    v = ds['v'][i].data  # v grid
    u = ds['u'][i].data  # u grid
    up = tracpy.op.resize(u, 0)  # psi grid
    vp = tracpy.op.resize(v, 1)  # psi grid

    # # filter velocities a little
    # up0 = up.copy()
    # up0[np.isnan(up0)] = 0  # fill nan's with 0s to use in filter
    # upf = scipy.ndimage.uniform_filter(up0, size=3)
    # vp0 = vp.copy()
    # vp0[np.isnan(vp0)] = 0  # fill nan's with 0s to use in filter
    # vpf = scipy.ndimage.uniform_filter(vp0, size=3)
    # # reinstate nan's
    # up = upf.copy()
    # up[up == 0] = np.nan
    # vp = vpf.copy()
    # vp[vp == 0] = np.nan

    # after vec array is initially filled, always add to last time entry
    if i>4:
        ii=4
    else:
        ii=i

    # get first calculation of zeros
    bottom = up[:-2,1:-1]  # inner psi grid, 1:-1,1:-1
    top = up[2:,1:-1]  # inner psi grid
    left = vp[1:-1,:-2]  # inner psi grid
    right = vp[1:-1,2:]  # inner psi grid
    test1 = (np.sign(bottom) + np.sign(top))  # inner psi grid
    test2 = (np.sign(left) + np.sign(right))  # inner psi grid

    # check not shear: these should equal 0
    test3 = np.sign(bottom) + np.sign(left)  # inner psi grid

    # check correct rotation direction (ccw): this should be positive
    test4 = np.sign(right)  # inner psi grid

    # set up a grid of vortex points instead of the where statement
    # inner psi grid
    vec[ii] = (~((test1==0) * (test2 == 0) * (test3 == 0) * (test4 > 0))).astype(int)

    # limit vec entries
    # too close to mask if within 40 cells in any direction
    # inner psi grid
    nearland = scipy.ndimage.uniform_filter(grid.mask_psi[1:-1,1:-1], size=40)

    # make all vec points near land 1 (not 0)
    # inner psi grid
    vec[ii][~np.isclose(nearland, 1)] = 1

    # also remove if outside a generous north GOM box
    # must be north of 24 and west of 84
    # inner psi grid
    vec[ii][grid.lon_psi[1:-1,1:-1] > -84] = 1  # remove pts east of 84W

    # remove points south of 24 N
    # inner psi grid
    vec[ii][grid.lat_psi[1:-1,1:-1] < 24] = 1

    # izeros is on inner psi grid
    izeros = np.where(vec[ii]==0)

    # number of zeros
    nzeros = izeros[0].size

    # initial value of newzeros
    nnewzeros = nzeros + 1
    atol = 0.005

    # strength is the magnitude of the velocity components
    # calculated for everywhere but only applies for the indices in np.where statement
    # inner psi grid
    strength0 = abs(up[:-2,1:-1]) + abs(up[2:,1:-1]) + abs(vp[1:-1,:-2]) + abs(vp[1:-1,2:])
    strength = strength0.copy()

    # keep running while new points being found
    while nnewzeros > nzeros:
        # update izeros (on inner psi grid)
        izeros = np.where(vec[ii]==0)
        # update number of zeros
        nzeros = izeros[0].size

        # select neighbors to 0s which are potential new 0s
        # i0 = izeros[0]-1 + izeros[0]+1 + izeros[0] + izeros[0]
        # these are all on inner psi grid
        i0 = np.concatenate((izeros[0]-1, izeros[0]+1, izeros[0], izeros[0]))
        i1 = np.concatenate((izeros[1], izeros[1], izeros[1]-1, izeros[1]+1))
        inewzeros = (i0,i1)

        # remove old zeros from new zeros with set, then put back into similar form
        # inewzeros is on inner psi grid
        inewzeros = set(zip(*inewzeros)) - set(zip(*izeros))

        # loop over new possible zeros to see if they fit
        for i0, i1 in inewzeros:

            # find where 0 is in relation to possible new vort pt
            if vec[ii][i0-1,i1] == 0:  # check for zero below i0, i1
                # find last 1, which would be closest to ind
                # this is the index for the bottom velocity to us
                ibottom = np.where(vec[ii][:i0,i1]==1)[0][-1]
                bottom = up[1:-1,1:-1][ibottom,i1]
            else:
                # bottom = up[1:-1,1:-1][i0-4:i0-1,i1]
                bottom = up[1:-1,1:-1][i0-1,i1]

            if vec[ii][i0+1,i1] == 0:  # check for zero above i0, i1
                itop = np.where(vec[ii][i0+1:,i1]==1)[0][0]
                top = up[1:-1,1:-1][i0+1+itop,i1]
            else:
                # top = up[1:-1,1:-1][i0+1:i0+4,i1]
                top = up[1:-1,1:-1][i0+1,i1]

            if vec[ii][i0,i1-1] == 0:  # check for zero to left of i0,i1
                ileft = np.where(vec[ii][i0,:i1]==1)[0][-1]
                left = vp[1:-1,1:-1][i0,ileft]
            else:
                left = vp[1:-1,1:-1][i0,i1-1]
                # left = vp[1:-1,1:-1][i0,i1-4:i1-1]

            if (vec[ii][i0,i1+1] == 0):  # check for zero to right of i0,i1 (possible new zero)
                iright = np.where(vec[ii][i0,i1+1:]==1)[0][0]
                right = vp[1:-1,1:-1][i0,i1+1+iright]
            else:
                right = vp[1:-1,1:-1][i0,i1+1]
                # right = vp[1:-1,1:-1][i0,i1+1:i1+4]

            # redo tests for this possible new zero
            # CHANGE TESTS TO BE FOR SPECIFIC SIGN AND CHECK MORE THAN ONE VALUE
            # test1 = ((bottom>0).any() & (top<0).any())  # inner psi grid
            # test1 = (np.sign(bottom) + np.sign(top))  # inner psi grid
            test1 = (np.isclose(0,bottom,atol=atol) or (bottom>0)) & \
                     (np.isclose(0,top,atol=atol) or (top<0))  # inner psi grid
            # test2 = ((left<0).any() & (right>0).any())  # inner psi grid
            # test2 = (np.sign(left) + np.sign(right))  # inner psi grid
            test2 = (np.isclose(0,left,atol=atol) or (left<0)) & \
                     (np.isclose(0,right,atol=atol) or (right>0))  # inner psi grid

            # check not shear: these should equal 0
            # test3 = (bottom>0).any() & (left<0).any()  # inner psi grid
            # test3 = np.sign(bottom) + np.sign(left)  # inner psi grid
            test3 = (np.isclose(0,bottom,atol=atol) or (bottom>0)) & \
                     (np.isclose(0,left,atol=atol) or (left<0))  # inner psi grid

            # check correct rotation direction (ccw): this should be positive
            # test4 = (right>0).any()  # inner psi grid
            # test4 = np.sign(right)  # inner psi grid
            test4 = (right>0) or (left<0) or (top<0) or (bottom>0)  # inner psi grid

            if test1 and test2 and test3 and test4:
            # if (test1 == 0) and (test2 == 0) and (test3 == 0) and (test4 > 0):
                vec[ii][i0,i1] = 0

        # update inewzeros with new for sure zeros from the for loop
        inewzeros = np.where(vec[ii]==0)

        # update value of newzeros with number of new points
        nnewzeros = inewzeros[0].size

        strength[vec[ii] == 0] += strength0[vec[ii] == 0]

    # sum over 5 time steps (15 hours)
    if i>4:
        # add across time steps and only save when always 0
        vecsum = np.ones(vec[0].shape)
        vecsum[vec.sum(axis=0) == 0] = 0
        # roll vec forward so can stick in next loop's vec entry to final time location
        vec = np.roll(vec, -1, axis=0)

        # collect contours
        cs = contour(grid.x_psi[1:-1,1:-1], grid.y_psi[1:-1,1:-1], vecsum, [0])
        plt.close(plt.gcf())

        # only keep cyclone size range cyclones
        paths = cs.collections[0].get_paths()

        psave = []  # to store accepted polygon coordinate tuples

        for path in paths:
        # path = paths[0]
            try:
                p = shapely.geometry.Polygon(path.vertices)
                assert p.area
            except:
                continue

            # calculate equivalent radius for polygon area (km)
            R = np.sqrt(p.area/1000**2/np.pi)

            # radius should be right size to be an eddy (30-140km?)
            if not (R>30) and (R<140):
                continue

            # nearest distance between centroid and exterior of polygon (km)
            D = p.exterior.distance(p.centroid)/1000

            # ratio of radius and centroid to boundary distance should be around 1
            ratio = R/D
            if not (ratio > 0.1) and (ratio < 2):
                continue
            #
            # skewed: 0.07 , 3.5
            # ok: 0.2, .7, .8

            # see if second to last pt is near last pt to see if
            # boundary was closed by shapely (km)
            x, y = p.exterior.xy
            p1 = shapely.geometry.Point(x[-1], y[-1])
            p2 = shapely.geometry.Point(x[-2], y[-2])
            if p1.distance(p2)/1000 > 20:
                continue
            #
            # too high: 150, 260, 32?
            # ok: 0, 3
            #

            # save whatever polygon points make it to this point
            # reduce number of points used to define polygon
            psave.append(p.simplify(10000).exterior.xy)

        # put into dataframe at corresponding time
        df.loc[pd.Timestamp(ds['ocean_time'][i].values),'Points'] = psave
        df.to_csv('points.csv')

    # figure()
    # pcolormesh(grid.lon_rho[1:-1,1:-1],
    #        grid.lat_rho[1:-1,1:-1], vec)
    # contour(grid.lon_psi[1:-1,1:-1],
    #        grid.lat_psi[1:-1,1:-1], vec, [1], colors=colors[i])
# quiver(grid.lon_psi[1:-1,1:-1][::10,::10],
#        grid.lat_psi[1:-1,1:-1][::10,::10],
#        up[1:-1,1:-1][::10,::10],
#        vp[1:-1,1:-1][::10,::10])

# # check another
# figure()
# pcolormesh(grid.lon_rho[1:-1,1:-1][266:286,129:149],
#        grid.lat_rho[1:-1,1:-1][266:286,129:149], vec[266:286,129:149])
# quiver(grid.lon_psi[1:-1,1:-1][266:286,129:149][::1,::1],
#        grid.lat_psi[1:-1,1:-1][266:286,129:149][::1,::1],
#        up[1:-1,1:-1][266:286,129:149][::1,::1],
#        vp[1:-1,1:-1][266:286,129:149][::1,::1])
