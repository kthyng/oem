'''
Script to run drifters at 1km initial spacing daily forward for 30 days.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import matplotlib.pyplot as plt
import tracpy
from datetime import datetime, timedelta
import glob
from tracpy.tracpy_class import Tracpy
import plots
import sinks

loc = 'http://terrebonne.tamu.edu:8080/thredds/dodsC/NcML/gom_roms_hycom'
grid_file = 'gom03_grd_N050_new.nc'
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid(grid_file, proj=proj)


def init(name, lonsink, latsink, sinkarrows):
    '''
    Initialization for the simulation.
    '''

    time_units = 'seconds since 1970-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 5 # interpolation forced between model output

    # Number of steps to divide model output for outputting drifter location
    N = 1

    # Number of days
    ndays = 30

    # This is a forward-moving simulation
    ff = 1

    # Time between outputs
    tseas = 3*3600 # 3 hours between outputs, in seconds, time between model outputs
    ah = 0. # old tracks: 5.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'
    zpar = 0
    zparuv = 0

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines.
    dostream = 0

    # Initialize Tracpy class
    tp = Tracpy(loc, grid=grid, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, dostream=dostream, savell=False, doperiodic=0,
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units, sinkarrows=sinkarrows)

    # initial separation distance of drifters, in meters, from sensitivity project
    dx = 100
    seedsfile = 'calcs/seeds_lon0_%2.2f_lat0_%2.2f.npz' % (abs(lonsink), latsink)
    if os.path.exists(seedsfile):
        seeds = np.load(seedsfile)
        lon0 = seeds['lon0']; lat0 = seeds['lat0']
        seeds.close()
    else:
        # Initial lon/lat locations for drifters
        # Start uniform array of drifters across domain using x,y coords
        llcrnrlon = lonsink - 0.1
        urcrnrlon = lonsink + 0.1
        llcrnrlat = latsink - 0.1
        urcrnrlat = latsink + 0.1
        xcrnrs, ycrnrs = tp.grid.proj([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
        X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], dx), np.arange(ycrnrs[0], ycrnrs[1], dx))
        lon0, lat0 = tp.grid.proj(X, Y, inverse=True)

        # save starting locations for future use
        np.savez(seedsfile, lon0=lon0, lat0=lat0)

    return tp, lon0, lat0


def run():

    # Make sure necessary directories exist
    os.makedirs('tracks', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    # speeds to use
    speeds = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    # center of starting locations/sink locations (same thing)
    sinklocs = np.array([[-95, 26.75], [-91, 26.75], [-88, 27.5], [-85.5, 25.5]])

    # plot sinks and region of interest
    plots.roi(grid, sinks=sinklocs)

    overallstartdate = datetime(2012, 1, 1, 0, 0)
    overallstopdate = datetime(2013, 1, 1, 0, 0)

    date = overallstartdate
    # Start from the beginning and add days on for loop
    # keep running until we hit the next month
    while date < overallstopdate:

        for sinkloc in sinklocs:

            lonsink, latsink = sinkloc

            for speed in speeds:

                # create sink velocity file for tracpy simulation
                iu, jv = sinks.create(lonsink, latsink, speed)

                basedir = date.isoformat()[0:13]
                os.makedirs('tracks/' + basedir, exist_ok=True)
                name = 'lon0_%2.2f_lat0_%2.2f_s_%2.2f' % (abs(sinkloc[0]), sinkloc[1], speed)

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc') and \
                    not os.path.exists('tracks/' + name + 'gc.nc'):

                    # Read in simulation initialization
                    tp, lon0, lat0 = init(basedir + '/' + name, lonsink, latsink, sinkarrows=[iu, jv])

                    # plot sink location, region of interest, arrows, and seed locations
                    plots.roi(grid, sinks=sinkloc, sinkarrows=[iu, jv],
                              seeds=[lon0, lat0])

                    # Run tracpy
                    lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

        # Increment by delta time for next loop, to move through more quickly
        date = date + timedelta(days=7)


if __name__ == "__main__":
    run()
