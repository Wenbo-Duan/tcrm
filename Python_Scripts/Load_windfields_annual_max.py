# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:07:01 2022

@author: Wenbo.Duan
"""
import os
import sys
import numpy as np
import netCDF4 as nc

import Utilities.lmomentFit as lmom

import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

import Utilities.nctools as nctools


def aggregateWindFields(inputPath, numSimulations, tilelimits):
    """
    Aggregate wind field data into annual maxima for use in fitting
    extreme value distributions.

    :param str inputPath: path to individual wind field files.
    :param int numSimulations: Number of simulated years of activity.

    """
    from glob import glob
    # log.info("Aggregating individual events to annual maxima")
    ysize = tilelimits[3] - tilelimits[2]
    xsize = tilelimits[1] - tilelimits[0]
    Vm = np.zeros((numSimulations, ysize, xsize), dtype='f')

    for year in range(numSimulations):
        filespec = pjoin(inputPath, "gust.*-%05d.nc"%year)
        fileList = glob(filespec)
        if len(fileList) == 0:
            # log.debug("No files for year: {0}".format(year))
            Vm[year, :, :] = np.zeros((ysize, xsize), dtype='f')
            continue

        Va = np.zeros((len(fileList), ysize, xsize), dtype='f')
        for n, f in enumerate(fileList):
            Va[n, :, :] = loadFile(f, tilelimits)

        Vm[year, :, :] = np.max(Va, axis=0)

    return Vm

def loadFile(filename, limits):
    """
    Load a subset of the data from the given file, with the extent
    of the subset specified in the `limits` tuple

    :param str filename: str full path to file to load.

    :param tuple limits: tuple of index limits of a tile.

    :returns: 2-D `numpy.ndarray` of wind speed values.

    """

    (xmin, xmax, ymin, ymax) = limits

    try:
        ncobj = nctools.ncLoadFile(filename)
        ncobj_vmax = nctools.ncGetVar(ncobj, 'vmax')
        data_subset = ncobj_vmax[ymin:ymax, xmin:xmax]
        ncobj.close()

        # if xmax < xmin or ymax < ymin:
        #     log.debug("max tile limits are not smaller than min")
            
        return data_subset

    except IOError:
        # log.debug('{0} file does not exist'.format(filename))
        raise
        
windfieldpath= 'G:/Advanced Technology & Research/Climate Analysis/079007-70 IiA Extreme Winds/tcrm/output/Wilmington_10by10_5000yr_1980_run3/windfield'
numSim = 5000 # number of simulation
tilelimits1=(0, 401, 0,401)

wind_annual_max = aggregateWindFields(windfieldpath, numSim, tilelimits1)

Target_location = [282.0,34.0]  # lon, lat Wilmington

grid_limit = [272, 292, 24, 44] # lon min, lon max, lat min, lat max Boston
grid_step = 0.05
lon_grid = np.linspace(grid_limit[0],grid_limit[1],int((grid_limit[1]-grid_limit[0])/grid_step)+1)
lat_grid = np.linspace(grid_limit[2],grid_limit[3],int((grid_limit[3]-grid_limit[2])/grid_step)+1)

ind_target_lon=np.where(lon_grid == Target_location[0])
ind_target_lat=np.where(lat_grid == Target_location[1])

# Target_location_ind = [200,238]

winddata_target = wind_annual_max[:,ind_target_lon,ind_target_lat]
winddata_target_list = winddata_target[:,0,0]



l1, l2, l3 = lmom.samlmu(winddata_target_list, 3)

t3=l3/l2

if (l2 <= 0.) or (np.abs(t3) >= 1.):
    # Reject points where the second l-moment is negative
    # or the ratio of the third to second is > 1, i.e. positive
    # skew.
    print("Invalid l-moments")
else:
    # Parameter estimation returns the location, scale and
    # shape parameters
    xmom = [l1, l2, t3]
    # Calculates the parameters of the distribution given its
    # L-moments. GEV distribution calculated.
    loc, scale, shp = np.array(lmom.pelgev(xmom))
    # We only store the values if the first parameter is
    # finite (i.e. the location parameter is finite)
    if not np.isfinite(loc):
        print("Invalid l-moments")
        


plt.scatter(range(len(winddata_target_list)), winddata_target_list)
plt.show()
positive_value = winddata_target_list[winddata_target_list>0]
plt.scatter(range(len(positive_value)), positive_value)
plt.show()