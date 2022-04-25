# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 07:01:41 2022

@author: Wenbo.Duan
"""

import os
import sys
import numpy as np
import netCDF4 as nc


from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

import Utilities.nctools as nctools



def loadFilesFromPath(inputPath, tilelimits):
    """
    Load wind field data for each subset into a 3-D array.

    :param str inputPath: str path to wind field files.

    :param tuple tilelimits: tuple of index limits of a tile.

    :returns: 3-D `numpy.narray` of wind field records.

    """

    fileList = os.listdir(inputPath)
    files = [pjoin(inputPath, f) for f in fileList]
    files = [f for f in files if os.path.isfile(f)]
    # log.debug("Loading data from %d files" % (len(files)))

    ysize = tilelimits[3] - tilelimits[2]
    xsize = tilelimits[1] - tilelimits[0]
    Vr = np.empty((len(files), ysize, xsize), dtype='f')

    for n, f in enumerate(sorted(files)):
        Vr[n,:,:] = loadFile(f, tilelimits)

    return Vr

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
        #     # log.debug("max tile limits are not smaller than min")
            
        return data_subset

    except IOError:
        # log.debug('{0} file does not exist'.format(filename))
        raise
        
windfieldpath= 'G:/Advanced Technology & Research/Climate Analysis/079007-70 IiA Extreme Winds/tcrm/output/Wilmington_10by10_5000yr_1980_run3/windfield'
tilelimits1=(0, 401, 0,401)
windfiled = loadFilesFromPath(windfieldpath, tilelimits1)


Target_location = [282.0,34.0] # lon, lat Wilmington
# Target_location = [289.0,42.4] # lon, lat Boston

# grid_limit = [279, 299, 31.5, 51.5] # lon min, lon max, lat min, lat max Boston
grid_limit = [272, 292, 24, 44] # lon min, lon max, lat min, lat max Boston
grid_step = 0.05
lon_grid = np.linspace(grid_limit[0],grid_limit[1],int((grid_limit[1]-grid_limit[0])/grid_step)+1)
lat_grid = np.linspace(grid_limit[2],grid_limit[3],int((grid_limit[3]-grid_limit[2])/grid_step)+1)

ind_target_lon=np.where(lon_grid == Target_location[0])
ind_target_lat=np.where(lat_grid == Target_location[1])

# Target_location_ind = [200,238]

winddata_target = windfiled[:,ind_target_lon,ind_target_lat]
winddata_target_list = winddata_target[:,0,0]

# # create NC file
# fn = 'C:/Users/wenbo.duan/OneDrive - Arup/03_Project/IiA/079007-70_Climate_change_Extreme_wind/TCRM note/Boston_10by10_5000yr_1980_Windfield.nc'
# ds = nc.Dataset(fn, 'w', format='NETCDF4')
# time = ds.createDimension('time', None)
# lat = ds.createDimension('lat', 400)
# lon = ds.createDimension('lon', 400)
# times = ds.createVariable('time', 'f4', ('time',))
# lats = ds.createVariable('lat', 'f4', ('lat',))
# lons = ds.createVariable('lon', 'f4', ('lon',))
# value = ds.createVariable('slp', 'f4', ('time', 'lat', 'lon',))
# value.units = 'm/s'
# value.description = fn
# value.long_name = 'sea level pressure' # add attributes
# times.long_name = 'Time'
# lats.long_name = 'Latitude'
# lons.long_name = 'Longitude'
# lats[:] = lat_grid
# lons[:] = lon_grid

# value[:, :, :] = windfiled
# ds.close()
