# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 21:31:14 2022

@author: Wenbo.Duan
"""
import Utilities.nctools as nctools
import netCDF4 as nc
import numpy as np
import os
import calendar
from datetime import date
import datetime

path = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\Climate_Model\\CMIP6 climate projections\\SSP4.5'
file = 'psl_day_HadGEM3-GC31-LL_ssp245_r1i1p1f3_gn_20600101-20601230_v20190908.nc'
ncobj = nctools.ncLoadFile(os.path.join(path,file))
data_2031_2060 = nctools.ncGetData(ncobj, 'psl')
data_lon = nctools.ncGetData(ncobj, 'lon')
data_lat = nctools.ncGetData(ncobj, 'lat')

# Numberofyear_2031_2060 = int(np.shape(data_2031_2060)[0] / 365)
# YearlyData_2031_2060 = np.empty((Numberofyear_2031_2060, 365, 64*128))
# for i in range(Numberofyear_2031_2060):
#     YearlyData_2031_2060[i,:,:] = data_2031_2060[i*365:365*i+365,:,:].reshape(365, 64*128)

# DailyMean_2031_2060 = np.mean(YearlyData_2031_2060,axis = 0)
# DailyMean_2031_2060_grid = DailyMean_2031_2060.reshape(365,64,128)
# DailyMean_2031_2060_grid = DailyMean_2031_2060_grid[: , ::-1,: ] # lat in decending order
# data_lat = data_lat[::-1 ]# lat in decending order

# # create NC file
# fn = 'C:/tcrm/MSLP/mslp.day.2031-2060.nc'
# ds = nc.Dataset(fn, 'w', format='NETCDF4')
# time = ds.createDimension('time', None)
# lat = ds.createDimension('lat', 64)
# lon = ds.createDimension('lon', 128)
# times = ds.createVariable('time', 'f4', ('time',))
# lats = ds.createVariable('lat', 'f4', ('lat',))
# lons = ds.createVariable('lon', 'f4', ('lon',))
# value = ds.createVariable('slp', 'f4', ('time', 'lat', 'lon',))
# value.units = 'Pascals'
# value.description = file
# value.long_name = 'sea level pressure' # add attributes
# times.long_name = 'Time'
# lats.long_name = 'Latitude'
# lons.long_name = 'Longitude'
# lats[:] = data_lat
# lons[:] = data_lon

# value[:, :, :] = DailyMean_2031_2060_grid

# ds.close()

# create NC file for 2060
DailyMean_2060_grid = data_2031_2060
# DailyMean_2060 = YearlyData_2031_2060[-1,:,:]
# DailyMean_2060_grid = DailyMean_2060.reshape(365,64,128)
DailyMean_2060_grid = DailyMean_2060_grid[: , ::-1,: ] # lat in decending order
data_lat = data_lat[::-1 ]# lat in decending order

fn = 'C:/tcrm/MSLP/mslp.day.2060.HadGEM3LL_SSP4.5.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
time = ds.createDimension('time', None)
lat = ds.createDimension('lat', len(data_lat))
lon = ds.createDimension('lon', len(data_lon))
times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('slp', 'f4', ('time', 'lat', 'lon',))
value.units = 'Pascals'
value.description = file
value.long_name = 'sea level pressure' # add attributes
times.long_name = 'Time'
lats.long_name = 'Latitude'
lons.long_name = 'Longitude'
lats[:] = data_lat
lons[:] = data_lon
value[:, :, :] = DailyMean_2060_grid
ds.close()


#load historical data from CanESM5
path_his = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\Climate_Model\\CMIP6 climate projections\\Daily-SLP-CanESM5_historical'
file_his = 'psl_day_CanESM5_historical_r1i1p1f1_gn_18500101-20141231_v20190429.nc'
ncobj_his = nctools.ncLoadFile(os.path.join(path_his,file_his))
data_his = nctools.ncGetData(ncobj_his, 'psl')
data_his_lon = nctools.ncGetData(ncobj_his, 'lon')
data_his_lat = nctools.ncGetData(ncobj_his, 'lat')


his_1981_2010 = data_his[(1981 - 1850) * 365 : (2011 - 1850) * 365,:,:]

Numberofyear_1981_2010 = int(np.shape(his_1981_2010)[0] / 365)
YearlyData_1981_2010 = np.empty((Numberofyear_1981_2010, 365, 64*128))
for i in range(Numberofyear_1981_2010):
    YearlyData_1981_2010[i,:,:] = his_1981_2010[i*365:365*i+365,:,:].reshape(365, 64*128)

DailyMean_1981_2010 = np.mean(YearlyData_1981_2010,axis = 0)
DailyMean_1981_2010_grid = DailyMean_1981_2010.reshape(365,64,128)
DailyMean_1981_2010_grid = DailyMean_1981_2010_grid[: , ::-1,: ] # lat in decending order
data_his_lat = data_his_lat[::-1 ]# lat in decending order

Delta_slp = DailyMean_2031_2060_grid - DailyMean_1981_2010_grid

# create NC file
fn = 'C:/tcrm/MSLP/Deltamslp.day.2031-2060.1981-2010.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
time = ds.createDimension('time', None)
lat = ds.createDimension('lat', 64)
lon = ds.createDimension('lon', 128)
times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('slp', 'f4', ('time', 'lat', 'lon',))
value.units = 'Pascals'
value.description = 'difference in mean daily slp from average of 2031 - 2060 to average of 1981 - 2010 '
value.long_name = 'sea level pressure' # add attributes
times.long_name = 'Time'
lats.long_name = 'Latitude'
lons.long_name = 'Longitude'
lats[:] = data_his_lat
lons[:] = data_his_lon

value[:, :, :] = Delta_slp

ds.close()