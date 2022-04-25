# -*- coding: utf-8 -*-
"""
Created on Feb 10 2022

@author: Jens.Bauer
"""

from pickle import FALSE

import Utilities.nctools as nctools
import numpy as np
import os
import calendar
from scipy.interpolate import griddata
import netCDF4 as nc

path = 'C:\Climate_data\CESM2 (USA)\SSP3_7.0_Daily_SLP'
file_list = ['psl_day_CESM2_ssp370_r4i1p1f1_gn_20310101-20401231_v20200528.nc', 
             'psl_day_CESM2_ssp370_r4i1p1f1_gn_20410101-20501231_v20200528.nc',
             'psl_day_CESM2_ssp370_r4i1p1f1_gn_20510101-20601231_v20200528.nc']

inputshape = (192, 288) # shape of data in projection files
YearlyData_2031_2060 = np.ones((30, 365,inputshape[0]*inputshape[1]))*-1


processed_years = 0
for file in file_list:
    ncobj = nctools.ncLoadFile(os.path.join(path,file))
    datainfile = nctools.ncGetData(ncobj, 'psl')
    data_lon = nctools.ncGetData(ncobj, 'lon')
    data_lat = nctools.ncGetData(ncobj, 'lat')
    print(datainfile.shape)

    Numberofyearinfile = int(np.shape(datainfile)[0] / 365)
    for i in range(Numberofyearinfile):
        YearlyData_2031_2060[processed_years,:,:] = datainfile[i*365:365*i+365,:,:].reshape(365, datainfile.shape[1]*datainfile.shape[2])
        processed_years += 1
    
# check if all of the array was filled with data
if (-1 in YearlyData_2031_2060): print("projection data NOT correct")

DailyMean_2031_2060 = np.mean(YearlyData_2031_2060,axis = 0)
DailyMean_2031_2060_grid = DailyMean_2031_2060.reshape(365,inputshape[0], inputshape[1])
DailyMean_2031_2060_grid = DailyMean_2031_2060_grid[: , ::-1,: ] # lat in decending order
data_lat = data_lat[::-1 ]# lat in decending order


#load historical data from CanESM5
path_his = 'C:\Climate_data\CESM2 (USA)\historical_data_Daily_SLP'
file_his_list = ['psl_day_CESM2_historical_r1i1p1f1_gn_18500101-20150101_v20190308.nc']

print("Import his data")

# load the three historical datasets
ncobj_his_1 = nctools.ncLoadFile(os.path.join(path_his,file_his_list[0]))
data_his_1 = nctools.ncGetData(ncobj_his_1, 'psl')
data_his_lon= nctools.ncGetData(ncobj_his_1, 'lon')
data_his_lat = nctools.ncGetData(ncobj_his_1, 'lat')


his_1981_2010 = np.ones((30*365, inputshape[0],inputshape[1]))*-1
# insert data 
his_1981_2010[:,:,:] = data_his_1[(1981 - 1850) * 365 : (2011 - 1850) * 365,:,:]


if (-1 in his_1981_2010): print("his data NOT correct")


Numberofyear_1981_2010 = int(np.shape(his_1981_2010)[0] / 365)
YearlyData_1981_2010 = np.empty((Numberofyear_1981_2010, 365, data_his_1.shape[1]*data_his_1.shape[2]))
for i in range(Numberofyear_1981_2010):
    YearlyData_1981_2010[i,:,:] = his_1981_2010[i*365:365*i+365,:,:].reshape(365, data_his_1.shape[1]*data_his_1.shape[2])

DailyMean_1981_2010 = np.mean(YearlyData_1981_2010,axis = 0)
DailyMean_1981_2010_grid = DailyMean_1981_2010.reshape(365,data_his_1.shape[1],data_his_1.shape[2])
DailyMean_1981_2010_grid = DailyMean_1981_2010_grid[: , ::-1,: ] # lat in decending order
data_his_lat = data_his_lat[::-1 ]# lat in decending order

Delta_slp = DailyMean_2031_2060_grid - DailyMean_1981_2010_grid


# get slp.day.ltm.nc data
ncobj = nctools.ncLoadFile('C:/tcrm/MSLP/slp.day.ltm.nc')
data = nctools.ncGetData(ncobj, 'slp')

data2 = data.reshape(1, 365*73*144).T

lat_grid = np.sort(np.arange(-90, 92.5, 2.5))[::-1]
lon_grid = np.arange(0,360,2.5)
day_grid = np.linspace(1,365,365)

grid1, grid2 = np.meshgrid(data_lat, data_lon, indexing='ij')

grid1_reshape = grid1.reshape(data_his_1.shape[1]*data_his_1.shape[2])
grid2_reshape = grid2.reshape(data_his_1.shape[1]*data_his_1.shape[2])
intp_point = np.vstack((grid1_reshape,grid2_reshape))
outgrid1, outgrid2 = np.meshgrid(lat_grid, lon_grid, indexing='ij')

Delta_slp_intp = np.empty((365,73,144))

for i in range(365):
    if i % 10 == 0: print(i)

    test_data = Delta_slp[i,:,:]

    intp_value = test_data.reshape(data_his_1.shape[1]*data_his_1.shape[2])

    outgrid1, outgrid2 = np.meshgrid(lat_grid, lon_grid, indexing='ij')

    Delta_slp_intp[i,:,:] = griddata(intp_point.T, intp_value, (outgrid1, outgrid2), method='linear')
    
    
Corrected_future_mslp = data + np.nan_to_num(Delta_slp_intp)


# create NC file
fn = 'C:/tcrm/MSLP/Corrected.CESM2.mslp.day.2031-2060.1981-2010.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
time = ds.createDimension('time', None)
lat = ds.createDimension('lat', 73)
lon = ds.createDimension('lon', 144)
times = ds.createVariable('time', 'f4', ('time',))
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('slp', 'f4', ('time', 'lat', 'lon',))
value.units = 'Pascals'
value.description = 'Porjected mean daily slp from average of 2031 - 2060, = delta value + historical average of 1981 - 2010 '
value.long_name = 'sea level pressure' # add attributes
times.long_name = 'Time'
lats.long_name = 'Latitude'
lons.long_name = 'Longitude'
lats[:] = lat_grid
lons[:] = lon_grid

value[:, :, :] = Corrected_future_mslp

ds.close()
# #=======================
# #=== Interpolate the data for plotting
# #=======================
# idx1 = np.logical_and(p[:,0] >= xmin, p[:,0] <= xmax)
# idx2 = np.logical_and(p[:,1] >= ymin, p[:,1] <= ymax)
# idxrange = np.logical_and(idx1,idx2)



# data_intrp = data[idxrange,:]
# p_intrp = p[idxrange,:]
# p_intrp = p_intrp[:,[0,1]]



# if variable_to_load == 'U':
# uCx = griddata(p_intrp, data_intrp[:,0], (xi, yi), method='linear')
# uCy = griddata(p_intrp, data_intrp[:,1], (xi, yi), method='linear')
# VCi = np.sqrt(uCx**2 + uCy**2)
# elif variable_to_load == 'k':
# VCi = griddata(p_intrp, data_intrp[:], (xi, yi), method='linear')
# VCi = VCi[:,:,0]
# data_grid1, data_grid2 = np.meshgrid(data_lat, data_lon, indexing='ij')