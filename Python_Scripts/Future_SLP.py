# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 14:41:00 2022

@author: Wenbo.Duan
"""
import Utilities.nctools as nctools
import numpy as np
import os
import calendar


path = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\Climate_Model\\CMIP6 climate projections\\SSP3-7.0-Daily-SLP-CanESM5_2030-2040'
file = 'psl_day_CanESM5_ssp370_r1i1p1f1_gn_20300101-20401231_v20190429.nc'
ncobj = nctools.ncLoadFile(os.path.join(path,file))
data_2030_2040 = nctools.ncGetData(ncobj, 'psl')

path = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\Climate_Model\\CMIP6 climate projections\\SSP3-7.0-Daily-SLP-CanESM5_2041-2049'
file = 'psl_day_CanESM5_ssp370_r1i1p1f1_gn_20410101-20491231_v20190429.nc'
ncobj = nctools.ncLoadFile(os.path.join(path,file))
data_2041_2049 = nctools.ncGetData(ncobj, 'psl')


Numberofyear_2030_2040 = int(np.shape(data_2030_2040)[0] / 365)
YearlyData_2030_2040 = np.empty((Numberofyear_2030_2040, 365, 64*128))
for i in range(Numberofyear_2030_2040):
    YearlyData_2030_2040[i,:,:] = data_2030_2040[i*365:365*i+365,:,:].reshape(365, 64*128)


Numberofyear_2041_2049 = int(np.shape(data_2041_2049)[0] / 365)
YearlyData_2041_2049 = np.empty((Numberofyear_2041_2049, 365, 64*128))
for i in range(Numberofyear_2041_2049):
    YearlyData_2041_2049[i,:,:] = data_2041_2049[i*365:365*i+365,:,:].reshape(365, 64*128)
    
YearlyData_2030_2049 = np.concatenate((YearlyData_2030_2040, YearlyData_2041_2049), axis=0)
DailyMean_2041_2049 = np.mean(YearlyData_2030_2049,axis = 0)

DailyMean_2041_2049_grid = DailyMean_2041_2049.reshape(365,64,128)
# reverse the order in axis 1 to match the order in slp.day.ltm.nc
DailyMean_2041_2049_grid = DailyMean_2041_2049_grid[: , ::-1,: ]

