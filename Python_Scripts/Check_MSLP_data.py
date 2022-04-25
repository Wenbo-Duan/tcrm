# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 15:06:02 2022

Check if data from slp.day.ltm.nc if just an average of value from each year
of 1981 - 2010
@author: Wenbo.Duan
"""
import Utilities.nctools as nctools
import numpy as np
import os
import calendar

# get slp.day.ltm.nc data
ncobj = nctools.ncLoadFile('C:/tcrm/MSLP/slp.day.ltm.nc')
data = nctools.ncGetData(ncobj, 'slp')

data2 = data.reshape(1, 365*73*144).T

lat_grid = np.sort(np.arange(-90, 92.5, 2.5))[::-1]
lon_grid = np.arange(0,360,2.5)
day_grid = np.linspace(1,365,365)

# get yearly data

Yearly_path = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\Mean Sea Level Pressure\\Yearly Data'

filelist = os.listdir(Yearly_path)
YearlyData = len(filelist)

YearlyData = np.empty((len(filelist), 365, 10512))

for n_file,file in enumerate(filelist):
    
    test = nctools.ncLoadFile(os.path.join(Yearly_path,file))

    dataj = nctools.ncGetData(test, 'slp')
    if calendar.isleap(int(file.split('.')[1])):
       dataj = np.delete(dataj,59,0)
    YearlyData[n_file,:,:] = dataj.reshape(365, 73*144)

DailyMean = np.mean(YearlyData,axis = 0)

DailyMean_grid = DailyMean.reshape(365,73,144)

diff = data-DailyMean_grid
print(np.max(np.absolute(diff)))
# df.to_csv('slp.day.ltm.csv')
# slpunits = getattr(ncobj.variables[var], 'units')