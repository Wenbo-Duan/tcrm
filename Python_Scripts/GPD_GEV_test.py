# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 08:32:29 2022

@author: Wenbo.Duan
"""

import pandas as pd
import numpy as np
from scipy.stats import genpareto, scoreatpercentile
import matplotlib.pyplot as plt
import Utilities.lmomentFit as lmom 


def gpdReturnLevel(intervals, mu, shape, scale, rate, npyr=365.25):
    """
    Calculate return levels for specified intervals for a distribution with
    the given threshold, scale and shape parameters.

    :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
              to evaluate return levels for.
    :param float mu: Threshold parameter (also called location).
    :param float shape: Shape parameter.
    :param float scale: Scale parameter.
    :param float rate: Rate of exceedances (i.e. number of observations greater
                       than `mu`, divided by total number of observations).
    :param float npyr: Number of observations per year.

    :returns: return levels for the specified recurrence intervals.

    """

    rp = mu + (scale / shape) * (np.power(intervals * npyr * rate, shape) - 1.)
    return rp


# xls = pd.ExcelFile('C:/tcrm/output/Boston windfield data.xlsx')
xls = pd.ExcelFile('C:/tcrm/output/Wilmington windfield data.xlsx')

data = pd.read_excel(xls, 'Max wind per track')
data2 = pd.read_excel(xls, 'Annual Max')
Max_wind_per_track = data.iloc[1:,1].values
Annual_Max = data2.iloc[1:,0].values
Annual_Max = pd.to_numeric(Annual_Max)
Max_wind_per_track = np.array(Max_wind_per_track)
# Annual_Max = np.array(Annual_Max)
print(len(Max_wind_per_track))
print(len(Annual_Max))

numsim = 5000 #
years = [200,300,500,700,1000,1700,2000,2500,3000] # Return Period
years = np.array(years)
# years = [200,300,500,700,1000,1700,2000,2500,3000] # Return Period


ASCE_RP= [50, 100, 300, 700, 1700, 3000]
Welmington = [47.8, 52.8, 60.35, 65.26, 69.29,71.52]
Boston = [41.1277,43.8099, 49.1744, 53.6448,57.2211,59.45]

#### GEV

l1, l2, l3 = lmom.samlmu(Annual_Max, 3) # find 3 l-moments
# t3 = L-skewness (Hosking 1990)
t3 = l3 / l2 
xmom = [l1, l2, t3]
# Calculates the parameters of the distribution given its
# L-moments. GEV distribution calculated.
loc, scal, shp = np.array(lmom.pelgev(xmom))
w = np.ones(len(years)) # Create empty return period array
yrspersim = 1

for i, t in enumerate(years): 
# if no valid fit was found, then there are no return period wind speeds
    if shp <= 0.:
        w[i] = 9999
# if a valid fit was found...
    else:
# Calculate wind speed for each return period
        w[i] = (np.transpose(loc + (scal / shp) * (1. - np.power(-1. * np.log(1. - (yrspersim / t)), shp))))

plt.plot(years, w, marker='o',linestyle = '-', label = 'GEV')
plt.plot(ASCE_RP, Welmington, marker='o',linestyle = '--', label = 'ASCE')
plt.grid()
plt.ylabel('Wind speed (m/s)')
plt.xlabel('Return period')
plt.title('GEV')
plt.legend()
plt.show()

#### GPD fitting
Max_wind_per_track=np.sort(Max_wind_per_track)
recs = Max_wind_per_track[Max_wind_per_track > 0]
# mu = scoreatpercentile(data, threshold)
# mu = scoreatpercentile(recs, threshold) # calculate threshold using obersavations with speed values - WD
mu = recs[-50] # calculate threshold using obersavations with speed values - WD


Rp = np.ones(len(years))

# Fill each day that a cyclone isn't recorded with zero so we get 
# the correct rate for the return periods
datafilled = np.zeros(int(numsim * 365.25))
datafilled[-len(Max_wind_per_track):] = Max_wind_per_track


rate = float(len(datafilled[datafilled > mu])) / float(len(datafilled))

# try to fit GPD distribution 

try:
    shape, location, scale = genpareto.fit(datafilled[datafilled > mu],
                                            floc = mu)
except:
    print("fitting error")

Rpeval = gpdReturnLevel(years, mu, shape, scale, rate)
    
plt.plot(years, Rpeval, marker='o',linestyle = '-', label = 'GPD')
plt.plot(ASCE_RP, Welmington, marker='o',linestyle = '--', label = 'ASCE')
plt.grid()
plt.ylabel('Wind speed (m/s)')
plt.xlabel('Return period')
plt.title('GPD')
plt.legend()
plt.show()