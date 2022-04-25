# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 16:11:08 2022

@author: Wenbo.Duan
"""

import pandas as pd
import numpy as np

from scipy import stats
from scipy.stats import genpareto
import matplotlib.pyplot as plt

# xls = pd.ExcelFile('C:/tcrm/output/Boston windfield data.xlsx')
xls = pd.ExcelFile('C:/tcrm/output/Wilmington windfield data.xlsx')

data = pd.read_excel(xls, 'Max wind per track')
data2 = pd.read_excel(xls, 'Annual Max')
Max_wind_per_track = data.iloc[1:,1].values
Annual_Max = data2.iloc[1:,0].values
Annual_Max = pd.to_numeric(Annual_Max)
# Max_wind_per_track = np.array(Max_wind_per_track)
# Annual_Max = np.array(Annual_Max)
# print(len(Max_wind_per_track))
print('Number of year:', len(Annual_Max))
years = [200,300,500,700,1000,1700,2000,2500,3000] # Return Period
years = np.array(years)

print(np.sort(Annual_Max))

print('Maximum wind speed recorded in year', np.argmax(Annual_Max, axis=0))
print('2nd largest wind speed recorded in year', Annual_Max.argsort()[-2])
RP = 1 - 1 / years

ASCE_RP= [50, 100, 300, 700, 1700, 3000]
Wilmington = [47.8, 52.8, 60.35, 65.26, 69.29,71.52]
Boston = [41.1277,43.8099, 49.1744, 53.6448,57.2211,59.45]
Miami = [57.66816, 62.13856, 70.18528, 75.54976, 81.36128, 84.49056]

cap = 200 # m/s, to remove outlier
print('Upper Bound Limit:', cap, 'm/s')
threshold = np.linspace(40, 60, 9)
WindSpeed = np.zeros((len(threshold),len(years)))
WindSpeed_gpd = np.zeros((len(threshold),len(years)))

for i, V_limit in enumerate(threshold):
    print('V_limit:', V_limit,'m/s')
    # rec = Annual_Max[(Annual_Max > V_limit) &  (Annual_Max < cap)]
    rec = Max_wind_per_track[(Max_wind_per_track > V_limit) &  (Max_wind_per_track < cap)]
    rec_sorted = np.sort(rec)
    # print(rec_sorted)
    
    ## Gumbel
    stats.gumbel_r.fit(rec_sorted) # sample size still 5000
    SD = np.std(rec_sorted)
    Mean = np.mean(rec_sorted)
    
    # Moment estimator for Gumbel
    # https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
    # 'Parameter Estimation'	
    Beta = np.sqrt(6) * SD / 3.1416 # scale
    Mu = Mean - 0.5772 * Beta # location
    # print(Mu)
    # print(Beta)
    # WindSpeed_1 = -1 * np.log(np.log(years / (years-1))) * Beta + Mu
    
    # WindSpeed_3 = -1 * np.log(-1 * np.log(1 - (1 / years))) * Beta + Mu
    # WindSpeed[i,:] = stats.gumbel_r.ppf(RP, loc=Mu, scale=Beta) # sample size still 5000
    
    ## GPD
    shape, location, scale = genpareto.fit(rec_sorted, floc = V_limit)
    if shape > 0:
        print('Warning: data follows Frechet distribution, Fat tail alert!')
            
    rate = len(rec_sorted) / len(Max_wind_per_track)
    WindSpeed[i,:] = V_limit + (scale / shape) * (np.power(years  * rate, shape) - 1.)
    
    plt.plot(years, WindSpeed[i,:], marker='o',linestyle = '-', label = 'Threshold %.1f m/s' %V_limit)
    
plt.plot(ASCE_RP, Wilmington, marker='o',linestyle = '--', label = 'ASCE')
plt.grid()
plt.xscale("log")
plt.ylabel('Wind speed (m/s)')
plt.xlabel('Return period')
plt.title('Gumbel')
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.show()

plt.scatter(range(len(Max_wind_per_track)), Max_wind_per_track)
# plt.yscale("log")
plt.ylabel('Wind speed (m/s)')
plt.title('Wilmington')
plt.show()
# positive_value = winddata_target_list[winddata_target_list>0]
# plt.scatter(range(len(positive_value)), positive_value)
# plt.show()
