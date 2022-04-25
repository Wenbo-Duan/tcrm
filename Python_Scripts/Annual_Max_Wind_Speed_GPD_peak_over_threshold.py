# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 16:11:08 2022

@author: Wenbo.Duan
"""

import pandas as pd
import numpy as np

from scipy.stats import genpareto,genextreme,linregress
import scipy.stats as stats
import matplotlib.pyplot as plt

# read wind speed data
# xls = pd.ExcelFile('C:/tcrm/output/Boston windfield data.xlsx')
xls = pd.ExcelFile('Wilmington wind speed data run 3.xlsx')

# data2 = pd.read_excel(xls, 'Annual Max')
# Annual_Max = data2.iloc[1:,0].values

data2 = pd.read_excel(xls, 'Max wind per track')
Annual_Max = data2.iloc[1:,1].values

Annual_Max = pd.to_numeric(Annual_Max)
print('Number of year:', len(Annual_Max))
print(np.sort(Annual_Max))
print('Maximum wind speed recorded in year', np.argmax(Annual_Max, axis=0))
print('2nd largest wind speed recorded in year', Annual_Max.argsort()[-2])

years = [100,200,300,500,700,1000,1700,2000,2500,3000] # Return Period
years = np.array(years)

P = 1 - 1 / years

# code value in m/s
ASCE_RP= [50, 100, 300, 700, 1700, 3000]
Wilmington = [47.8, 52.8, 60.35, 65.26, 69.29,71.52]
# Boston = [41.1277,43.8099, 49.1744, 53.6448,57.2211,59.45]
# Miami = [57.66816, 62.13856, 70.18528, 75.54976, 81.36128, 84.49056]

cap = 100 # m/s, upper bound
print('Upper Bound Limit:', cap, 'm/s')
threshold = np.linspace(40, 60, 9) # lower bound
WindSpeed = np.zeros((len(threshold),len(years)))

Dist = 'GPD' # 'GEV' 'GPD' 'Gumbel'
# Calculate empirical mean excess function
# ref: file:///C:/Users/wenbo.duan/OneDrive%20-%20Arup/03_Project/IiA/079007-70_Climate_change_Extreme_wind/TCRM%20note/Python_Scripts/POT%20approach/Chapter4.pdf
# page 6
emef = np.zeros(len(Annual_Max))
for n, x_n in enumerate(Annual_Max):
    diff = Annual_Max - Annual_Max[n]
    diff = diff[diff > 0]
    emef[n] = np.mean(diff)
    
# plt.scatter(Annual_Max, emef)
# # plt.yscale("log")
# plt.grid()
# plt.xlabel('Threholds (m/s)')
# plt.ylabel('Empirical Mean Excess Function')
# plt.show()


for i, V_limit in enumerate(threshold):
    print('V_limit:', V_limit,'m/s')
    rec = Annual_Max[(Annual_Max > V_limit) & (Annual_Max < cap)] # truncate data
    rec_sorted = np.sort(rec)
    
    # # Calculate empirical mean excess function
    # # ref: file:///C:/Users/wenbo.duan/OneDrive%20-%20Arup/03_Project/IiA/079007-70_Climate_change_Extreme_wind/TCRM%20note/Python_Scripts/POT%20approach/Chapter4.pdf
    # # page 6
    # emef[i] = np.mean(rec_sorted - V_limit)
    
    if Dist == 'Gumbel':
        ### Gumbel
        rank = np.linspace(1,len(rec_sorted),len(rec_sorted))
        inverse_rank =  np.sort(len(Annual_Max) + 1 - rank)
        CDF = 1 - inverse_rank / (len(Annual_Max) + 1)
        Var = -1 * np.log(-1 * np.log(1 - CDF))
        slope, intercept, r, p, se = linregress(Var,rec_sorted)
        WindSpeed[i,:] = slope * (-1 * np.log(-1 * np.log(1 - 1 / years))) + intercept
        # np.sum((rec_sorted - np.mean(rec_sorted)) * (Var - np.mean(Var))) / np.sum((rec_sorted - np.mean(rec_sorted)) ** 2)
        # plt.scatter(rec_sorted, Var)
        # plt.show()
    
    # # stats.gumbel_r.fit(rec_sorted) # sample size still 5000
    # SD = np.std(rec_sorted)
    # Mean = np.mean(rec_sorted)
    # # Moment estimator for Gumbel
    # # https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
    # # 'Parameter Estimation'	
    # Beta = np.sqrt(6) * SD / 3.1416 # scale
    # Mu = Mean - 0.5772 * Beta # location
    # # print(Mu)
    # # print(Beta)
    # # WindSpeed_1 = -1 * np.log(np.log(years / (years-1))) * Beta + Mu
    
    # # WindSpeed_3 = -1 * np.log(-1 * np.log(1 - (1 / years))) * Beta + Mu
    # WindSpeed[i,:] = stats.gumbel_r.ppf(P, loc=Mu, scale=Beta) # sample size still 5000
    
    if Dist == 'GPD':
        #### GPD generalized Pareto distribution
        shape, location, scale = genpareto.fit(rec_sorted, floc = V_limit)
        
        if shape > 0:
            print('Warning: Positive shape factor, Fat tail alert!')
        if shape > -0.01 and shape < 0:
            print('Warning: Shape factor close to 0')
            
        rate = len(rec_sorted) / len(Annual_Max)
        WindSpeed[i,:] = V_limit + (scale / shape) * (np.power(years  * rate, shape) - 1.)
    
    if Dist == 'GEV':
        #### GEV Generalized extreme value distribution
        shape, location, scale = genextreme.fit(rec_sorted, floc = V_limit) 
        if shape > 0:
            print('Warning: data fitted to Frechet distribution, Fat tail alert!')
        for ii, t in enumerate(years): 
            WindSpeed[i,ii] = (np.transpose(location + (scale / shape) * (1. - np.power(-1. * np.log(1. - (1 / t)), shape))))
       
     
    plt.plot(years, WindSpeed[i,:], marker='o',linestyle = '-', label = 'Threshold %.1f m/s' %V_limit)
    
plt.plot(ASCE_RP, Wilmington, marker='o',linestyle = '--', label = 'ASCE')
plt.grid()
plt.xscale("log")
plt.ylabel('Wind speed (m/s)')
plt.xlabel('Return period')
if Dist == 'Gumbel':
    plt.title('Gumbel distribution fit, Cap = %.1f m/s' %cap)
if Dist == 'GPD':
    plt.title('Generalized Pareto distribution fit, Cap = %.1f m/s' %cap)
if Dist == 'GEV':
    plt.title('Generalized extreme value distribution fit, Cap = %.1f m/s' %cap)
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.show()


# plt.scatter(range(len(Annual_Max)), Annual_Max)
# # plt.yscale("log")
# plt.ylabel('Wind speed (m/s)')
# plt.title('Wilmington')
# plt.show()
# positive_value = winddata_target_list[winddata_target_list>0]
# plt.scatter(range(len(positive_value)), positive_value)
# plt.show()
