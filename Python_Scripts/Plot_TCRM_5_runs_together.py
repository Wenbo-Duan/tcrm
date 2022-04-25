# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 15:13:00 2022

@author: Wenbo.Duan
"""

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

from Utilities.track import ncReadTrackData, Track
from os.path import join as pjoin
import folium
import numpy as np
import random
from scipy.stats import genpareto,genextreme,linregress
import scipy.stats as stats
import time

start = time.time()
def loadTracks(trackfile):
    """
    Read tracks from a track .nc file and return a list of :class:`Track`
    objects.

    This calls the function `ncReadTrackData` to parse the track .nc
    file.

    :param str trackfile: the track data filename.

    :return: list of :class:`Track` objects.

    """

    tracks = ncReadTrackData(trackfile)
    return tracks



# Folder_list = {'Boston_10by10_5000yr_Current_run1',
#                 'Boston_10by10_5000yr_Current_run2',
#                 'Boston_10by10_5000yr_Current_run3',
#                 'Boston_10by10_5000yr_Current_run4',
#                 'Boston_10by10_5000yr_Current_run5'}
Folder_list = {'Boston_10by10_5000yr_2060_CanESM5_SSP85',
                'Boston_10by10_5000yr_2060_CESM2_SSP85',
                'Boston_10by10_5000yr_2060_CSM2_SSP85',
                'Boston_10by10_5000yr_2060_Earth3Veg_SSP85',
                'Boston_10by10_5000yr_2060_HadGEM3LL_SSP85'}

ASCE_RP = [100, 300, 700, 1700, 3000]

ASCE_WindSpeed = {'9020': [98,	109,	120,	131,	136],
'11922':	[119,	136,	147,	156,	161], #Wrightsville Beach
'11992': [117,	132,	143,	150,	155],
'12054': [121,	137,	148,	154,	158],
'12068': [115,	133,	149,	159,	166],
'12527':	[102,	115,	126,	136,	142], 
'8764':	[98,	110,	120,	128,	133],# LOGAN
'8570':	[122,	138,	150,	160,	170],
'8567':	[138,	156,	167,	179,	186]} # Miami airport

LocationID = 8764

years = [100,200,300,500,700,1000,1700,2000,2500,3000] # Return Period
years = np.array(years)

Wind_Hazard = np.zeros((len(years),len(Folder_list)))
numsim = 5000

for i_run,folder_name in enumerate(Folder_list):
    print(folder_name)
    # database_path = 'G:/group_data/Advanced Technology & Research/Climate Analysis/079007-70 IiA Extreme Winds/tcrm/output/Miami_10by10_5000yr_1980_run2/hazard.db'
    # database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\'
    database_path = 'C:\\Users\\wenbo.duan\\OneDrive - Arup\\03_Project\\IiA\\079007-70_Climate_change_Extreme_wind\\TCRM note\\out\\'
    file = 'hazard.db'
    
    tblHazard = pd.read_pickle(database_path + folder_name + '\\' + 'tblHazard.pkl')
    
    # tblWindSpeed = pd.read_pickle(database_path + folder_name + '\\' + 'tblWindSpeed.pkl')
    
    
    Hazard_at_Location = tblHazard[tblHazard['locId'] == LocationID]
    # df_temp = Hazard_at_Location.loc[:,['wspd']]* 2.23694 
    # Hazard_at_Location[:,['wspd']] = df_temp
    Wind_Hazard[:,i_run] = Hazard_at_Location['wspd']* 2.23694 # convert to mph
    
    # Wind_Hazard[:,i_run] = Hazard_at_Location['wspd']



# # ax = Hazard_at_Location.plot(x ='returnPeriod', y='wspd', kind = 'line')
# plt.boxplot(Wind_Hazard[[0,2,4,6,9],:].T,positions=years[[0,2,4,6,9]])
# # ax =  plt.gca()
# plt.plot(years,np.median(Wind_Hazard,axis=1),label = 'TCRM', color='black')
yerr = [ abs(np.min(Wind_Hazard,axis=1) - np.median(Wind_Hazard,axis=1)),np.max(Wind_Hazard,axis=1) - np.median(Wind_Hazard,axis=1)]
yerr = np.column_stack(yerr)
plt.errorbar(years[[0,2,4,6,9]], np.median(Wind_Hazard[[0,2,4,6,9],:],axis=1), yerr=yerr[[0,2,4,6,9],:].T,elinewidth=2,capsize=5, capthick = 2)
plt.plot(ASCE_RP,  [number for number in ASCE_WindSpeed[str(LocationID)] ] ,linestyle = '--', label = 'ASCE', color='red')

# ax.plot(years, Wind_Hazard,linestyle = '--', label = 'ASCE', color='red')
ax =  plt.gca()
ax.set_xscale('log')
ax.set_ylabel('Wind Speed (mph)')
ax.set_xlabel('Return Period (year)')
ax.set_xticks([100,300,1000,3000])
ax.set_xticklabels(['100','300','1000','3000'])
ax.set_xlim(100,3000)
ax.set_ylim(0,225)
ax.grid(visible=True, which='major', axis='both',linestyle='--')
ax.legend(['ASCE', 'TCRM'],bbox_to_anchor=(1.05, 1.0), loc='upper left')
# # plt.set_xticks([200,300,500,700,1000,1700,2000,2500,3000])
plt.savefig(database_path + '\\' + 'Boston 2060 ssp8.5 TCRM_vs_ASCE ID%i.png' %LocationID,bbox_inches='tight')
plt.show()


plt.plot(ASCE_RP,  [number for number in ASCE_WindSpeed[str(LocationID)] ] ,linestyle = '--', label = 'ASCE', color='red')
plt.plot(years[[0,2,4,6,9]],  Wind_Hazard[[0,2,4,6,9],2])
plt.plot(years[[0,2,4,6,9]],  Wind_Hazard[[0,2,4,6,9],3])
plt.plot(years[[0,2,4,6,9]],  Wind_Hazard[[0,2,4,6,9],1])
plt.plot(years[[0,2,4,6,9]],  Wind_Hazard[[0,2,4,6,9],0])
plt.plot(years[[0,2,4,6,9]],  Wind_Hazard[[0,2,4,6,9],4])
ax =  plt.gca()
ax.set_xscale('log')
ax.set_ylabel('Wind Speed (mph)')
ax.set_xlabel('Return Period (year)')
ax.set_xticks([100,300,1000,3000])
ax.set_xticklabels(['100','300','1000','3000'])
ax.set_xlim(100,3000)
ax.set_ylim(0,225)
ax.grid(visible=True, which='major', axis='both',linestyle='--')
ax.legend(['ASCE', 'CESM2','HadGEM3-GC31-LL','CanESM5','Earth3-Veg-LL','CSM2'],bbox_to_anchor=(1.05, 1.0), loc='upper left')
# ax.legend(['ASCE', 'Earth3-Veg-LL', 'CanESM5','CESM2','HadGEM3-GC31-LL','CSM2'],bbox_to_anchor=(1.05, 1.0), loc='upper left')

# # plt.set_xticks([200,300,500,700,1000,1700,2000,2500,3000])
plt.savefig(database_path + '\\' + 'Boston 2060 ssp8.5 GCMs_vs_ASCE ID%i.png' %LocationID,bbox_inches='tight')
plt.show()
