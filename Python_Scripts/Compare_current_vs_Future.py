# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:59:02 2022

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


database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980'

database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_run2'

database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_run3'

database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP'

database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP_CESM2'

database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP_GFDL-ESM4'
database_path = ['G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980',

'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_run2',

 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_run3',

 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP',

'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP_CESM2',
 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Wilmington_10by10_5000yr_1980_future_MSLP_GFDL-ESM4']
ASCE_RP = [100, 300, 700, 1700, 3000]

ASCE_WindSpeed = {'9020': [98,	109,	120,	131,	136],
'11922':	[119,	136,	147,	156,	161],
'11992': [117,	132,	143,	150,	155],
'12054': [121,	137,	148,	154,	158],
'12068': [115,	133,	149,	159,	166],
'12527':	[102,	115,	126,	136,	142]}
LocationID = 12068
plt.plot(ASCE_RP,  [number / 2.23694 for number in ASCE_WindSpeed[str(LocationID)] ] , marker='o',linestyle = '--', label = 'ASCE', color='red')

for i,path in enumerate(database_path):
    tblHazard = pd.read_pickle(path + '\\' + 'tblHazard.pkl')
    
    # tblWindSpeed = pd.read_pickle(path + '\\' + 'tblWindSpeed.pkl')
    
    
    Hazard_at_Location = tblHazard[tblHazard['locId'] == LocationID]
    
    Hazard_at_Location.plot(x ='returnPeriod', y='wspd', kind = 'line', marker='o')
    
    plt.xscale("log")
    plt.ylabel('Wind speed (m/s)')
    plt.xlabel('Return period')
    plt.legend(['TCRM', 'ASCE'],bbox_to_anchor=(1.05, 1.0), loc='upper left')
    # plt.set_xticks([200,300,500,700,1000,1700,2000,2500,3000])
    # plt.savefig(database_path + '\\' + 'Present vs Future ID%i.png' %LocationID,bbox_inches='tight')
    plt.show()
