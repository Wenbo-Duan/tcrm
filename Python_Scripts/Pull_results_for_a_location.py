# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 12:25:35 2022

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



# database_path = 'G:/group_data/Advanced Technology & Research/Climate Analysis/079007-70 IiA Extreme Winds/tcrm/output/Wilmington_10by10_5000yr_1980_run2/hazard.db'
database_path = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Boston_10by10_5000yr_Current_run1'
file = 'hazard.db'

# trackpath= 'C:/tcrm/output/Boston_10by10_5000yr_1980/tracks/'


# database_path = 'C:\\tcrm\\input\\ibtracs.since1980.list.v04r00.csv'

# con = sqlite3.connect(database_path +'\\'+file)
# cur = con.cursor()
# cur.execute('SELECT * FROM tblHazard')
# tblHazard = pd.DataFrame(cur.fetchall())
# tblHazard.columns = [x[0] for x in cur.description]
# tblHazard.to_pickle(database_path + '\\' + 'tblHazard.pkl')

tblHazard = pd.read_pickle(database_path + '\\' + 'tblHazard.pkl')

# cur.execute('SELECT * FROM tblTracks')
# tblTracks = pd.DataFrame(cur.fetchall())
# tblTracks.columns = [x[0] for x in cur.description]
# tblTracks.to_pickle(database_path + '\\' + 'tblTracks.pkl')

# cur.execute('SELECT * FROM tblWindSpeed')
# tblWindSpeed = pd.DataFrame(cur.fetchall())
# tblWindSpeed.columns = [x[0] for x in cur.description]
# tblWindSpeed.to_pickle(database_path + '\\' + 'tblWindSpeed.pkl')

tblWindSpeed = pd.read_pickle(database_path + '\\' + 'tblWindSpeed.pkl')

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

Hazard_at_Location = tblHazard[tblHazard['locId'] == LocationID]
# df_temp = Hazard_at_Location.loc[:,['wspd']]* 2.23694 
# Hazard_at_Location[:,['wspd']] = df_temp
Hazard_at_Location['wspd'] = Hazard_at_Location['wspd']* 2.23694 # convert to mph
ax = Hazard_at_Location.plot(x ='returnPeriod', y='wspd', kind = 'line')

ax.plot(ASCE_RP,  [number for number in ASCE_WindSpeed[str(LocationID)] ] ,linestyle = '--', label = 'ASCE', color='red')
ax.set_xscale('log')
ax.set_ylabel('Wind speed (mph)')
ax.set_xlabel('Return Period (year)')
ax.set_xticks([100,300,1000,3000])
ax.set_xticklabels(['100','300','1000','3000'])
ax.set_xlim(100,3000)
ax.set_ylim(0,225)
ax.grid(visible=True, which='major', axis='both',linestyle='--')
ax.legend(['TCRM', 'ASCE'],bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.set_xticks([200,300,500,700,1000,1700,2000,2500,3000])
plt.savefig(database_path + '\\' + 'TCRM_vs_ASCE ID%i.png' %LocationID,bbox_inches='tight')
plt.show()


WindSpeed_at_Location = tblWindSpeed[tblWindSpeed['locId'] == LocationID]

WindSpeed_value = WindSpeed_at_Location['wspd'].values
WindSpeed_value = pd.to_numeric(WindSpeed_value)
WindSpeed_value = WindSpeed_value * 2.23694 # convert to mph
years = [100,200,300,500,700,1000,1700,2000,2500,3000] # Return Period
years = np.array(years)

P = 1 - 1 / years
cap = 450 # mph, upper bound
print('Upper Bound Limit:', cap, 'mph')
threshold = np.linspace(67.1, 111.85, 9) # lower bound 
WindSpeed_RP = np.zeros((len(threshold),len(years)))
numsim = 5000
Dist = 'GPD' # 'GEV' 'GPD' 'Gumbel'
# Calculate empirical mean excess function
# ref: file:///C:/Users/wenbo.duan/OneDrive%20-%20Arup/03_Project/IiA/079007-70_Climate_change_Extreme_wind/TCRM%20note/Python_Scripts/POT%20approach/Chapter4.pdf
# page 6
emef = np.zeros(len(WindSpeed_value))
for n, x_n in enumerate(WindSpeed_value):
    diff = WindSpeed_value - WindSpeed_value[n]
    diff = diff[diff > 0]
    emef[n] = np.mean(diff)
    
# plt.scatter(WindSpeed_value, emef)
# # plt.yscale("log")
# plt.grid()
# plt.xlabel('Threholds (m/s)')
# plt.ylabel('Empirical Mean Excess Function')
# plt.show()


for i, V_limit in enumerate(threshold):
    print('V_limit:', V_limit,'mph')
    rec = WindSpeed_value[(WindSpeed_value > V_limit) & (WindSpeed_value < cap)] # truncate data
    rec_sorted = np.sort(rec)
    
    # # Calculate empirical mean excess function
    # # ref: file:///C:/Users/wenbo.duan/OneDrive%20-%20Arup/03_Project/IiA/079007-70_Climate_change_Extreme_wind/TCRM%20note/Python_Scripts/POT%20approach/Chapter4.pdf
    # # page 6
    # emef[i] = np.mean(rec_sorted - V_limit)
    
    if Dist == 'Gumbel':
        ### Gumbel
        rank = np.linspace(1,len(rec_sorted),len(rec_sorted))
        inverse_rank =  np.sort(len(WindSpeed_value) + 1 - rank)
        CDF = 1 - inverse_rank / (len(WindSpeed_value) + 1)
        Var = -1 * np.log(-1 * np.log(1 - CDF))
        slope, intercept, r, p, se = linregress(Var,rec_sorted)
        WindSpeed_RP[i,:] = slope * (-1 * np.log(-1 * np.log(1 - 1 / years))) + intercept
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
        datafilled = np.zeros(int(numsim * 365.25))
        datafilled[-len(WindSpeed_value):] = WindSpeed_value
        # log.debug("The length of the filled data is {0}".format(len(datafilled)))
    
        rate = float(len(datafilled[datafilled > V_limit])) / float(len(datafilled))
        # print(rate)    
        # rate2 = len(rec_sorted) / len(WindSpeed_value)
        WindSpeed_RP[i,:] = V_limit + (scale / shape) * (np.power(years  * rate * 365.25, shape) - 1.)
        
    if Dist == 'GEV':
        #### GEV Generalized extreme value distribution
        shape, location, scale = genextreme.fit(rec_sorted, floc = V_limit) 
        if shape > 0:
            print('Warning: data fitted to Frechet distribution, Fat tail alert!')
        for ii, t in enumerate(years): 
            WindSpeed_RP[i,ii] = (np.transpose(location + (scale / shape) * (1. - np.power(-1. * np.log(1. - (1 / t)), shape))))
       
     
    plt.plot(years, WindSpeed_RP[i,:],linestyle = '-', label = 'Threshold %.1f mph' %V_limit)
ax2 =  plt.gca()
ax2.plot(ASCE_RP,  [number for number in ASCE_WindSpeed[str(LocationID)] ],linestyle = '--', label = 'ASCE', color='red')
# plt.grid()
# plt.xscale("log")
# plt.ylabel('Wind speed (mph)')
# plt.xlabel('Return period')
# plt.ylim([0, 220])
# plt.xlim([100, 3000])

ax2.set_xscale('log')
ax2.set_ylabel('Wind speed (mph)')
ax2.set_xlabel('Return Period (year)')
ax2.set_xticks([100,300,1000,3000])
ax2.set_xticklabels(['100','300','1000','3000'])
ax2.set_xlim(100,3000)
ax2.set_ylim(0,225)
ax2.grid(visible=True, which='major', axis='both',linestyle='--')

if Dist == 'Gumbel':
    ax2.set_title('Gumbel distribution fit, Cap = %.1f mph' %cap)
if Dist == 'GPD':
    ax2.set_title('Generalized Pareto distribution fit, Cap = %.1f mph' %cap)
if Dist == 'GEV':
    ax2.set_title('Generalized extreme value distribution fit, Cap = %.1f mph' %cap)
ax2.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.savefig(database_path + '\\' + 'POT fit TCRM ID%i.png' %LocationID,bbox_inches='tight')
plt.show()

# find and plot tracks produce largest 10 wind speed
Top10 = WindSpeed_at_Location.sort_values(by=['wspd'],ascending=False).head(10)

my_map = folium.Map(location=[40, -71],zoom_start = 5)

for i_track, track_info in enumerate(Top10['eventId']):
    print(track_info)
    i, year = track_info.split('-')
    trackname = 'tracks.' + year + '.nc' # Boston 2nd largest wind speed
    
    trackfile = database_path + '\\tracks\\' + trackname

    trackiter = loadTracks(trackfile)

    Track_Data = trackiter[int(i)]

    # Track.Longitude
    
    # Coordinate= [Track.Latitude  Track.Longitude -360]
    Coordinate = np.vstack((Track_Data.Latitude , Track_Data.Longitude -360)).T
    # Coordinate = np.vstack((Track_Data.Latitude[::3] , Track_Data.Longitude[::3] -360)).T
    
    hexadecimal = ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
    folium.PolyLine(Coordinate,  popup=track_info, color=hexadecimal, weight=2.5, opacity=1).add_to(my_map)
    
# domain Boston
Track_generation_domain=folium.PolyLine(locations=[[8, 263-360],[8, 344-360],[47, 344-360],[47, 263-360],[8, 263-360]],weight=2)
my_map.add_child(Track_generation_domain)

Wind_domain=folium.PolyLine(locations=[[36.5, -76],[36.5, -66],[46.5, -66],[46.5,-76],[36.5, -76]],weight=2, color = 'red')
my_map.add_child(Wind_domain)

# # domain Wilmington
# Track_generation_domain=folium.PolyLine(locations=[[8, 262-360],[8, 344-360],[39, 344-360],[39, 262-360],[8, 262-360]],weight=2)
# my_map.add_child(Track_generation_domain)


# Wind_domain=folium.PolyLine(locations=[[29, -83],[29, -73],[39, -73],[39,-83],[29, -83]],weight=2, color = 'red')
# my_map.add_child(Wind_domain)

my_map.save(database_path + '\\' + "Top 50 Tracks 1hr ID%i.html" %LocationID)

end = time.time()
print("The time of execution of above program is :", end-start)