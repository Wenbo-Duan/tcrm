# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 15:26:30 2022

@author: Wenbo.Duan
"""
from os.path import join as pjoin
import logging
import math
import numpy as np
from Utilities.files import flLoadFile
import pandas as pd
import folium
import random

IBTRAC_path = 'C:\\tcrm\\input\\ibtracs.since1980.list.v04r00.csv'

col_list = ["SID", "SEASON",'BASIN','LAT','LON', 'WMO_WIND','USA_WIND','WMO_PRES','USA_PRES','USA_ROCI','USA_RMW']
IBTRAC = pd.read_csv(IBTRAC_path,header=0, usecols=col_list,keep_default_na=False)
IBTRAC = IBTRAC.iloc[1: , :]
IBTRAC = IBTRAC.reset_index(drop=True)
# outputPath = 'C:\\tcrm\\output\\Boston_10by10_5000yr_1980'
outputPath = 'G:\\Advanced Technology & Research\\Climate Analysis\\079007-70 IiA Extreme Winds\\tcrm\\output\\Boston_10by10_5000yr_Current_run1'

cyclone_tracks = flLoadFile(pjoin(outputPath, 'process','cyclone_tracks'),'%', ',')

index = cyclone_tracks[:, 0]
lons = cyclone_tracks[:, 1]
lats = cyclone_tracks[:, 2]
# #boston
# wf_domain = {'xMin': 284, 'xMax': 294, 'yMin': 36.5, 'yMax': 46.5}
# tg_domain = {'xMin': 284, 'xMax': 294, 'yMin': 36.5, 'yMax': 46.5}
# #Wilmingtong
# wf_domain = {'xMin': 277, 'xMax': 287, 'yMin': 29, 'yMax': 39}
# tg_domain = {'xMin': 277, 'xMax': 287, 'yMin': 29, 'yMax': 39}
# #Miami
wf_domain = {'xMin': 275, 'xMax': 285, 'yMin': 20, 'yMax': 30}
tg_domain = {'xMin': 275, 'xMax': 285, 'yMin': 20, 'yMax': 30}
track_limits = {'xMin':9999, 'xMax':-9999,
                'yMin':9999, 'yMax':-9999}
saved_track_idx = []
count = 0

for i, [idx, lon, lat] in enumerate(zip(index, lons, lats)):
   
    if idx == 1:
        # Reset cyclone lon/lon limits
        track_limits = {'xMin':9999, 'xMax':-9999,
                        'yMin':9999, 'yMax':-9999}

    track_limits['xMin'] = min(track_limits['xMin'], lon)
    track_limits['xMax'] = max(track_limits['xMax'], lon)
    track_limits['yMin'] = min(track_limits['yMin'], lat)
    track_limits['yMax'] = max(track_limits['yMax'], lat)
        # if TC pass through the wind field domain - WD Comment
    if (wf_domain['xMin'] <= lon <= wf_domain['xMax']) & \
       (wf_domain['yMin'] <= lat <= wf_domain['yMax']):
           
           saved_track_idx.append(i)
           tg_domain['xMin'] = min(tg_domain['xMin'],
                                        track_limits['xMin'])
           tg_domain['xMax'] = max(tg_domain['xMax'],
                                        track_limits['xMax'])
           tg_domain['yMin'] = min(tg_domain['yMin'],
                                        track_limits['yMin'])
           tg_domain['yMax'] = max(tg_domain['yMax'],
                                        track_limits['yMax'])
    # Extend domain to closest integer lat/lon value
    tg_domain['xMin'] = math.floor(tg_domain['xMin'])
    tg_domain['xMax'] = math.ceil(tg_domain['xMax'])
    tg_domain['yMin'] = math.floor(tg_domain['yMin'])
    tg_domain['yMax'] = math.ceil(tg_domain['yMax'])


basin_check = IBTRAC['BASIN'][np.array(saved_track_idx)].unique()
SID_list = IBTRAC['SID'][np.array(saved_track_idx)].unique() # Unique ID of filtered tracks

my_map = folium.Map(location=[25, -80],zoom_start = 5) #Boston [40, -71], Miami [25, -80]

my_map_marker = folium.Map(location=[25, -80],zoom_start = 5)

for i_track, SID_track in enumerate(SID_list):
    
    Track_data = IBTRAC[IBTRAC['SID'].str.contains(SID_track)]
    
    # First_track = IBTRAC[IBTRAC['SID'].str.contains(SID_list[0])]

    lon_track = Track_data['LON']
    lat_track = Track_data['LAT']
    lon_track = lon_track.to_numpy().astype(float)
    lat_track = lat_track.to_numpy().astype(float)
    coordinate = np.vstack((lat_track,lon_track)).T


    hexadecimal = ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
    folium.PolyLine(coordinate,  popup=SID_track, color=hexadecimal, weight=2.5, opacity=1).add_to(my_map)
    
    folium.Marker(location=[lat_track[0],lon_track[0]], popup=SID_track).add_to(my_map_marker)
    
# {'xMin': 263, 'xMax': 344, 'yMin': 8, 'yMax': 47}

# domain Boston
my_box=folium.PolyLine(locations=[[8, 263-360],[8, 344-360],[47, 344-360],[47, 263-360],[8, 263-360]],weight=2)
my_map.add_child(my_box)

my_box=folium.PolyLine(locations=[[31.5, 279-360],[31.5, 299-360],[51.5, 299-360],[51.5, 279-360],[31.5, 279-360]],weight=2)
my_map.add_child(my_box)

# my_box=folium.PolyLine(locations=[[36.5, -76],[36.5, -66],[46.5, -66],[46.5,-76],[36.5, -76]],weight=2, color = 'red')
# my_map.add_child(my_box)

# # domain Wilmington
# Track_generation_domain=folium.PolyLine(locations=[[8, 262-360],[8, 344-360],[39, 344-360],[39, 262-360],[8, 262-360]],weight=2)
# my_map.add_child(Track_generation_domain)


# Wind_domain=folium.PolyLine(locations=[[29, -83],[29, -73],[39, -73],[39,-83],[29, -83]],weight=2, color = 'red')
# my_map.add_child(Wind_domain)

# # domain Miami
# Track_generation_domain=folium.PolyLine(locations=[[8, 263-360],[8, 340-360],[39, 340-360],[39, 263-360],[8, 263-360]],weight=2)
# my_map.add_child(Track_generation_domain)


# Wind_domain=folium.PolyLine(locations=[[20, -85],[20, -75],[30, -75],[30,-85],[20, -85]],weight=2, color = 'red')
# my_map.add_child(Wind_domain)

my_map.save(outputPath + "\\" + "Historical Tracks.html")
my_map_marker.save(outputPath + "\\" + "Historical hurricanes initial positions.html")


# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
# # Lambert Conformal Conic map.
# m = Basemap(llcrnrlon=-100.,llcrnrlat=0.,urcrnrlon=-20.,urcrnrlat=57.,
#             projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60.,
#             resolution ='l',area_thresh=1000.)

# m.plot(xx,yy,color='k')
# # draw coastlines, meridians and parallels.
# m.drawcoastlines()
# m.drawcountries()
# m.drawmapboundary(fill_color='#99ffff')
# m.fillcontinents(color='#cc9966',lake_color='#99ffff')
# m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
# m.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
# plt.title('Atlantic Hurricane Tracks (Storms Reaching Category 4, 1851-2004)')
# plt.show()