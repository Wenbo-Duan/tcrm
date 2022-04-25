# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 13:17:13 2022

@author: Wenbo.Duan
"""

from Utilities.track import ncReadTrackData, Track
from os.path import join as pjoin
import folium
import numpy as np

import os



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

# read all tracks that passes through region
windfieldpath= 'C:/tcrm/output/Boston_10by10_5000yr_1980/windfield/'
arr = os.listdir(windfieldpath)
index_track = np.zeros(len(arr))
index_year = []
for i, file in enumerate(arr):
    gust,index,nc = file.split('.')
    i_track,i_year = index.split('-')
    index_track[i] = i_track
    index_year.append(i_year)
    
    

trackpath= 'C:/tcrm/output/Boston_10by10_5000yr_1980/tracks/'


# trackname = 'tracks.02246.nc' # Boston largest wind speed
trackname = 'tracks.03277.nc' # Boston 2nd largest wind speed
for i, year in enumerate(index_year):
    
    trackname = 'tracks.' + year + '.nc' # Boston 2nd largest wind speed
    
    trackfile = pjoin(trackpath,trackname)

    trackiter = loadTracks(trackfile)

    SubTrack = index_track[i]  # Boston largest wind speed

    # SubTrack = 6  # Boston 2nd largest wind speed

    Track_Data = trackiter[int(SubTrack)]

    # Track.Longitude
    
    # Coordinate= [Track.Latitude  Track.Longitude -360]
    Coordinate = np.vstack((Track_Data.Latitude , Track_Data.Longitude -360)).T
    
    m = folium.Map(location=[42.4, -71])
    
    # n_step = len(Track.Longitude)
    my_PolyLine=folium.PolyLine(locations=Coordinate,weight=5)
    m.add_child(my_PolyLine)
    [279, 299, 31.5, 51.5]
    
    my_box=folium.PolyLine(locations=[[31.5, 279-360],[31.5, 299-360],[51.5, 299-360],[51.5, 279-360],[31.5, 279-360]],weight=2)
    m.add_child(my_box)
    
    my_box=folium.PolyLine(locations=[[36.5, -76],[36.5, -66],[46.5, -66],[46.5,-76],[36.5, -76]],weight=2, color = '#ff4433')
    m.add_child(my_box)
    
    for n_step in range(len(Track_Data.Longitude)):
        folium.Circle(
            radius=Track_Data.rMax[n_step]*1000,
            location=[Track_Data.Latitude[n_step], Track_Data.Longitude[n_step] - 360],
            popup=Track_Data.TimeElapsed[n_step],
            color="crimson",
            fill=False,
            ).add_to(m)
        
        
    
    m.save("Boston_all_track_map.html")