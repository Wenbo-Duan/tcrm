# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 21:56:01 2022

@author: Wenbo.Duan
"""

import numpy as np
Valid_Speed = np.zeros((100,100))
for i in range(100):
    for j in range(100):
        Valid_Speed[i,j] = np.count_nonzero(Vr[:,i,j]>0)
        
np.min(Valid_Speed)


count = 0
# 
for i in range(0,len(Vr)):
    if np.count_nonzero(Vr[i,:,:]) == 0:
        count += 1
        
# x,y,z=np.shape(Vr)
# Vr_2D = Vr.reshape(x, y*z)

np.save("Miami_Vr", Vr)

# a = np.load('Wilmington_Vr.npy')