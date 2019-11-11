#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:50:15 2019

Code containing functions that involve working on quadrants of the image, 
including the quadrants function itself

@author: ppxee
"""

def quadrants(initdata,sem):
    
    ira = initdata['X_IMAGE_'+sem]
    idec = initdata['Y_IMAGE_'+sem]

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec < middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec >= middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad4data = initdata[mask1*mask2]
    
    return quad1data, quad2data, quad3data, quad4data