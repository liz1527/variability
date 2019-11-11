#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:50:15 2019

Code containing functions that involve dividing up or restricting the area that
is being studied (e.g. slice into quads or restrict to chandra)

@author: ppxee
"""

def quadrants(initdata,sem):
    ''' Function that splits the field up into the 4 quadrants covered by the 
    4 camera on WFCAM
    
    Inputs:
        initdata = intial table of data for the whole field
        sem = the semester that you want to divide along. This shouldn't make 
              much of a difference so usually 05B is used but is here just in 
              case I need to be more specific at some point.
    
    Outputs:
        quad1data = data for the top left quad (according to XY coords)
        quad2data = data for the top right quad (according to XY coords)
        quad3data = data for the bottom left quad (according to XY coords)
        quad4data = data for the bottom right quad (according to XY coords)
    '''
    
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


def chandra_only(tbdata):
    ''' Function that restricts the objects included in analysis to only 
    those within the Chandra footprint 
    Input:
        tbdata = original catalogue of data 
    Output:
        newtbdata = new catalogue of data which only includes objects within 
                    the chandra footprint 
    '''
    ### Restrict objects to those in the Chandra field ###
    mask1 = tbdata['DELTA_J2000_05B'] < -4.93 #max Dec
    mask2 = tbdata['DELTA_J2000_05B'] > -5.403 #min Dec
    mask3 = tbdata['ALPHA_J2000_05B'] < 34.72 #max RA
    mask4 = tbdata['ALPHA_J2000_05B'] > 34.07 #min RA
    
    mask = mask1 * mask2 * mask3 * mask4
    
    newtbdata = tbdata[mask]
    
    return newtbdata