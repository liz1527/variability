#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:50:15 2019

Code containing functions that involve dividing up or restricting the area that
is being studied (e.g. slice into quads or restrict to chandra)

@author: ppxee
"""

def quadrants(initdata,sem='05B'):
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
    mask2 = idec >= middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec < middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad4data = initdata[mask1*mask2]    
    
    return quad1data, quad2data, quad3data, quad4data


def chandra_only(tbdata, sem='05B'):
    ''' Function that restricts the objects included in analysis to only 
    those within the Chandra footprint 
    Input:
        tbdata = original catalogue of data 
        sem = semester or month on which the split should be done
    Output:
        newtbdata = new catalogue of data which only includes objects within 
                    the chandra footprint 
    '''
    if sem == '': # allows variable tables to be limited when this is entered
        ### Restrict objects to those in the Chandra field ###
        mask1 = tbdata['Dec'] < -4.93 #max Dec
        mask2 = tbdata['Dec'] > -5.403 #min Dec
        mask3 = tbdata['RA'] < 34.72 #max RA
        mask4 = tbdata['RA'] > 34.07 #min RA
    else:
        ### Restrict objects to those in the Chandra field ###
        mask1 = tbdata['DELTA_J2000_'+sem] < -4.93 #max Dec
        mask2 = tbdata['DELTA_J2000_'+sem] > -5.403 #min Dec
        mask3 = tbdata['ALPHA_J2000_'+sem] < 34.72 #max RA
        mask4 = tbdata['ALPHA_J2000_'+sem] > 34.07 #min RA
    
    mask = mask1 * mask2 * mask3 * mask4
    
    newtbdata = tbdata[mask]
    
    return newtbdata


def remove_edges(tbdata, sem='05B'):
    ''' Function that masks the edges of the image as there were a large number
    of variables in this area originally so must be too much noise here.
    Inputs:
        tbdata = orginal table of data covering the whole field
        sem = semester or month on which the split should be done
    Ouputs:
        tbdata = table of date with objects need the edges removed
    '''
    
    ### Set X limits ###
    x = tbdata['X_IMAGE_'+sem]
    xmask1 = x > 1000
    xmask2 = x < 24000
    xmask = xmask1 * xmask2.astype(bool)
    tbdata = tbdata[xmask]
    
    ### Set Y limits ###
    y = tbdata['Y_IMAGE_'+sem]
    ymask1 = y > 2000
    ymask2 = y < 24500
    ymask = ymask1 * ymask2.astype(bool)
    tbdata = tbdata[ymask]
    
    return tbdata