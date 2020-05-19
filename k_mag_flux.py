#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:40:22 2019

Module to contain the functions to create flux or magnitude arrays in K
Also contains functions to create error arrays using both SExtractor and our
self-calibrated method.

Created as part of a restucture of my code to make functions easier to find

@author: ppxee
"""
import numpy as np #for handling arrays
import field_funcs # to allow error creation by quadrant

def flux_stacks(tbdata, aper=5):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
    aper -= 1 # to make it zero indexed
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,aper], #tbdata['FLUX_APER_06B'][:,aper],
                tbdata['FLUX_APER_07B'][:,aper], tbdata['FLUX_APER_08B'][:,aper],
                tbdata['FLUX_APER_09B'][:,aper], tbdata['FLUX_APER_10B'][:,aper], 
                tbdata['FLUX_APER_11B'][:,aper], tbdata['FLUX_APER_12B'][:,aper]]), axis=1)
    return flux

def fluxerr_stacks(tbdata, aper=5):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
    aper -= 1 # to make it zero indexed
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,aper],#tbdata['FLUXERR_APER_06B'][:,aper],
                tbdata['FLUXERR_APER_07B'][:,aper],tbdata['FLUXERR_APER_08B'][:,aper],
                tbdata['FLUXERR_APER_09B'][:,aper],tbdata['FLUXERR_APER_10B'][:,aper], 
                tbdata['FLUXERR_APER_11B'][:,aper],tbdata['FLUXERR_APER_12B'][:,aper]], axis=1)
    return fluxerr

def month_flux_stacks(tbdata, aper=4):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
    aper -= 1 # to make it zero indexed
    months = ['sep05','oct05','nov05','dec05', 'jan06', #'dec06', 
              'jan07', 'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 
              'jul09', 'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 
              'feb10', 'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', #'feb11', 
              'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
              'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
#    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
#          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
#          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
#          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
#          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
#          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    for month in months:
        if month == 'sep05':
            flux = tbdata['FLUX_APER_'+month][:,aper]
        else:
            flux = np.vstack([flux, tbdata['FLUX_APER_'+month][:,aper]])
    return np.transpose(flux)

def month_fluxerr_stacks(tbdata, aper=4):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 40 columns containing flux values for each year '''
    aper -= 1 # to make it zero indexed
    months = ['sep05','oct05','nov05','dec05', 'jan06', #'dec06', 
              'jan07', 'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 
              'jul09', 'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 
              'feb10', 'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', #'feb11', 
              'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
              'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
#    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
#          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
#          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
#          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
#          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
#          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    for month in months:
        if month == 'sep05':
            fluxerr = tbdata['FLUXERR_APER_'+month][:,aper]
        else:
            fluxerr = np.vstack([fluxerr, tbdata['FLUXERR_APER_'+month][:,aper]])
    return np.transpose(fluxerr)

def mag_stacks(tbdata, aper=5):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        mag = an array with 8 columns containing mag values for each year '''
    aper -= 1 # to make it zero indexed
    mag = np.stack(([tbdata['MAG_APER_05B'][:,aper], #tbdata['MAG_APER_06B'][:,aper],
                tbdata['MAG_APER_07B'][:,aper], tbdata['MAG_APER_08B'][:,aper],
                tbdata['MAG_APER_09B'][:,aper], tbdata['MAG_APER_10B'][:,aper], 
                tbdata['MAG_APER_11B'][:,aper], tbdata['MAG_APER_12B'][:,aper]]), axis=1)
    return mag

def magerr_stacks(tbdata, aper=5):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing mag error values for each
        year '''
    aper -= 1 # to make it zero indexed
    magerr = np.stack([tbdata['MAGERR_APER_05B'][:,aper],#tbdata['MAGERR_APER_06B'][:,aper],
                tbdata['MAGERR_APER_07B'][:,aper],tbdata['MAGERR_APER_08B'][:,aper],
                tbdata['MAGERR_APER_09B'][:,aper],tbdata['MAGERR_APER_10B'][:,aper], 
                tbdata['MAGERR_APER_11B'][:,aper],tbdata['MAGERR_APER_12B'][:,aper]], axis=1)
    return magerr

def month_mag_stacks(tbdata, aper=4):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
    aper -= 1 # to make it zero indexed
    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    for month in months:
        if month == 'sep05':
            mag = tbdata['MAG_APER_'+month][:,aper]
        else:
            mag = np.vstack([mag, tbdata['MAG_APER_'+month][:,aper]])
    return mag

def month_magerr_stacks(tbdata, aper=4):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the specified aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
    aper -= 1 # to make it zero indexed
    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    for month in months:
        if month == 'sep05':
            magerr = tbdata['MAGERR_APER_'+month][:,aper]
        else:
            magerr = np.vstack([magerr, tbdata['MAGERR_APER_'+month][:,aper]])
    return magerr

def flux1_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 0.67 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,0], #tbdata['FLUX_APER_06B'][:,0],
                tbdata['FLUX_APER_07B'][:,0], tbdata['FLUX_APER_08B'][:,0],
                tbdata['FLUX_APER_09B'][:,0], tbdata['FLUX_APER_10B'][:,0], 
                tbdata['FLUX_APER_11B'][:,0], tbdata['FLUX_APER_12B'][:,0]]), axis=1)
    return flux

def fluxerr1_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 0.67 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,0],#tbdata['FLUXERR_APER_06B'][:,0],
                tbdata['FLUXERR_APER_07B'][:,0],tbdata['FLUXERR_APER_08B'][:,0],
                tbdata['FLUXERR_APER_09B'][:,0],tbdata['FLUXERR_APER_10B'][:,0], 
                tbdata['FLUXERR_APER_11B'][:,0],tbdata['FLUXERR_APER_12B'][:,0]], axis=1)
    return fluxerr

def flux2_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 1 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,1], #tbdata['FLUX_APER_06B'][:,1],
                tbdata['FLUX_APER_07B'][:,1], tbdata['FLUX_APER_08B'][:,1],
                tbdata['FLUX_APER_09B'][:,1], tbdata['FLUX_APER_10B'][:,1], 
                tbdata['FLUX_APER_11B'][:,1], tbdata['FLUX_APER_12B'][:,1]]), axis=1)
    return flux

def fluxerr2_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 1 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,1],#tbdata['FLUXERR_APER_06B'][:,1],
                tbdata['FLUXERR_APER_07B'][:,1],tbdata['FLUXERR_APER_08B'][:,1],
                tbdata['FLUXERR_APER_09B'][:,1],tbdata['FLUXERR_APER_10B'][:,1], 
                tbdata['FLUXERR_APER_11B'][:,1],tbdata['FLUXERR_APER_12B'][:,1]], axis=1)
    return fluxerr

def flux3_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 1.5 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,2],# tbdata['FLUX_APER_06B'][:,2],
                tbdata['FLUX_APER_07B'][:,2], tbdata['FLUX_APER_08B'][:,2],
                tbdata['FLUX_APER_09B'][:,2], tbdata['FLUX_APER_10B'][:,2], 
                tbdata['FLUX_APER_11B'][:,2], tbdata['FLUX_APER_12B'][:,2]]), axis=1)
    return flux

def fluxerr3_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 1.5 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,2],#tbdata['FLUXERR_APER_06B'][:,2],
                tbdata['FLUXERR_APER_07B'][:,2],tbdata['FLUXERR_APER_08B'][:,2],
                tbdata['FLUXERR_APER_09B'][:,2],tbdata['FLUXERR_APER_10B'][:,2], 
                tbdata['FLUXERR_APER_11B'][:,2],tbdata['FLUXERR_APER_12B'][:,2]], axis=1)
    return fluxerr

def flux4_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 2 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,3], #tbdata['FLUX_APER_06B'][:,3],
                tbdata['FLUX_APER_07B'][:,3], tbdata['FLUX_APER_08B'][:,3],
                tbdata['FLUX_APER_09B'][:,3], tbdata['FLUX_APER_10B'][:,3], 
                tbdata['FLUX_APER_11B'][:,3], tbdata['FLUX_APER_12B'][:,3]]), axis=1)
    return flux

def fluxerr4_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 2 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,3],#tbdata['FLUXERR_APER_06B'][:,3],
                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
                tbdata['FLUXERR_APER_09B'][:,3],tbdata['FLUXERR_APER_10B'][:,3], 
                tbdata['FLUXERR_APER_11B'][:,3],tbdata['FLUXERR_APER_12B'][:,3]], axis=1)
    return fluxerr

def flux5_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 3 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,4], #tbdata['FLUX_APER_06B'][:,4],
                tbdata['FLUX_APER_07B'][:,4], tbdata['FLUX_APER_08B'][:,4],
                tbdata['FLUX_APER_09B'][:,4], tbdata['FLUX_APER_10B'][:,4], 
                tbdata['FLUX_APER_11B'][:,4], tbdata['FLUX_APER_12B'][:,4]]), axis=1)
    return flux

def fluxerr5_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 3 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,4],#tbdata['FLUXERR_APER_06B'][:,4],
                tbdata['FLUXERR_APER_07B'][:,4],tbdata['FLUXERR_APER_08B'][:,4],
                tbdata['FLUXERR_APER_09B'][:,4],tbdata['FLUXERR_APER_10B'][:,4], 
                tbdata['FLUXERR_APER_11B'][:,4],tbdata['FLUXERR_APER_12B'][:,4]], axis=1)
    return fluxerr


def create_error_array(flux, sigtb, tbdata):
    ''' Function that creates an error array from sigma values calculated from
    data variations within flux bins 
    *** This is for an out of date error determination ***
    Inputs:
        flux = array of flux values to find errors for
        sigtb = table that defines the error values within named bins
        tbdata = overall mag-flux table of the data
    Outputs:
        flux = array of flux that matches errarr
        errarr = array of corrected flux errors
        tbdata = overall table with any null values removed
        '''
    bins = np.array(sigtb.colnames).astype(int)
    avgflux = np.nanmean(flux, axis=1)
    errarr = np.empty(np.shape(flux))
    for n, lower in enumerate(bins):
        if n == 0:
            #remove lower values
            mask = avgflux< lower
            nanarr = np.full_like(sigtb[str(lower)], np.nan)
            flux[mask,:] = nanarr
            errarr[mask,:] = nanarr
        mask1 = avgflux>lower 
        if n != len(bins)-1:
            mask2 = avgflux<bins[n+1]
        else:
            #remove lower values
            mask = avgflux > 6309573 #max bin
            nanarr = np.full_like(sigtb[str(lower)], np.nan)
            flux[mask,:] = nanarr
            errarr[mask,:] = nanarr
            mask2 = avgflux < 6309573 #max bin
        
        mask = mask1*mask2.astype(bool)
        errarr[mask,:] = sigtb[str(lower)]
        
    # remove nans
    mask = ~np.isnan(flux).any(axis=1)
    flux = flux[mask]
    errarr = errarr[mask]
    tbdata = tbdata[mask]
    return flux, errarr, tbdata

def create_quad_error_array(sigtb, tbdata, aper=5, quadoutput=False):
    ''' Function that creates an error array from sigma values calculated from
    data variations within flux, epoch and quadrant bins - i.e. the error 
    analysis used in my first paper 
    
    Inputs:
        sigtb = table of error values, column names indicate which bin the  
                value applies to.
        tbdata = table of data the array is being created for
        aper = aperture to create flux arrays for, default is 5 which is the 
               3 arcsec aperture
        quadoutput = bool, indicate whether to output the data in 4 separate
                     arrays corresponding to the four quadrants. Default is 
                     false so output is as single arrays and tables.
    Outputs:
        If quadoutput = False (default):
            flux = array of flux values that matches errarr
            erarr = array of calibrated error values
            newtbdata = overall data table with any null values removed so it 
                        matches flux and errarr in length
        If quadoutput = True:
            quadflux = dictionary with 4 arrays inside that are for flux arrays
                       for each quadrant of the image. Dict key is the quadrant
                       number.
            quaderr = dictionary with 4 arrays inside that are for error arrays
                      for each quadrant of the image. Dict key is the quadrant
                      number.
            newquaddata = dictionary with 4 tables inside that are for the 
                          overall data tables for each quadrant of the image. 
                          Dict key is the quadrant number.
    '''
    bins = np.array(sigtb.colnames)
    binarr = np.empty([4,int(len(bins)/4)])
    k = 0
    l = 0
    m = 0
    n = 0
    ### Create binedge array for each quadrant ###
    for bin in bins:
        if bin[0] == '1':
            binarr[0,k] = int(bin[2:])
            k+=1
        elif bin[0] == '2':
            binarr[1,l] = int(bin[2:])
            l+=1
        elif bin[0] == '3':
            binarr[2,m] = int(bin[2:])
            m+=1
        elif bin[0] == '4':
            binarr[3,n] = int(bin[2:])
            n+=1
    
    ### Set up empty arrays for data ###
    flux = np.array([])
    errarr = np.array([])
    newtbdata = []
    quadflux = {}
    quaderr = {}
    newquaddata = {}
    
    ### Get quadrant data ###
    quaddata = field_funcs.quadrants(tbdata, '05B')
    for n, qdata in enumerate(quaddata):
        ### create flux stacks and find average
        qflux = flux_stacks(qdata, aper)
        avgflux = np.nanmean(qflux, axis=1)
        qerrarr = np.empty(np.shape(qflux))
        
        ### Find values within bins and assign correct sigma ###
        for m, lower in enumerate(binarr[n,0:-1]):
            if m == 0:
                #remove lower values
                mask = avgflux< lower
                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
                qflux[mask,:] = nanarr
                qerrarr[mask,:] = nanarr
            mask1 = avgflux>lower 
            if m != len(binarr)-1:
                mask2 = avgflux<binarr[n,m+1]
            else:
                #remove lower values
                mask = avgflux > 6309573 #max bin
                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
                qflux[mask,:] = nanarr
                qerrarr[mask,:] = nanarr
                mask2 = avgflux < 6309573 #max bin
            
            mask = mask1*mask2.astype(bool)
            qerrarr[mask,:] = sigtb[str(n+1)+' '+str(int(lower))]
            
        # remove nans
        mask = ~np.isnan(qflux).any(axis=1)
        qflux = qflux[mask]
        qerrarr = qerrarr[mask]
        qdata = qdata[mask]
        
        if quadoutput == False:
            ### Define full arrays to use in the rest of the analysis ###
            if n == 0:
                flux = np.copy(qflux)
                errarr = np.copy(qerrarr)
                newtbdata = np.copy(qdata)
            else:
                flux = np.vstack((flux, np.copy(qflux)))
                errarr = np.vstack((errarr, np.copy(qerrarr)))
                newtbdata = np.hstack((newtbdata, np.copy(qdata)))
        
        else: #this means wants quadrant output
            quadflux[n] = qflux
            quaderr[n] = qerrarr
            newquaddata[n] = qdata
    
    if quadoutput == False:
        if np.isin('X-ray', newtbdata.dtype.names):
            newtbdata['X-ray'][newtbdata['X-ray']==70] = False 
            newtbdata['X-ray'][newtbdata['X-ray']==84] = True
            newtbdata['X-ray'] = newtbdata['X-ray'].astype(bool)
        return flux, errarr, newtbdata
    else:
        return quadflux, quaderr, newquaddata

def create_quad_error_array_month(sigtb, tbdata, aper=4, quadoutput=False):
    ''' Function that creates an error array from sigma values calculated from
    data variations within flux, epoch and quadrant bins - i.e. the error 
    analysis used in my first paper - this version is for month stacks
    
    Inputs:
        sigtb = table of error values, column names indicate which bin the  
                value applies to.
        tbdata = table of data the array is being created for
        aper = aperture to create flux arrays for, default is 4 which is the 
               2 arcsec aperture
        quadoutput = bool, indicate whether to output the data in 4 separate
                     arrays corresponding to the four quadrants. Default is 
                     false so output is as single arrays and tables.
    Outputs:
        If quadoutput = False (default):
            flux = array of flux values that matches errarr
            erarr = array of calibrated error values
            newtbdata = overall data table with any null values removed so it 
                        matches flux and errarr in length
        If quadoutput = True:
            quadflux = dictionary with 4 arrays inside that are for flux arrays
                       for each quadrant of the image. Dict key is the quadrant
                       number.
            quaderr = dictionary with 4 arrays inside that are for error arrays
                      for each quadrant of the image. Dict key is the quadrant
                      number.
            newquaddata = dictionary with 4 tables inside that are for the 
                          overall data tables for each quadrant of the image. 
                          Dict key is the quadrant number.
    '''
    bins = np.array(sigtb.colnames)
    binarr = np.empty([4,int(len(bins)/4)])
    k = 0
    l = 0
    m = 0
    n = 0
    ### Create binedge array for each quadrant ###
    for bin in bins:
        if bin[0] == '1':
            binarr[0,k] = int(bin[2:])
            k+=1
        elif bin[0] == '2':
            binarr[1,l] = int(bin[2:])
            l+=1
        elif bin[0] == '3':
            binarr[2,m] = int(bin[2:])
            m+=1
        elif bin[0] == '4':
            binarr[3,n] = int(bin[2:])
            n+=1
    
    ### Set up empty arrays for data ###
    flux = np.array([])
    errarr = np.array([])
    newtbdata = []
    quadflux = {}
    quaderr = {}
    newquaddata = {}
    
    ### Get quadrant data ###
#    print('Splitting quadrants')
    quaddata = field_funcs.quadrants(tbdata, 'sep05')
#    print('Split quadrants')
    for n, qdata in enumerate(quaddata):
        ### create flux stacks and find average
        qflux = month_flux_stacks(qdata, aper)
        avgflux = np.nanmean(qflux, axis=1)
        qerrarr = np.empty(np.shape(qflux))
        
        ### Find values within bins and assign correct sigma ###
        for m, lower in enumerate(binarr[n,0:-1]):
            if m == 0:
                #remove lower values
                mask = avgflux< lower
                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
                qflux[mask,:] = nanarr
                qerrarr[mask,:] = nanarr
            mask1 = avgflux>lower 
            if m != len(binarr)-1:
                mask2 = avgflux<binarr[n,m+1]
            else:
                #remove lower values
                mask = avgflux > 6309573 #max bin
                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
                qflux[mask,:] = nanarr
                qerrarr[mask,:] = nanarr
                mask2 = avgflux < 6309573 #max bin
            
            mask = mask1*mask2.astype(bool)
            qerrarr[mask,:] = sigtb[str(n+1)+' '+str(int(lower))]
            
        # remove nans
        mask = ~np.isnan(qflux).any(axis=1)
        qflux = qflux[mask]
        qerrarr = qerrarr[mask]
        qdata = qdata[mask]
        
        if quadoutput == False:
            ### Define full arrays to use in the rest of the analysis ###
#            print('stacking tables')
            if n == 0:
                flux = np.copy(qflux)
                errarr = np.copy(qerrarr)
                newtbdata = np.copy(qdata)
            else:
                flux = np.vstack((flux, np.copy(qflux)))
                errarr = np.vstack((errarr, np.copy(qerrarr)))
                newtbdata = np.hstack((newtbdata, np.copy(qdata)))
#            print('tables stacked')
        
        else: #this means wants quadrant output
            quadflux[n] = qflux
            quaderr[n] = qerrarr
            newquaddata[n] = qdata
    
    if quadoutput == False:
        if np.isin('X-ray_1', newtbdata.dtype.names):
            newtbdata['X-ray_1'][newtbdata['X-ray_1']==70] = False 
            newtbdata['X-ray_1'][newtbdata['X-ray_1']==84] = True
            newtbdata['X-ray_1'] = newtbdata['X-ray_1'].astype(bool)
            newtbdata['X-ray_2'][newtbdata['X-ray_2']==70] = False 
            newtbdata['X-ray_2'][newtbdata['X-ray_2']==84] = True
            newtbdata['X-ray_2'] = newtbdata['X-ray_2'].astype(bool)
        return flux, errarr, newtbdata
    else:
        return quadflux, quaderr, newquaddata
