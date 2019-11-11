#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:10:31 2019

Module containing functions that plot lightcurves

@author: ppxee
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table #for handling tables

def lightcurve4(ob, fitsdata)  :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 4 pixels across (i.e. 2 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
        
    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- invalid object number')
        return
    
    #Create arrays of flux values and error values
    flux = np.array([obdata['FLUX_APER_05B'][:,4],obdata['FLUX_APER_06B'][:,4],
                    obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
                    obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
                    obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
    
    fluxerr = np.array([obdata['FLUXERR_APER_05B'][:,4],obdata['FLUXERR_APER_06B'][:,4],
                    obdata['FLUXERR_APER_07B'][:,4],obdata['FLUXERR_APER_08B'][:,4],
                    obdata['FLUXERR_APER_09B'][:,4],obdata['FLUXERR_APER_10B'][:,4], 
                    obdata['FLUXERR_APER_11B'][:,4],obdata['FLUXERR_APER_12B'][:,4]])
    
    # normalise and correct for seeing
    flux = psf_correct(flux, flux, 'mean')
    fluxerr = psf_correct(flux, fluxerr, 'mean')
    avgflux =np.mean(flux)
    normflux = flux/avgflux
    normerr = fluxerr/avgflux
    
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    #Plot graph in new figure
    plt.figure()
    plt.xticks(t, years)
    plt.errorbar(t, normflux, yerr=normerr, fmt = 'ro')
    plt.xlabel('Semester')
    plt.ylabel('Flux of object')
    plt.title('Light curve for object number %i' % ob)
    return

def lightcurve5(ob, fitsdata)  :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
        
    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- '+str(ob)+' is an invalid object number')
        return
    
    #Create arrays of flux values and error values
    flux = np.array([obdata['MAG_APER_05B'][:,4],obdata['MAG_APER_06B'][:,4],
                    obdata['MAG_APER_07B'][:,4],obdata['MAG_APER_08B'][:,4],
                    obdata['MAG_APER_09B'][:,4],obdata['MAG_APER_10B'][:,4], 
                    obdata['MAG_APER_11B'][:,4],obdata['MAG_APER_12B'][:,4]])
    
    fluxerr = np.array([obdata['MAGERR_APER_05B'][:,4],obdata['MAGERR_APER_06B'][:,4],
                    obdata['MAGERR_APER_07B'][:,4],obdata['MAGERR_APER_08B'][:,4],
                    obdata['MAGERR_APER_09B'][:,4],obdata['MAGERR_APER_10B'][:,4], 
                    obdata['MAGERR_APER_11B'][:,4],obdata['MAGERR_APER_12B'][:,4]])

    
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    #Plot graph in new figure
    plt.figure()
    plt.xticks(t, years)
    plt.errorbar(t, flux, yerr=fluxerr, fmt = 'ro')
    plt.xlabel('Semester')
    plt.ylabel('K-band magnitude of object')
    plt.title('Light curve for object number %i' % ob)
    return

def lightcurveflux5(ob, fitsdata, corrected=False, new_fig=True) :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
    sigtb = Table.read('quad_epoch_sigma_table_extra_clean.fits')

    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- invalid object number')
        return
    
    #Create arrays of flux values and error values
    flux = np.array([obdata['FLUX_APER_05B'][:,4],#obdata['FLUX_APER_06B'][:,4],
                    obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
                    obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
                    obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
    flux, fluxerr, obdata = create_quad_error_array(sigtb, obdata)
#    if corrected == False:
#        fluxerr = np.array([obdata['FLUXERR_APER_05B'][:,4],obdata['FLUXERR_APER_06B'][:,4],
#                            obdata['FLUXERR_APER_07B'][:,4],obdata['FLUXERR_APER_08B'][:,4],
#                            obdata['FLUXERR_APER_09B'][:,4],obdata['FLUXERR_APER_10B'][:,4], 
#                            obdata['FLUXERR_APER_11B'][:,4],obdata['FLUXERR_APER_12B'][:,4]])
#    else:
#        fluxerr = np.array([obdata['FLUXERR_APER_05B_CORR'],obdata['FLUXERR_APER_06B_CORR'],
#                            obdata['FLUXERR_APER_07B_CORR'],obdata['FLUXERR_APER_08B_CORR'],
#                            obdata['FLUXERR_APER_09B_CORR'],obdata['FLUXERR_APER_10B_CORR'], 
#                            obdata['FLUXERR_APER_11B_CORR'],obdata['FLUXERR_APER_12B_CORR']])
    
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    x = [1,3,4,5,6,7,8]
    print('Plotting lightcurve')
    
    #Plot graph in new figure
    if new_fig == True:
        plt.figure()
    plt.xticks(t, years)
    plt.errorbar(x, flux, yerr=fluxerr, fmt = 'ro')
    plt.xlabel('Semester')
    plt.ylabel('K-band flux of object')
    plt.title('Light curve for object number %i' % ob)
    return
