#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:06:44 2017

First attempt to plot a light curve from yearly stack data.

@author: ppxee
"""

#Import required libraries
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
plt.close('all') #close any open plots


#ob = np.array([93634, 123294, 135921, 201338, 236901]) #array of clear anomalies from 3 arcsec
#ob = np.array([293821, 135921, 48769, 44190, 236901]) #array of clear(ish) anomalies from 2 arcsec
#ob = np.array([164877, 118496, 81012]) #array of clear anomalies in the chandra agn 2 arcsec

# Open the fits file and get data
combined = fits.open('combined_mag_flux_table.fits')
tbdata = combined[1].data
def lightcurve(ob)  :
    for n in range(0, np.size(ob)):
        
        #Get data for a single object in the field
        mask = tbdata['NUMBER_05B'] == ob[n]
        obdata = tbdata[mask]
        
        #Create array of mag values
#        mag = np.array([obdata['MAG_APER_4_05B'],obdata['MAG_APER_4_06B'],
#                        obdata['MAG_APER_4_07B'],obdata['MAG_APER_4_08B'],
#                        obdata['MAG_APER_4_09B'],obdata['MAG_APER_4_10B'], 
#                        obdata['MAG_APER_4_11B'],obdata['MAG_APER_4_12B']])
#        
#        magerr = np.array([obdata['MAGERR_APER_4_05B'],obdata['MAGERR_APER_4_06B'],
#                        obdata['MAGERR_APER_4_07B'],obdata['MAGERR_APER_4_08B'],
#                        obdata['MAGERR_APER_4_09B'],obdata['MAGERR_APER_4_10B'], 
#                        obdata['MAGERR_APER_4_11B'],obdata['MAGERR_APER_4_12B']])
        
        #Create array of flux values
        flux = np.array([obdata['FLUX_APER_4_05B'],obdata['FLUX_APER_4_06B'],
                        obdata['FLUX_APER_4_07B'],obdata['FLUX_APER_4_08B'],
                        obdata['FLUX_APER_4_09B'],obdata['FLUX_APER_4_10B'], 
                        obdata['FLUX_APER_4_11B'],obdata['FLUX_APER_4_12B']])
        
        fluxerr = np.array([obdata['FLUXERR_APER_4_05B'],obdata['FLUXERR_APER_4_06B'],
                        obdata['FLUXERR_APER_4_07B'],obdata['FLUXERR_APER_4_08B'],
                        obdata['FLUXERR_APER_4_09B'],obdata['FLUXERR_APER_4_10B'], 
                        obdata['FLUXERR_APER_4_11B'],obdata['FLUXERR_APER_4_12B']])
        
        #set up time variable for plot
        t = np.linspace(1, 8, num=8)
        years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
        
        #plot points with errorbars for both mag and flux
        #plt.figure(1)
        #plt.xticks(t, years)
        #plt.errorbar(t, mag, yerr=magerr, fmt = 'bo')
        #plt.xlabel('Semester')
        #plt.ylabel('Magnitude')
        #plt.title('Light curve using magnitude')
        #
        plt.figure(n+1)
        plt.xticks(t, years)
        plt.errorbar(t, flux, yerr=fluxerr, fmt = 'ro')
        plt.xlabel('Semester')
        plt.ylabel('Flux of object')
        plt.title('Light curve using flux')

lightcurve(203034)