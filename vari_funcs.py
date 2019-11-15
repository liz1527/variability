#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:32:41 2017

Code containing all functions that will be useful in codes looking at 
variability and plotting light curves.

This should be imported at the start of the majority of my codes.

Updated on 27/2/18 to change the way the mag_flux tables are read and to bring 
in line with updates done on vari_funcs_no06 while 06 was not being used.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
from scipy.stats import chisquare

plt.style.use('default')
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Read in sub-modules ###
import k_mag_flux #for flux and magnitude arrays in K
import j_mag_flux #for flux and magnitude arrays in J
import h_mag_flux #for flux and magnitude arrays in H

import field_funcs #for dividing or restricting the area of the field 

import lightcurve_funcs #for plotting lightcurve in various forms

import correction_funcs #for correcting and testing seeing/PSFs

import vary_stats #for variability statistics

import flux_funcs #for changes to flux/mag stacks

import selection_plot_funcs #for plotting main selection plot + onclicks

### Define miscellanseous functions that don't fit in a sub-module ###
def get_z(tbdata):
    ''' Get z with spectroscopic data where possible 
    Inputs:
        tbdata = table of data for all the objects we want z for (should contain
                 DR11 data in order to get z)
    Outputs:
        z = array of redshifts with spectroscopic redshifts where possible and 
            z_p where spec is not available.
    '''
    z = tbdata['z_spec']
    z[z==-1] = tbdata['z_p'][z==-1]
    return z

def get_jansky_flux(tbdata):
    ''' Function to get flux in Janskys from the DR11 AB magnitude of an object 
    Inputs:
        tbdata = table of DR11 data for the objects you want Janksy fluxes for
    Outputs:
        meanflux = flux in Janksky calculated from the AB magnitude in 2 arcsec
                   aperture from teh DR11 catalogue.
    '''
    meanmag = tbdata['KMAG_20']
    meanflux = 10**(23-((meanmag+48.6)/2.5))
    return meanflux

### Define Functions ###
#def flux_stacks(tbdata, aper=5):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#    aper -= 1 # to make it zero indexed
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,aper], #tbdata['FLUX_APER_06B'][:,aper],
#                tbdata['FLUX_APER_07B'][:,aper], tbdata['FLUX_APER_08B'][:,aper],
#                tbdata['FLUX_APER_09B'][:,aper], tbdata['FLUX_APER_10B'][:,aper], 
#                tbdata['FLUX_APER_11B'][:,aper], tbdata['FLUX_APER_12B'][:,aper]]), axis=1)
#    return flux
#
#def fluxerr_stacks(tbdata, aper=5):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#    aper -= 1 # to make it zero indexed
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,aper],#tbdata['FLUXERR_APER_06B'][:,aper],
#                tbdata['FLUXERR_APER_07B'][:,aper],tbdata['FLUXERR_APER_08B'][:,aper],
#                tbdata['FLUXERR_APER_09B'][:,aper],tbdata['FLUXERR_APER_10B'][:,aper], 
#                tbdata['FLUXERR_APER_11B'][:,aper],tbdata['FLUXERR_APER_12B'][:,aper]], axis=1)
#    return fluxerr
#
#def flux1_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,0], #tbdata['FLUX_APER_06B'][:,0],
#                tbdata['FLUX_APER_07B'][:,0], tbdata['FLUX_APER_08B'][:,0],
#                tbdata['FLUX_APER_09B'][:,0], tbdata['FLUX_APER_10B'][:,0], 
#                tbdata['FLUX_APER_11B'][:,0], tbdata['FLUX_APER_12B'][:,0]]), axis=1)
#    return flux
#
#def fluxerr1_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,0],#tbdata['FLUXERR_APER_06B'][:,0],
#                tbdata['FLUXERR_APER_07B'][:,0],tbdata['FLUXERR_APER_08B'][:,0],
#                tbdata['FLUXERR_APER_09B'][:,0],tbdata['FLUXERR_APER_10B'][:,0], 
#                tbdata['FLUXERR_APER_11B'][:,0],tbdata['FLUXERR_APER_12B'][:,0]], axis=1)
#    return fluxerr
#
#def jflux1_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,0], tbdata['FLUX_APER_06B'][:,0],
#                tbdata['FLUX_APER_07B'][:,0], tbdata['FLUX_APER_08B'][:,0],
#                tbdata['FLUX_APER_09B'][:,0], tbdata['FLUX_APER_10B'][:,0], 
#                tbdata['FLUX_APER_11B'][:,0], tbdata['FLUX_APER_12B'][:,0]], axis=1)
#    return flux
#
#def jfluxerr1_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,0],tbdata['FLUXERR_APER_06B'][:,0],
#                tbdata['FLUXERR_APER_07B'][:,0],tbdata['FLUXERR_APER_08B'][:,0],
#                tbdata['FLUXERR_APER_09B'][:,0], tbdata['FLUXERR_APER_10B'][:,0], 
#                tbdata['FLUXERR_APER_11B'][:,0], tbdata['FLUXERR_APER_12B'][:,0]], axis=1)
#    return fluxerr
#
#def flux2_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,1], #tbdata['FLUX_APER_06B'][:,1],
#                tbdata['FLUX_APER_07B'][:,1], tbdata['FLUX_APER_08B'][:,1],
#                tbdata['FLUX_APER_09B'][:,1], tbdata['FLUX_APER_10B'][:,1], 
#                tbdata['FLUX_APER_11B'][:,1], tbdata['FLUX_APER_12B'][:,1]]), axis=1)
#    return flux
#
#def fluxerr2_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,1],#tbdata['FLUXERR_APER_06B'][:,1],
#                tbdata['FLUXERR_APER_07B'][:,1],tbdata['FLUXERR_APER_08B'][:,1],
#                tbdata['FLUXERR_APER_09B'][:,1],tbdata['FLUXERR_APER_10B'][:,1], 
#                tbdata['FLUXERR_APER_11B'][:,1],tbdata['FLUXERR_APER_12B'][:,1]], axis=1)
#    return fluxerr
#
#def jflux2_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,1], tbdata['FLUX_APER_06B'][:,1],
#                tbdata['FLUX_APER_07B'][:,1], tbdata['FLUX_APER_08B'][:,1],
#                tbdata['FLUX_APER_09B'][:,1], tbdata['FLUX_APER_10B'][:,1], 
#                tbdata['FLUX_APER_11B'][:,1], tbdata['FLUX_APER_12B'][:,1]], axis=1)
#    return flux
#
#def jfluxerr2_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,1],tbdata['FLUXERR_APER_06B'][:,1],
#                tbdata['FLUXERR_APER_07B'][:,1],tbdata['FLUXERR_APER_08B'][:,1],
#                tbdata['FLUXERR_APER_09B'][:,1], tbdata['FLUXERR_APER_10B'][:,1], 
#                tbdata['FLUXERR_APER_11B'][:,1], tbdata['FLUXERR_APER_12B'][:,1]], axis=1)
#    return fluxerr
#
#def flux3_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,2],# tbdata['FLUX_APER_06B'][:,2],
#                tbdata['FLUX_APER_07B'][:,2], tbdata['FLUX_APER_08B'][:,2],
#                tbdata['FLUX_APER_09B'][:,2], tbdata['FLUX_APER_10B'][:,2], 
#                tbdata['FLUX_APER_11B'][:,2], tbdata['FLUX_APER_12B'][:,2]]), axis=1)
#    return flux
#
#def fluxerr3_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,2],#tbdata['FLUXERR_APER_06B'][:,2],
#                tbdata['FLUXERR_APER_07B'][:,2],tbdata['FLUXERR_APER_08B'][:,2],
#                tbdata['FLUXERR_APER_09B'][:,2],tbdata['FLUXERR_APER_10B'][:,2], 
#                tbdata['FLUXERR_APER_11B'][:,2],tbdata['FLUXERR_APER_12B'][:,2]], axis=1)
#    return fluxerr
#
#def jflux3_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,2], tbdata['FLUX_APER_06B'][:,2],
#                tbdata['FLUX_APER_07B'][:,2], tbdata['FLUX_APER_08B'][:,2],
#                tbdata['FLUX_APER_09B'][:,2], tbdata['FLUX_APER_10B'][:,2], 
#                tbdata['FLUX_APER_11B'][:,2], tbdata['FLUX_APER_12B'][:,2]], axis=1)
#    return flux
#
#def jfluxerr3_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,2],tbdata['FLUXERR_APER_06B'][:,3],
#                tbdata['FLUXERR_APER_07B'][:,2],tbdata['FLUXERR_APER_08B'][:,2],
#                tbdata['FLUXERR_APER_09B'][:,2], tbdata['FLUXERR_APER_10B'][:,2], 
#                tbdata['FLUXERR_APER_11B'][:,2], tbdata['FLUXERR_APER_12B'][:,2]], axis=1)
#    return fluxerr
#
#def flux4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,3], #tbdata['FLUX_APER_06B'][:,3],
#                tbdata['FLUX_APER_07B'][:,3], tbdata['FLUX_APER_08B'][:,3],
#                tbdata['FLUX_APER_09B'][:,3], tbdata['FLUX_APER_10B'][:,3], 
#                tbdata['FLUX_APER_11B'][:,3], tbdata['FLUX_APER_12B'][:,3]]), axis=1)
#    return flux
#
#def fluxerr4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 3 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,3],#tbdata['FLUXERR_APER_06B'][:,3],
#                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
#                tbdata['FLUXERR_APER_09B'][:,3],tbdata['FLUXERR_APER_10B'][:,3], 
#                tbdata['FLUXERR_APER_11B'][:,3],tbdata['FLUXERR_APER_12B'][:,3]], axis=1)
#    return fluxerr
#
#def jflux4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,3], tbdata['FLUX_APER_06B'][:,3],
#                tbdata['FLUX_APER_07B'][:,3], tbdata['FLUX_APER_08B'][:,3],
#                tbdata['FLUX_APER_09B'][:,3], tbdata['FLUX_APER_10B'][:,3], 
#                tbdata['FLUX_APER_11B'][:,3], tbdata['FLUX_APER_12B'][:,3]], axis=1)
#    return flux
#
#def jfluxerr4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,3],tbdata['FLUXERR_APER_06B'][:,3],
#                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
#                tbdata['FLUXERR_APER_09B'][:,3], tbdata['FLUXERR_APER_10B'][:,3], 
#                tbdata['FLUXERR_APER_11B'][:,3], tbdata['FLUXERR_APER_12B'][:,3]], axis=1)
#    return fluxerr
#
#def hflux4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_06B'][:,3],
#                tbdata['FLUX_APER_07B'][:,3], tbdata['FLUX_APER_08B'][:,3],
#                tbdata['FLUX_APER_09B'][:,3], tbdata['FLUX_APER_10B'][:,3], 
#                tbdata['FLUX_APER_11B'][:,3]]), axis=1)#, tbdata['FLUX_APER_12B'][:,3]
#    return flux
#
#def hfluxerr4_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_06B'][:,3],
#                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
#                tbdata['FLUXERR_APER_09B'][:,3],tbdata['FLUXERR_APER_10B'][:,3], 
#                tbdata['FLUXERR_APER_11B'][:,3]], axis=1)#,tbdata['FLUXERR_APER_12B'][:,3]
#    return fluxerr
#
#def flux5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,4], #tbdata['FLUX_APER_06B'][:,4],
#                tbdata['FLUX_APER_07B'][:,4], tbdata['FLUX_APER_08B'][:,4],
#                tbdata['FLUX_APER_09B'][:,4], tbdata['FLUX_APER_10B'][:,4], 
#                tbdata['FLUX_APER_11B'][:,4], tbdata['FLUX_APER_12B'][:,4]]), axis=1)
#    return flux
#
#def fluxerr5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,4],#tbdata['FLUXERR_APER_06B'][:,4],
#                tbdata['FLUXERR_APER_07B'][:,4],tbdata['FLUXERR_APER_08B'][:,4],
#                tbdata['FLUXERR_APER_09B'][:,4],tbdata['FLUXERR_APER_10B'][:,4], 
#                tbdata['FLUXERR_APER_11B'][:,4],tbdata['FLUXERR_APER_12B'][:,4]], axis=1)
#    return fluxerr
#
#def jflux5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,4], tbdata['FLUX_APER_06B'][:,4],
#                tbdata['FLUX_APER_07B'][:,4], tbdata['FLUX_APER_08B'][:,4],
#                tbdata['FLUX_APER_09B'][:,4], tbdata['FLUX_APER_10B'][:,4], 
#                tbdata['FLUX_APER_11B'][:,4], tbdata['FLUX_APER_12B'][:,4]], axis=1)
#    return flux
#
#def jfluxerr5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,4],tbdata['FLUXERR_APER_06B'][:,4],
#                tbdata['FLUXERR_APER_07B'][:,4],tbdata['FLUXERR_APER_08B'][:,4],
#                tbdata['FLUXERR_APER_09B'][:,4], tbdata['FLUXERR_APER_10B'][:,4], 
#                tbdata['FLUXERR_APER_11B'][:,4], tbdata['FLUXERR_APER_12B'][:,4]], axis=1)
#    return fluxerr
#
#def hflux5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_06B'][:,4], #tbdata['FLUX_APER_06B'][:,4],
#                tbdata['FLUX_APER_07B'][:,4], tbdata['FLUX_APER_08B'][:,4],
#                tbdata['FLUX_APER_09B'][:,4], tbdata['FLUX_APER_10B'][:,4], 
#                tbdata['FLUX_APER_11B'][:,4], tbdata['FLUX_APER_12B'][:,4]]), axis=1)
#    return flux
#
#def hfluxerr5_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_06B'][:,4],#tbdata['FLUXERR_APER_06B'][:,4],
#                tbdata['FLUXERR_APER_07B'][:,4],tbdata['FLUXERR_APER_08B'][:,4],
#                tbdata['FLUXERR_APER_09B'][:,4],tbdata['FLUXERR_APER_10B'][:,4], 
#                tbdata['FLUXERR_APER_11B'][:,4],tbdata['FLUXERR_APER_12B'][:,4]], axis=1)
#    return fluxerr
#def flux6_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['FLUX_APER_05B'][:,5], #tbdata['FLUX_APER_06B'][:,4],
#                tbdata['FLUX_APER_07B'][:,5], tbdata['FLUX_APER_08B'][:,5],
#                tbdata['FLUX_APER_09B'][:,5], tbdata['FLUX_APER_10B'][:,5], 
#                tbdata['FLUX_APER_11B'][:,5], tbdata['FLUX_APER_12B'][:,5]]), axis=1)
#    return flux
#
#def fluxerr6_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,5],#tbdata['FLUXERR_APER_06B'][:,4],
#                tbdata['FLUXERR_APER_07B'][:,5],tbdata['FLUXERR_APER_08B'][:,5],
#                tbdata['FLUXERR_APER_09B'][:,5],tbdata['FLUXERR_APER_10B'][:,5], 
#                tbdata['FLUXERR_APER_11B'][:,5],tbdata['FLUXERR_APER_12B'][:,5]], axis=1)
#    return fluxerr
#
#def jflux6_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture data for each 
#    epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack([tbdata['FLUX_APER_05B'][:,5], tbdata['FLUX_APER_06B'][:,5],
#                tbdata['FLUX_APER_07B'][:,5], tbdata['FLUX_APER_08B'][:,5],
#                tbdata['FLUX_APER_09B'][:,5], tbdata['FLUX_APER_10B'][:,5], 
#                tbdata['FLUX_APER_11B'][:,5], tbdata['FLUX_APER_12B'][:,5]], axis=1)
#    return flux
#
#def jfluxerr6_stacks(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,5],tbdata['FLUXERR_APER_06B'][:,5],
#                tbdata['FLUXERR_APER_07B'][:,5], tbdata['FLUXERR_APER_08B'][:,5],
#                tbdata['FLUXERR_APER_09B'][:,5], tbdata['FLUXERR_APER_10B'][:,5], 
#                tbdata['FLUXERR_APER_11B'][:,5], tbdata['FLUXERR_APER_12B'][:,5]], axis=1)
#    return fluxerr
#
#def fluxerr5_stacks_corr(tbdata):
#    ''' Function that takes a catalogue of flux data from the sextracor output
#    and makes a np array containing only the 2 arcsec aperture error data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['FLUXERR_APER_05B_CORR'],#tbdata['FLUXERR_APER_06B_CORR'],
#                tbdata['FLUXERR_APER_07B_CORR'],tbdata['FLUXERR_APER_08B_CORR'],
#                tbdata['FLUXERR_APER_09B_CORR'],tbdata['FLUXERR_APER_10B_CORR'], 
#                tbdata['FLUXERR_APER_11B_CORR'],tbdata['FLUXERR_APER_12B_CORR']], axis=1)
#    return fluxerr
#     
#def mag1_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['MAG_APER_05B'][:,0], #tbdata['MAG_APER_06B'][:,0],
#                tbdata['MAG_APER_07B'][:,0], tbdata['MAG_APER_08B'][:,0],
#                tbdata['MAG_APER_09B'][:,0], tbdata['MAG_APER_10B'][:,0], 
#                tbdata['MAG_APER_11B'][:,0], tbdata['MAG_APER_12B'][:,0]]), axis=1)
#    return flux
#
#def magerr1_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture error 
#    data for each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,0],#tbdata['MAGERR_APER_06B'][:,0],
#                tbdata['MAGERR_APER_07B'][:,0],tbdata['MAGERR_APER_08B'][:,0],
#                tbdata['MAGERR_APER_09B'][:,0],tbdata['MAGERR_APER_10B'][:,0], 
#                tbdata['MAGERR_APER_11B'][:,0],tbdata['MAGERR_APER_12B'][:,0]], axis=1)
#    return fluxerr
#    
#def mag2_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['MAG_APER_05B'][:,1], #tbdata['MAG_APER_06B'][:,1],
#                tbdata['MAG_APER_07B'][:,1], tbdata['MAG_APER_08B'][:,1],
#                tbdata['MAG_APER_09B'][:,1], tbdata['MAG_APER_10B'][:,1], 
#                tbdata['MAG_APER_11B'][:,1], tbdata['MAG_APER_12B'][:,1]]), axis=1)
#    return flux
#
#def magerr2_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture error 
#    data for each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,1],#tbdata['MAGERR_APER_06B'][:,1],
#                tbdata['MAGERR_APER_07B'][:,1],tbdata['MAGERR_APER_08B'][:,1],
#                tbdata['MAGERR_APER_09B'][:,1],tbdata['MAGERR_APER_10B'][:,1], 
#                tbdata['MAGERR_APER_11B'][:,1],tbdata['MAGERR_APER_12B'][:,1]], axis=1)
#    return fluxerr
#   
#def mag3_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['MAG_APER_05B'][:,2], #tbdata['MAG_APER_06B'][:,2],
#                tbdata['MAG_APER_07B'][:,2], tbdata['MAG_APER_08B'][:,2],
#                tbdata['MAG_APER_09B'][:,2], tbdata['MAG_APER_10B'][:,2], 
#                tbdata['MAG_APER_11B'][:,2], tbdata['MAG_APER_12B'][:,2]]), axis=1)
#    return flux
#
#def magerr3_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture error 
#    data for each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,2],#tbdata['MAGERR_APER_06B'][:,2],
#                tbdata['MAGERR_APER_07B'][:,2],tbdata['MAGERR_APER_08B'][:,2],
#                tbdata['MAGERR_APER_09B'][:,2],tbdata['MAGERR_APER_10B'][:,2], 
#                tbdata['MAGERR_APER_11B'][:,2],tbdata['MAGERR_APER_12B'][:,2]], axis=1)
#    return fluxerr
#
#def mag4_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture data for 
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['MAG_APER_05B'][:,3], #tbdata['MAG_APER_06B'][:,3],
#                tbdata['MAG_APER_07B'][:,3], tbdata['MAG_APER_08B'][:,3],
#                tbdata['MAG_APER_09B'][:,3], tbdata['MAG_APER_10B'][:,3], 
#                tbdata['MAG_APER_11B'][:,3], tbdata['MAG_APER_12B'][:,3]]), axis=1)
#    return flux
#
#def magerr4_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextracor 
#    output and makes a np array containing only the 3 arcsec aperture error 
#    data for each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,3],#tbdata['MAGERR_APER_06B'][:,3],
#                tbdata['MAGERR_APER_07B'][:,3],tbdata['MAGERR_APER_08B'][:,3],
#                tbdata['MAGERR_APER_09B'][:,3],tbdata['MAGERR_APER_10B'][:,3], 
#                tbdata['MAGERR_APER_11B'][:,3],tbdata['MAGERR_APER_12B'][:,3]], axis=1)
#    return fluxerr
#
#def mag5_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextractor 
#    output and makes a np array containing only the 2 arcsec aperture data for
#    each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        flux = an array with 8 columns containing flux values for each year '''
#        
#    flux = np.stack(([tbdata['MAG_APER_05B'][:,4],# tbdata['MAG_APER_06B'][:,4],
#                tbdata['MAG_APER_07B'][:,4], tbdata['MAG_APER_08B'][:,4],
#                tbdata['MAG_APER_09B'][:,4], tbdata['MAG_APER_10B'][:,4], 
#                tbdata['MAG_APER_11B'][:,4], tbdata['MAG_APER_12B'][:,4]]), axis=1)
#    return flux
#
#def magerr5_stacks(tbdata):
#    ''' Function that takes a catalogue of magnitude data from the sextractor 
#    output and makes a np array containing only the 2 arcsec aperture error 
#    data for each epoch
#    Input:
#        tbdata = the original combined catalogue of flux data 
#    Output:
#        fluxerr = an array with 8 columns containing flux error values for each
#        year '''
#        
#    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,4],#tbdata['MAGERR_APER_06B'][:,4],
#                tbdata['MAGERR_APER_07B'][:,4],tbdata['MAGERR_APER_08B'][:,4],
#                tbdata['MAGERR_APER_09B'][:,4],tbdata['MAGERR_APER_10B'][:,4], 
#                tbdata['MAGERR_APER_11B'][:,4],tbdata['MAGERR_APER_12B'][:,4]], axis=1)
#    return fluxerr 

#def create_error_array(flux, sigtb, tbdata):
#    ''' Function that creates an error array from sigma values calculated from
#    data variations within flux bins '''
#    bins = np.array(sigtb.colnames).astype(int)
#    avgflux = np.nanmean(flux, axis=1)
#    errarr = np.empty(np.shape(flux))
#    for n, lower in enumerate(bins):
#        if n == 0:
#            #remove lower values
#            mask = avgflux< lower
#            nanarr = np.full_like(sigtb[str(lower)], np.nan)
#            flux[mask,:] = nanarr
#            errarr[mask,:] = nanarr
#        mask1 = avgflux>lower 
#        if n != len(bins)-1:
#            mask2 = avgflux<bins[n+1]
#        else:
#            #remove lower values
#            mask = avgflux > 6309573 #max bin
#            nanarr = np.full_like(sigtb[str(lower)], np.nan)
#            flux[mask,:] = nanarr
#            errarr[mask,:] = nanarr
#            mask2 = avgflux < 6309573 #max bin
#        
#        mask = mask1*mask2.astype(bool)
#        errarr[mask,:] = sigtb[str(lower)]
#        
#    # remove nans
#    mask = ~np.isnan(flux).any(axis=1)
#    flux = flux[mask]
#    errarr = errarr[mask]
#    tbdata = tbdata[mask]
#    return flux, errarr, tbdata
#
#def create_quad_error_array(sigtb, tbdata, aper=5, quadoutput=False):
#    ''' Function that creates an error array from sigma values calculated from
#    data variations within flux bins '''
#    bins = np.array(sigtb.colnames)
#    binarr = np.empty([4,int(len(bins)/4)])
#    k = 0
#    l = 0
#    m = 0
#    n = 0
#    ### Create binedge array for each quadrant ###
#    for bin in bins:
#        if bin[0] == '1':
#            binarr[0,k] = int(bin[2:])
#            k+=1
#        elif bin[0] == '2':
#            binarr[1,l] = int(bin[2:])
#            l+=1
#        elif bin[0] == '3':
#            binarr[2,m] = int(bin[2:])
#            m+=1
#        elif bin[0] == '4':
#            binarr[3,n] = int(bin[2:])
#            n+=1
#    
#    ### Set up empty arrays for data ###
#    flux = np.array([])
#    errarr = np.array([])
#    newtbdata = []
#    quadflux = {}
#    quaderr = {}
#    newquaddata = {}
#    
#    ### Get quadrant data ###
#    quaddata = quadrants(tbdata, '05B')
#    for n, qdata in enumerate(quaddata):
#        ### create flux stacks and find average
#        qflux = k_mag_flux.flux_stacks(qdata, aper)
##        if aper == 6:
##            qflux = flux6_stacks(qdata)
##        elif aper == 5:
##            qflux = flux5_stacks(qdata)
##        elif aper == 4:
##            qflux = flux4_stacks(qdata)
##        elif aper == 3:
##            qflux = flux3_stacks(qdata)
##        elif aper == 2:
##            qflux = flux2_stacks(qdata)
##        elif aper == 1:
##            qflux = flux1_stacks(qdata)
#        avgflux = np.nanmean(qflux, axis=1)
#        qerrarr = np.empty(np.shape(qflux))
#        
#        ### Find values within bins and assign correct sigma ###
#        for m, lower in enumerate(binarr[n,0:-1]):
#            if m == 0:
#                #remove lower values
#                mask = avgflux< lower
#                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
#                qflux[mask,:] = nanarr
#                qerrarr[mask,:] = nanarr
#            mask1 = avgflux>lower 
#            if m != len(binarr)-1:
#                mask2 = avgflux<binarr[n,m+1]
#            else:
#                #remove lower values
#                mask = avgflux > 6309573 #max bin
#                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
#                qflux[mask,:] = nanarr
#                qerrarr[mask,:] = nanarr
#                mask2 = avgflux < 6309573 #max bin
#            
#            mask = mask1*mask2.astype(bool)
#            qerrarr[mask,:] = sigtb[str(n+1)+' '+str(int(lower))]
#            
#        # remove nans
#        mask = ~np.isnan(qflux).any(axis=1)
#        qflux = qflux[mask]
#        qerrarr = qerrarr[mask]
#        qdata = qdata[mask]
#        
#        if quadoutput == False:
#            ### Define full arrays to use in the rest of the analysis ###
#            if n == 0:
#                flux = np.copy(qflux)
#                errarr = np.copy(qerrarr)
#                newtbdata = np.copy(qdata)
#            else:
#                flux = np.vstack((flux, np.copy(qflux)))
#                errarr = np.vstack((errarr, np.copy(qerrarr)))
#                newtbdata = np.hstack((newtbdata, np.copy(qdata)))
#        
#        else: #this means wants quadrant output
#            quadflux[n] = qflux
#            quaderr[n] = qerrarr
#            newquaddata[n] = qdata
#    
#    if quadoutput == False:
#        if np.isin('X-ray', newtbdata.dtype.names):
#            newtbdata['X-ray'][newtbdata['X-ray']==70] = False 
#            newtbdata['X-ray'][newtbdata['X-ray']==84] = True
#            newtbdata['X-ray'] = newtbdata['X-ray'].astype(bool)
#        return flux, errarr, newtbdata
#    else:
#        return quadflux, quaderr, newquaddata
#
#def create_quad_error_array_J(sigtb, tbdata, aper=5, quadoutput=False):
#    ''' Function that creates an error array from sigma values calculated from
#    data variations within flux bins '''
#    bins = np.array(sigtb.colnames)
#    binarr = np.empty([4,int(len(bins)/4)])
#    k = 0
#    l = 0
#    m = 0
#    n = 0
#    ### Create binedge array for each quadrant ###
#    for bin in bins:
#        if bin[0] == '1':
#            binarr[0,k] = int(bin[2:])
#            k+=1
#        elif bin[0] == '2':
#            binarr[1,l] = int(bin[2:])
#            l+=1
#        elif bin[0] == '3':
#            binarr[2,m] = int(bin[2:])
#            m+=1
#        elif bin[0] == '4':
#            binarr[3,n] = int(bin[2:])
#            n+=1
#    
#    ### Set up empty arrays for data ###
#    flux = np.array([])
#    errarr = np.array([])
#    newtbdata = []
#    quadflux = {}
#    quaderr = {}
#    newquaddata = {}
#    
#    ### Get quadrant data ###
#    quaddata = quadrants(tbdata, '05B')
#    for n, qdata in enumerate(quaddata):
#        ### create flux stacks and find average
#        if aper == 6:
#            qflux = jflux6_stacks(qdata)
#        elif aper == 5:
#            qflux = jflux5_stacks(qdata)
#        elif aper == 4:
#            qflux = jflux4_stacks(qdata)
#        elif aper == 3:
#            qflux = jflux3_stacks(qdata)
#        elif aper == 2:
#            qflux = jflux2_stacks(qdata)
#        elif aper == 1:
#            qflux = jflux1_stacks(qdata)
#        avgflux = np.nanmean(qflux, axis=1)
#        qerrarr = np.empty(np.shape(qflux))
#        
#        ### Find values within bins and assign correct sigma ###
#        for m, lower in enumerate(binarr[n,0:-1]):
#            if m == 0:
#                #remove lower values
#                mask = avgflux< lower
#                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
#                qflux[mask,:] = nanarr
#                qerrarr[mask,:] = nanarr
#            mask1 = avgflux>lower 
#            if m != len(binarr)-1:
#                mask2 = avgflux<binarr[n,m+1]
#            else:
#                #remove lower values
#                mask = avgflux > 6309573 #max bin
#                nanarr = np.full_like(sigtb[str(n+1)+' '+str(int(lower))], np.nan)
#                qflux[mask,:] = nanarr
#                qerrarr[mask,:] = nanarr
#                mask2 = avgflux < 6309573 #max bin
#            
#            mask = mask1*mask2.astype(bool)
#            qerrarr[mask,:] = sigtb[str(n+1)+' '+str(int(lower))]
#            
#        # remove nans
#        mask = ~np.isnan(qflux).any(axis=1)
#        qflux = qflux[mask]
#        qerrarr = qerrarr[mask]
#        qdata = qdata[mask]
#        
#        if quadoutput == False:
#            ### Define full arrays to use in the rest of the analysis ###
#            if n == 0:
#                flux = np.copy(qflux)
#                errarr = np.copy(qerrarr)
#                newtbdata = np.copy(qdata)
#            else:
#                flux = np.vstack((flux, np.copy(qflux)))
#                errarr = np.vstack((errarr, np.copy(qerrarr)))
#                newtbdata = np.hstack((newtbdata, np.copy(qdata)))
#        
#        else: #this means wants quadrant output
#            quadflux[n] = qflux
#            quaderr[n] = qerrarr
#            newquaddata[n] = qdata
#    
#    if quadoutput == False:
#        if np.isin('X-ray', newtbdata.dtype.names):
#            newtbdata['X-ray'][newtbdata['X-ray']==70] = False 
#            newtbdata['X-ray'][newtbdata['X-ray']==84] = True
#            newtbdata['X-ray'] = newtbdata['X-ray'].astype(bool)
#        return flux, errarr, newtbdata
#    else:
#        return quadflux, quaderr, newquaddata

#def lightcurve4(ob, fitsdata)  :
#    ''' Function that plots the light curve of an object in terms of its flux 
#    in an aperture 4 pixels across (i.e. 2 arcsec in diameter) 
#    Inputs:
#        ob = the ID of the object that you want the lightcurve from
#        fitsdata = the original catalogue of data that the curve will be 
#                    plotted from 
#    Output:
#        None '''
#        
#    #Get data for the object called
#    mask = fitsdata['NUMBER_05B'] == ob
#    obdata = fitsdata[mask]
#    if not obdata: #Reject if no object number matches the input value
#        print('error- invalid object number')
#        return
#    
#    #Create arrays of flux values and error values
#    flux = np.array([obdata['FLUX_APER_05B'][:,4],obdata['FLUX_APER_06B'][:,4],
#                    obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
#                    obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
#                    obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
#    
#    fluxerr = np.array([obdata['FLUXERR_APER_05B'][:,4],obdata['FLUXERR_APER_06B'][:,4],
#                    obdata['FLUXERR_APER_07B'][:,4],obdata['FLUXERR_APER_08B'][:,4],
#                    obdata['FLUXERR_APER_09B'][:,4],obdata['FLUXERR_APER_10B'][:,4], 
#                    obdata['FLUXERR_APER_11B'][:,4],obdata['FLUXERR_APER_12B'][:,4]])
#    
#    # normalise and correct for seeing
#    flux = psf_correct(flux, flux, 'mean')
#    fluxerr = psf_correct(flux, fluxerr, 'mean')
#    avgflux =np.mean(flux)
#    normflux = flux/avgflux
#    normerr = fluxerr/avgflux
#    
#    #set up time variable for plot
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    
#    #Plot graph in new figure
#    plt.figure()
#    plt.xticks(t, years)
#    plt.errorbar(t, normflux, yerr=normerr, fmt = 'ro')
#    plt.xlabel('Semester')
#    plt.ylabel('Flux of object')
#    plt.title('Light curve for object number %i' % ob)
#    return
#
#def lightcurve5(ob, fitsdata)  :
#    ''' Function that plots the light curve of an object in terms of its flux 
#    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
#    Inputs:
#        ob = the ID of the object that you want the lightcurve from
#        fitsdata = the original catalogue of data that the curve will be 
#                    plotted from 
#    Output:
#        None '''
#        
#    #Get data for the object called
#    mask = fitsdata['NUMBER_05B'] == ob
#    obdata = fitsdata[mask]
#    if not obdata: #Reject if no object number matches the input value
#        print('error- '+str(ob)+' is an invalid object number')
#        return
#    
#    #Create arrays of flux values and error values
#    flux = np.array([obdata['MAG_APER_05B'][:,4],obdata['MAG_APER_06B'][:,4],
#                    obdata['MAG_APER_07B'][:,4],obdata['MAG_APER_08B'][:,4],
#                    obdata['MAG_APER_09B'][:,4],obdata['MAG_APER_10B'][:,4], 
#                    obdata['MAG_APER_11B'][:,4],obdata['MAG_APER_12B'][:,4]])
#    
#    fluxerr = np.array([obdata['MAGERR_APER_05B'][:,4],obdata['MAGERR_APER_06B'][:,4],
#                    obdata['MAGERR_APER_07B'][:,4],obdata['MAGERR_APER_08B'][:,4],
#                    obdata['MAGERR_APER_09B'][:,4],obdata['MAGERR_APER_10B'][:,4], 
#                    obdata['MAGERR_APER_11B'][:,4],obdata['MAGERR_APER_12B'][:,4]])
#
#    
#    #set up time variable for plot
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    
#    #Plot graph in new figure
#    plt.figure()
#    plt.xticks(t, years)
#    plt.errorbar(t, flux, yerr=fluxerr, fmt = 'ro')
#    plt.xlabel('Semester')
#    plt.ylabel('K-band magnitude of object')
#    plt.title('Light curve for object number %i' % ob)
#    return
#
#def lightcurveflux5(ob, fitsdata, corrected=False, new_fig=True) :
#    ''' Function that plots the light curve of an object in terms of its flux 
#    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
#    Inputs:
#        ob = the ID of the object that you want the lightcurve from
#        fitsdata = the original catalogue of data that the curve will be 
#                    plotted from 
#    Output:
#        None '''
#    sigtb = Table.read('quad_epoch_sigma_table_extra_clean.fits')
#
#    #Get data for the object called
#    mask = fitsdata['NUMBER_05B'] == ob
#    obdata = fitsdata[mask]
#    if not obdata: #Reject if no object number matches the input value
#        print('error- invalid object number')
#        return
#    
#    #Create arrays of flux values and error values
#    flux = np.array([obdata['FLUX_APER_05B'][:,4],#obdata['FLUX_APER_06B'][:,4],
#                    obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
#                    obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
#                    obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
#    flux, fluxerr, obdata = create_quad_error_array(sigtb, obdata)
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
#    
#    #set up time variable for plot
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    x = [1,3,4,5,6,7,8]
#    print('Plotting lightcurve')
#    
#    #Plot graph in new figure
#    if new_fig == True:
#        plt.figure()
#    plt.xticks(t, years)
#    plt.errorbar(x, flux, yerr=fluxerr, fmt = 'ro')
#    plt.xlabel('Semester')
#    plt.ylabel('K-band flux of object')
#    plt.title('Light curve for object number %i' % ob)
#    return
#
#def chandra_only(tbdata):
#    ''' Function that restricts the objects included in analysis to only 
#    those within the Chandra footprint 
#    Input:
#        tbdata = original catalogue of data 
#    Output:
#        newtbdata = new catalogue of data which only includes objects within 
#                    the chandra footprint '''
#    ### Restrict objects to those in the Chandra field ###
#    mask1 = tbdata['DELTA_J2000_05B'] < -4.93 #max Dec
#    mask2 = tbdata['DELTA_J2000_05B'] > -5.403 #min Dec
#    mask3 = tbdata['ALPHA_J2000_05B'] < 34.72 #max RA
#    mask4 = tbdata['ALPHA_J2000_05B'] > 34.07 #min RA
#    mask = mask1 * mask2 * mask3 * mask4
#    newtbdata = tbdata[mask]
#    return newtbdata
    
##Calculate varience for each row
#def sigmasq(flux, baseerr):
#    ''' Function that calculates the excess varience value for each row in an 
#    array 
#    Inputs:
#        flux = array of fluxes from objects in a number of epochs 
#        baseerr = array of errors that the mean error should be calculated from
#    Output:
#        sig = array of excess variance values for every object '''
#    avgflux = np.nanmean(flux, axis=1)
##    meanerr = np.nanmean(baseerr, axis=1)
##    meanerrsq = np.square(meanerr)
#    N = np.size(flux, axis=1)
#    numobs = np.size(flux, axis=0)
#    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
#    sigsum = np.nansum(sig, axis=1)
#    normsig = sigsum/(N)
#    return normsig
#
#def normsigmasq(flux, baseerr):
#    ''' Function that calculates the excess varience value for each row in an 
#    array 
#    Inputs:
#        flux = array of fluxes from objects in a number of epochs 
#        baseerr = array of errors that the mean error should be calculated from
#    Output:
#        sig = array of excess variance values for every object '''
#    if np.shape(flux) == (0,8):
#        return np.array([])
#    avgflux = np.nanmean(flux, axis=1)
#    N = np.size(flux, axis=1)
#    numobs = np.size(flux, axis=0)
#    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
#    sigsum = np.nansum(sig, axis=1)
#    normsig = sigsum/(N)
#    return normsig
#
#def maximum_likelihood_fig(testmag, testmagerr, meanmag, posvar):
#
#    # Calculate likelihood curve
#    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
#            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
#            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
#    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
##    Lstd = np.std(L)
##    idx = (np.abs(L+Lstd).argmin())
##    sigstd = posvar[idx]
##    err = np.abs(sig+sigstd)
#    err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
#    plt.figure()
#    plt.plot(posvar, L)
#    plt.vlines(sig, np.min(L), np.max(L))
#    plt.vlines(sig+err, np.min(L), np.max(L))
#    return sig, err
#
#def maximum_likelihood(testmag, testmagerr, meanmag, posvar, n=None, printn=10):
#    if n != None and n%printn == 0:
#        print(n)
#    # Calculate likelihood curve
#    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
#            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
#            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
#    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
#    if np.sum(L) == 0:
#        return sig, np.nan
#    else:
#        err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
#        return sig, err


#
#def fluxbin(min, max, flux, tbdata):
#    ''' Separate a flux array into a bin depending on the average flux of the
#    object 
#    Inputs:
#        min = the minimum flux value for the bin
#        max = the maximum flux value for the bin
#        flux = the array for fluxes for the object
#    Output:
#        fluxbin = the array of fluxes for objects whose average flux is within 
#                    the limits of the bin '''
#    avgflux = np.mean(flux, axis=1)
#    bina = avgflux >= min
#    binb = avgflux < max
#    bin = bina*binb
#    fluxbin = flux[bin]
#    tbdata = tbdata[bin]
#    return fluxbin, tbdata
#
#def fluxbinerr(min, max, flux, fluxerr):
#    ''' Separate a flux array into a bin depending on the average flux of the
#    object 
#    Inputs:
#        min = the minimum flux value for the bin
#        max = the maximum flux value for the bin
#        flux = the array for fluxes for the object
#        fluxerr = the array for flux errors for the object
#    Output:
#        fluxbin = the array of fluxes for objects whose average flux is within 
#                    the limits of the bin 
#        fluxerrbin = the array of flux errors for objects whose average flux is 
#                    within the limits of the bin '''
#    avgflux = np.mean(flux, axis=1)
#    bina = avgflux >= min
#    binb = avgflux < max
#    bin = bina*binb
#    fluxbin = flux[bin]
#    fluxerrbin = fluxerr[bin]
#    return fluxbin, fluxerrbin
#
#def avg_lightcurve(avgfluxarr, errarr=[], shape='o', size=None, label=None):
#    ''' Plot light curves for a provided set of fluxes rather than extracting
#    a certain objects fluxes from a larger array.
#    Input:
#        avgfluxarr = array of 8 semester flux values
#    Output:
#        ax = axes handle so the plot can be added to if necessary '''
#        
#    #Check the array has values
#    if math.isnan(avgfluxarr[0]):
#        print('Array is empty')
#        return
#    #set up time variable for plot
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    
#    #Plot graph in new figure
##    plt.figure()
#    ax = plt.axes()
#    plt.xticks(t, years)
#    if errarr == []:
#        ax.plot(t, avgfluxarr, shape, markersize=size, label=label)
#    else:
#        ax.errorbar(t, avgfluxarr, errarr, fmt=shape, label=label)
#    plt.xlabel('Semester')
#    plt.ylabel('Flux of object')
#    plt.title('Average light curve')
##    plt.ylim(ymin = 6.3, ymax=7.2)
#    return ax
#
#def normalise_flux(flux):
#    ''' Normalise each objects flux to its average value
#    Input:
#        flux = array of object flux values 
#    Output:
#        array of object flux values normalised to the average flux of the 
#        object '''
#    avgflux = np.mean(flux, axis=1)
#    return flux / avgflux[:,None]
#
#def normalise_flux_and_errors(flux, fluxerr):
#    ''' Normalise each objects flux to its average value
#    Input:
#        flux = array of object flux values 
#    Output:
#        array of object flux values normalised to the average flux of the 
#        object '''
#    avgflux = np.mean(flux, axis=1)
#    flux = flux/avgflux[:,None]
#    fluxerr = fluxerr/avgflux[:,None]
#    return flux, fluxerr
#
#def normalise_mag(mag):
#    ''' Normalise each objects flux to its average value
#    Input:
#        flux = array of object flux values 
#    Output:
#        array of object flux values normalised to the average flux of the 
#        object '''
#    avgflux = np.mean(mag, axis=1)
#    diff = avgflux - 1
#    return mag - diff[:,None]
#
#def no99(fluxn, tbdata):
#    fluxn[fluxn == 99] = np.nan
#    mask = ~np.isnan(fluxn).any(axis=1)
#    fluxn = fluxn[mask]
#    tbdata = tbdata[mask]
#    return fluxn, tbdata
#
#def noneg(fluxn, tbdata):
#    fluxn[fluxn <= 0] = np.nan
#    mask = ~np.isnan(fluxn).any(axis=1)
#    fluxn = fluxn[mask]
#    tbdata = tbdata[mask]
#    return fluxn, tbdata
#
#def fluxlim(fluxn, tbdata, lim=3527):
#    fluxn, tbdata = noneg(fluxn, tbdata) #remove negative values first
#    avgflux = np.mean(fluxn, axis=1) #find average so only remove those whose average is less than limit
#    fluxn[avgflux <= lim] = np.nan
#    mask = ~np.isnan(fluxn).any(axis=1)
#    fluxn = fluxn[mask]
#    tbdata = tbdata[mask]
#    return fluxn, tbdata
#
#def semfluxlim(fluxn, tbdata):
##    fluxn, tbdata = noneg(fluxn, tbdata) #remove negative values first
#    limits = [2920.74, 6917.62, 4023.65, 5613.71, 
#              1985.52, 2725.65, 2111.91, 1915.96]
#    for n, lim in enumerate(limits):
#        fluxn[:,n][fluxn[:,n] < lim] = np.nan        
#    mask = ~np.isnan(fluxn).any(axis=1)
#    fluxn = fluxn[mask]
#    tbdata = tbdata[mask]
#    return fluxn, tbdata
#
#def onpickchanonly(event):
#    ''' Function that plots the lightcurve of an object when it is clicked on 
#    the vairiability v flux plot '''
#        
#    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
#    tbdata = combined[1].data
#    tbdata = chandra_only(tbdata)
#    
#    ### remove values that are +/-99 ###
#    fluxn = mag5_stacks(tbdata)
#    fluxn[fluxn == 99] = np.nan
#    mask = ~np.isnan(fluxn).any(axis=1)
#    tbdata = tbdata[mask]
#    
#    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
#    if len(ob) > 1: #Reject selection if more than one object has been selected
#        print('Too many objects selected')
#        return
#    print('Object identified')
#    lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function
#
#def onpick(event):
#    ''' Function that plots the lightcurve of an object when it is clicked on 
#    the vairiability v flux plot '''
#        
#    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
#    tbdata = combined[1].data
#    
#    ### remove values that are +/-99 ###
#    fluxn = mag5_stacks(tbdata)
#    fluxn, tbdata = no99(fluxn, tbdata)
#    
#    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
#    if len(ob) > 1: #Reject selection if more than one object has been selected
#        print('Too many objects selected')
#        return
#    print('Object identified')
#    lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function
# 
#def onpickflux(event):
#    ''' Function that plots the lightcurve of an object when it is clicked on 
#    the vairiability v flux plot '''
#    print('Click registered')
#    tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
#    sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')
#    tbdata = remove_edges(tbdata)
#
#    ### remove values that are +/-99 ###
#    fluxn = flux5_stacks(tbdata)
#    fluxn, tbdata = noneg(fluxn, tbdata)
#    flux, fluxerr, tbdata = create_quad_error_array(sigtb, tbdata)
#    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
#
#    ### reset X-ray column as messed up by stacking ###
#    tbdata['X-ray'][tbdata['X-ray']==70] = False 
#    tbdata['X-ray'][tbdata['X-ray']==84] = True
#    
#    if len(ob) > 1: #Reject selection if more than one object has been selected
#        print('Too many objects selected')
#        return
#    print('Object identified')
#    
##    flux,fluxerr = normalise_flux_and_errors(flux, fluxerr)
#    
#    print('Finding object values')
#    obflux = np.reshape(flux[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
#    obfluxerr = np.reshape(fluxerr[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
#    obdata = tbdata[tbdata['NUMBER_05B']==ob]
#    print(obfluxerr)
#    
#    #set up time variable for plot
#    print('setting up plot')
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    x = [1,3,4,5,6,7,8]
##    chisq =my_chisquare_err(obflux, obfluxerr)
#    
#    plt.figure()
#    if obdata['X-ray'] == True:
#        print('Plotting x-ray')
#        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='r')
#    else:
#        print('Plotting non x-ray')
#        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='b')
##        print('Plotted non x-ray')
#    plt.xlabel('Semester')
#    plt.ylabel('Flux')
#    plt.title('Lightcurve of Object '+str(obdata['NUMBER_05B']))#+' '+r' $\chi^{2} = $'+str(round(chisq, 2)))
#    plt.xticks(t, years)
#    plt.tight_layout()
#    print('I got here')
##    lightcurveflux5(ob, tbdata) #Plot the lightcurve from the lightcurve function   
# 
#def onpickflux_2arcsec(event):
#    ''' Function that plots the lightcurve of an object when it is clicked on 
#    the vairiability v flux plot '''
#    print('Click registered')
#    tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
#    sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
#    tbdata = remove_edges(tbdata)
#
#    ### remove values that are +/-99 ###
#    fluxn = flux4_stacks(tbdata)
#    fluxn, tbdata = noneg(fluxn, tbdata)
##    fluxn, tbdata = remove_low_flux(fluxn, tbdata)
#    flux, fluxerr, tbdata = create_quad_error_array(sigtb, tbdata, aper=4)
#    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
#
#    ### reset X-ray column as messed up by stacking ###
#    tbdata['X-ray'][tbdata['X-ray']==70] = False 
#    tbdata['X-ray'][tbdata['X-ray']==84] = True
#    
#    if len(ob) > 1: #Reject selection if more than one object has been selected
#        print('Too many objects selected')
#        return
#    print('Object identified')
#    
##    flux,fluxerr = normalise_flux_and_errors(flux, fluxerr)
#    
#    print('Finding object values')
#    obflux = np.reshape(flux[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
#    obfluxerr = np.reshape(fluxerr[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
#    obdata = tbdata[tbdata['NUMBER_05B']==ob]
#    print(obfluxerr)
#    
#    #set up time variable for plot
#    print('setting up plot')
#    t = np.linspace(1, 8, num=8)
#    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#    x = [1,3,4,5,6,7,8]
##    chisq =my_chisquare_err(obflux, obfluxerr)
#    
#    plt.figure()
#    if obdata['X-ray'] == True:
#        print('Plotting x-ray')
#        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='r')
#    else:
#        print('Plotting non x-ray')
#        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='b')
##        print('Plotted non x-ray')
#    plt.xlabel('Semester')
#    plt.ylabel('Flux')
#    plt.title('Lightcurve of Object '+str(obdata['NUMBER_05B']))#+' '+r' $\chi^{2} = $'+str(round(chisq, 2)))
#    plt.xticks(t, years)
#    plt.tight_layout()
#    print('I got here')
##    lightcurveflux5(ob, tbdata) #Plot the lightcurve from the lightcurve function   

#
#def my_chisquare_err(flux, fluxerr):
##    flux, fluxerr = normalise_flux_and_errors(flux, fluxerr)
#    meanflux = np.nanmean(flux, axis=1)
#    top = np.square(flux-meanflux[:,None])
#    bot = np.square(fluxerr)
#    chi = np.nansum(top/bot, axis=1)
#    return chi
#
#def flux_variability_plot(flux, fluxchan, plottype, flux2 = [], fluxchan2 = [],
#                          fluxerr = [], fluxerr2 = [], starflux=[], starfluxerr=[],
#                          comparison = False, normalised = False, stars=False,
#                          psfcorrect=False, chanerr = [], chanerr2 = [], scale=''):
#    ''' Function to plot the variability vs mean flux plot using either the
#    MAD statistic or the Excess Variance statistic 
#    If using MAD, specify plottype = 'mad' and supply just the flux and chandra 
#    flux arrays
#    If using Excess Variance, specify plottype = 'excess' and supply the error
#    arrays as well as the flus ones
#    If you with to have a comparison plot of either type make sure to specify 
#    the greyed out one first and then the coloured one 
#    Inputs:
#        flux = array of flux values for UDS objects
#        fluxchan = array of flux values for chandra objects
#        plottype = either 'mad' or 'excess' which defines what statistic is
#                    used to quantify variability
#        flux2 = optional second array of UDS objects if comparing two sets 
#                (e.g. a before and after) this array will be the top coloured
#                layer.
#        fluxchan2 = optional second array of chandra objects 
#        fluxerr = optional array of flux errors for objects, used if 
#                    'excess' is specified
#        fluxerr2 = optional second array of flux errors for objects, used 
#                    if 'excess' and 'comparision = True' are specified
#        starflux = optional array of fluxes for stars
#        starfluxerr = optional array of flux errors for objects, used if 
#                        'excess' and 'stars=True' are specified
#        comparison = True or False (default) depending on if a comparison plot
#                        is required 
#        normalised = True or False (default) depending if the fluxes should be
#                        normalised to the objects average flux 
#        stars = True or False (default) depending on if stars should be added
#                to the plot
#    Output:
#        fig = figure handle to allow clicking for light curves to be enabled if
#                required '''
#    
#    fig = plt.figure(figsize=[8,8])
#    avgfluxperob = np.nanmean(flux, axis=1) #for UDS
#    avgfluxchanperob = np.nanmean(fluxchan, axis=1) #for non-stellar chandra
#    if stars == True:
#        savgfluxperob = np.nanmean(starflux, axis=1) #for stars
#
#    ### Check if normalisation is true and normalise if necessary ###
#    if normalised == True:
#        if plottype == 'mad':
#    #        flux = normalise_mag(flux)
#    #        fluxchan = normalise_mag(fluxchan) 
#    #        if stars == True:
#    #            starflux = normalise_mag(starflux)
#            flux = normalise_flux(flux)
#            fluxchan = normalise_flux(fluxchan) 
#            if stars == True:
#                starflux = normalise_flux(starflux)
#        else:
#            flux, fluxerr = normalise_flux_and_errors(flux, fluxerr)
#            fluxchan, chanerr = normalise_flux_and_errors(fluxchan, chanerr)
#            if stars == True:
#                starflux, starfluxerr = normalise_flux_and_errors(starflux, starfluxerr)
#    if psfcorrect == True:
#        fluxchan = psf_correct_mag(flux, fluxchan, 'median')
#        if stars == True:
#            starflux = psf_correct_mag(flux, starflux, 'median')
#        flux = psf_correct_mag(flux, flux, 'median')
#    ### Find out which plot type is specified and calculate appropriate statistic ###
#    if plottype == 'mad':
#        vary = median_absolute_deviation(flux, axis=1)
#        varychan = median_absolute_deviation(fluxchan, axis=1)
#        if stars == True:
#            varystar = median_absolute_deviation(starflux, axis=1)
#        plt.ylabel('MAD')
#    elif plottype == 'excess':
#        vary = normsigmasq(flux, fluxerr)
#        varychan = normsigmasq(fluxchan, chanerr)
#        if stars == True:
#            varystar = normsigmasq(starflux, starfluxerr)
#        plt.ylabel('Excess Variance')
#    elif plottype == 'chisq':
#        vary = my_chisquare_err(flux, fluxerr)
#        varychan = my_chisquare_err(fluxchan, chanerr)
#        if stars == True:
#            varystar = my_chisquare_err(starflux, starfluxerr)
#        plt.ylabel('Chi Squared')
#    elif plottype == 'var':
#        vary = np.var(flux, axis=1, ddof=1)
#        varychan = np.var(fluxchan, axis=1, ddof=1)
#        if stars == True:
#            varystar = np.var(starflux, axis=1, ddof=1)
#        plt.ylabel('Variance')
#    else:
#        print('Invalid plottype') #returns if unrecognised value is entered
#        return
#    
#    ### Plot the variability v mean as appropriate ###
#    if comparison == True:     
#        
#        avgfluxperob2 = np.nanmean(flux2, axis=1) #for UDS
#        avgfluxchanperob2 = np.nanmean(fluxchan2, axis=1) #for non-stellar chandra
#
#        if plottype == 'mad':
#            if normalised == True:
#    #            flux2 = normalise_mag(flux2)
#    #            fluxchan2 = normalise_mag(fluxchan2)     
#                flux2 = normalise_flux(flux2)
#                fluxchan2 = normalise_flux(fluxchan2) 
#            varycorr = median_absolute_deviation(flux2, axis=1)
#            varychancorr = median_absolute_deviation(fluxchan2, axis=1)
#        elif plottype == 'excess':
#            if normalised == True:
#                flux2, fluxerr2 = normalise_flux_and_errors(flux2, fluxerr2)
#                fluxchan2, chanerr2 = normalise_flux_and_errors(fluxchan2, chanerr2)
#            varycorr = normsigmasq(flux2, fluxerr2)
#            varychancorr = normsigmasq(fluxchan2, chanerr2)
#
#        ### plot varibility v flux graph for original in gray ###
#        plt.plot(avgfluxperob, vary, '+', color =  'tab:gray', label='Galaxy', alpha = 0.5) 
#        plt.plot(avgfluxchanperob, varychan, 'o', color =  'tab:gray',  mfc = 'none', markersize = 10,
#                 label='X-ray detected', alpha = 0.5) #no picker as will be selected in the UDS point
#        
#        ### plot varibility v flux graph for new in colour ###
#        line, = plt.plot(avgfluxperob2, varycorr, 'b+', label='Galaxy Convolved', picker=2, alpha = 0.5) #label tells 
#        # the legend what to write, picker is the readius of pixels which a user must 
#        # click within to call the object
#        plt.plot(avgfluxchanperob2, varychancorr, 'ro', mfc = 'none', markersize = 10,
#                 label='X-ray detected', alpha = 0.5) #no picker as will be selected in the UDS point
#    else:
#        if stars==True:
#            plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
#                     label='DR11 Star') 
#        line, = plt.plot(avgfluxperob, vary, 'b+', label='Galaxy', picker=2)
#        plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
#                 label='X-ray detected') #no picker as will be selected in the UDS point
#
#        
#    ### Apply required plot charateristics ###
#    plt.xscale('log')
#    if scale == 'log':
#        plt.yscale('log')
#    elif scale == 'symlog':
#        plt.yscale('symlog', linthreshy=0.0001)
#        
#    plt.ylim(3e-2,3e4)
#    plt.xlim(8e1, 1e7)
##    plt.xlim(13,26)
#    plt.xlabel('Mean Flux')
##    plt.xlabel('Mean Magnitude')
#    plt.legend()
##    plt.gca().invert_xaxis()
#    return fig, vary


#def psf_correct_flux(baseflux, initflux, avgtype):
#    ''' Function that applies a fractional PSF correction to initflux, based on
#    the constant required to correct an epochs average flux to the overall 
#    average flux in the baseflux array.
#    i.e. if doing a correction based on all objects for objects in the chandra
#    data then baseflux = flux and initflux = fluxchan.
#    If doing a correction based on the stellar fluxes for all objects in the 
#    field then baseflux = sflux and initflux = flux. 
#    Basetype must be defined as either 'all' or 'star' so that the mean value
#    can be used for all objects and the median value can be used for stellar 
#    objects 
#    Inputs:
#        baseflux = flux array to base the corrections on
#        initflux = flux array to be corrected
#        basetype = either 'mean' or 'median' which dictate if the mean or 
#                    median value is used for the correction 
#    Output:
#        Flux array with values crudely corrected for differences in seeing
#        (average flux should now be the same for each epoch). '''
#    if avgtype == 'mean':
#        avgfluxperepoch = np.mean(baseflux, axis=0)#for UDS
#        avgflux = np.mean(baseflux)
#        const = avgflux/avgfluxperepoch
#    elif avgtype == 'median':
#        avgfluxperepoch = np.median(baseflux, axis=0)#for UDS
#        avgflux = np.median(baseflux)
#        const = avgflux/avgfluxperepoch
#    else:
#        print('Invalid basetype')
#        return
#    return initflux * const[None,:]
#
#def psf_correct_mag(basemag, initmag, avgtype):
#    ''' Function that applies a fractional PSF correction to initflux, based on
#    the constant required to correct an epochs average flux to the overall 
#    average flux in the baseflux array.
#    i.e. if doing a correction based on all objects for objects in the chandra
#    data then baseflux = flux and initflux = fluxchan.
#    If doing a correction based on the stellar fluxes for all objects in the 
#    field then baseflux = sflux and initflux = flux. 
#    Basetype must be defined as either 'all' or 'star' so that the mean value
#    can be used for all objects and the median value can be used for stellar 
#    objects 
#    Inputs:
#        baseflux = flux array to base the corrections on
#        initflux = flux array to be corrected
#        basetype = either 'mean' or 'median' which dictate if the mean or 
#                    median value is used for the correction 
#    Output:
#        Flux array with values crudely corrected for differences in seeing
#        (average flux should now be the same for each epoch). '''
#    if avgtype == 'mean':
#        avgfluxperepoch = np.mean(basemag, axis=0)#for UDS
#        avgflux = np.mean(basemag)
#        const = avgflux-avgfluxperepoch
#    elif avgtype == 'median':
#        avgfluxperepoch = np.median(basemag, axis=0)#for UDS
#        avgflux = np.median(basemag)
#        const = avgflux-avgfluxperepoch
#    else:
#        print('Invalid basetype')
#        return
#    return initmag + const[None,:]
#
#def err_correct(flux, fluxerr, fluxnew):
#    ''' Function that applies a correction to the array of error values that 
#    matches the correction applied to the corresponding array of fluxes.
#    Inputs:
#        flux = initial flux array before any corrections were applied
#        fluxerr = initial flux err array
#        fluxcorr = array of fluxes that have been corrected
#    Output:
#        Flux error array with values crudely corrected '''
#        
#    return fluxnew * (fluxerr/flux)
#
#def err_correct_flux(oldflux, fluxerr):
#    ''' Function that applies a correction to the array of error values that 
#    matches the correction applied to the corresponding array of fluxes.
#    Inputs:
#        flux = initial flux array before any corrections were applied
#        fluxerr = initial flux err array
#        fluxcorr = array of fluxes that have been corrected
#    Output:
#        Flux error array with values crudely corrected '''
#    avgflux = np.mean(oldflux, axis=1)
#    return fluxerr/avgflux[:,None]
#
#def mod_z_score(arr):
#    medx = np.median(arr)
#    mad = median_absolute_deviation(arr)
#    zvalues = np.array([(0.6745*(x-medx))/mad for x in arr])
#    return zvalues
#
#def find_outliers(flux, tbdata, bins, threshold=6):
#    ### Bin data ###
#    allmodz = []
#    tbnew = np.recarray([0], dtype=tbdata.dtype, formats=tbdata.formats)
#    for n, binedge in enumerate(bins):
#        if n==np.size(bins)-1:
#            break
#        fluxbin1, tbbin1 = fluxbin(binedge, bins[n+1], flux, tbdata)
#        #calulate mad values in the bins
#        mad1 = median_absolute_deviation(fluxbin1, axis=1) #for UDS
#        modz = mod_z_score(mad1)
#        tbnew = np.hstack([tbnew, tbbin1])
#        allmodz = np.append(allmodz, modz)
#    outliers = allmodz>=threshold
#    return outliers, tbnew, allmodz


def no_ticks():
    ''' Function that removes all ticks from axes, useful on image plots '''
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
#
#def radial_profile(data, center):
#    y, x = np.indices((data.shape)) #create coordinate grid
#    r = np.sqrt((x - center[0])**2 + (y - center[1])**2) #get radius values for grid
#    r = r.astype(np.int)
#
#    tbin = np.bincount(r.ravel(), data.ravel()) # counts number of times value
#                                                # of radius occurs in the psf
#                                                # weighted by the data
#    nr = np.bincount(r.ravel()) # counts number of radii values in psf
#    radialprofile = tbin / nr # as weighted is r*data then get profile by 
#                              # dividing by unweighted counts of r values.
#    return radialprofile 
#
#def plot_median_line(fluxn, tbdata, statistic='MAD',createplot=True):
#    bins = np.array([13, 15])
#    bins = np.append(bins, np.arange(16,24,0.2))
#    bins = np.append(bins, [24])
#    
#    bins = 10**((30-bins)/2.5)
#    bins = np.flip(bins, axis=0)
#    #bins = bins[16:44] #because of flux limit
#    
#    ### Bin data ###
#    allmedstat = np.array([])
#    for n, binedge in enumerate(bins):
#    #    print(binedge)
#        if n==np.size(bins)-1:
#            break
#        mag, bindata = fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
#        if statistic == 'excess':
#            magerr = fluxerr5_stacks(bindata) #make error array
#            nmag, nmagerr = normalise_flux_and_errors(mag, magerr)
#        else:
#            nmag = normalise_flux(mag)
#        
#        if statistic == 'std':
#            binstat = np.std(nmag, axis=1)
#        elif statistic == 'excess':
#            binstat = normsigmasq(nmag, nmagerr)
#        elif statistic == 'MAD':
#            binstat = median_absolute_deviation(nmag, axis=1)
#        elif statistic == 'var':
#            binstat = np.var(nmag, axis=1, ddof=1)
#        else:
#            print('Unrecognised statistic entered')
#            return
#        statmed = np.nanmedian(binstat)
#        allmedstat = np.append(allmedstat, statmed)
#    
#    if createplot==True:
#        plt.plot(bins[0:42], allmedstat, 'k--')
#    return bins, allmedstat
#
#def plot_median_line_stars(fluxn, tbdata, sflux, sdata, statistic='MAD',createplot=True):
#    bins = np.arange(13,25,0.2)
#    bins = np.append(bins, [25])
#    
#    bins = 10**((30-bins)/2.5)
#    bins = np.flip(bins, axis=0)
#    #bins = bins[16:44] #because of flux limit
#    
#    ### Bin data ###
#    allmedstat = np.array([])
#    for n, binedge in enumerate(bins):
#    #    print(binedge)
#        if n==np.size(bins)-1:
#            break
#        
#        # bin both set of data
#        gmag, bindata = fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
#        smag, sbindata = fluxbin(binedge, bins[n+1], sflux, sdata) #bindata
#        
#        # combine into one array
#        mag = np.vstack((gmag, smag))
#        
#        if statistic == 'excess':
#            gmagerr = fluxerr5_stacks(bindata) #make error array
#            smagerr = fluxerr5_stacks(sdata)
#            magerr = np.vstack((gmagerr, smagerr)) #combine
#            
#            nmag, nmagerr = normalise_flux_and_errors(mag, magerr)
#        else:
#            nmag = normalise_flux(mag)
#        
#        if statistic == 'std':
#            binstat = np.std(nmag, axis=1)
#        elif statistic == 'excess':
#            binstat = normsigmasq(nmag, nmagerr)
#        elif statistic == 'MAD':
#            binstat = median_absolute_deviation(nmag, axis=1)
#        elif statistic == 'var':
#            binstat = np.var(nmag, axis=1, ddof=1)
#        else:
#            print('Unrecognised statistic entered')
#            return
#        statmed = np.nanmedian(binstat)
#        allmedstat = np.append(allmedstat, statmed)
#        
#    if createplot==True:
#        plt.plot(bins[0:np.size(bins)-1], allmedstat, 'k--')
#    return bins, allmedstat
#
#def plot_median_line_stars_J(fluxn, tbdata, sflux, sdata, statistic='MAD'):
#    bins = np.arange(13,28,0.2)
#    bins = np.append(bins, [28,29,30])
#    
#    bins = 10**((30-bins)/2.5)
#    bins = np.flip(bins, axis=0)
#    #bins = bins[16:44] #because of flux limit
#    
#    ### Bin data ###
#    allmedstat = np.array([])
#    for n, binedge in enumerate(bins):
#    #    print(binedge)
#        if n==np.size(bins)-1:
#            break
#        
#        # bin both set of data
#        gmag, bindata = fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
#        smag, sbindata = fluxbin(binedge, bins[n+1], sflux, sdata) #bindata
#        
#        # combine into one array
#        mag = np.vstack((gmag, smag))
#        
#        if statistic == 'excess':
#            gmagerr = jfluxerr4_stacks(bindata) #make error array
#            smagerr = jfluxerr4_stacks(sdata)
#            magerr = np.vstack((gmagerr, smagerr)) #combine
#            
#            nmag, nmagerr = normalise_flux_and_errors(mag, magerr)
#        else:
#            nmag = normalise_flux(mag)
#        
#        if statistic == 'std':
#            binstat = np.std(nmag, axis=1)
#        elif statistic == 'excess':
#            binstat = normsigmasq(nmag, nmagerr)
#        elif statistic == 'MAD':
#            binstat = median_absolute_deviation(nmag, axis=1)
#        elif statistic == 'var':
#            binstat = np.var(nmag, axis=1, ddof=1)
#        else:
#            print('Unrecognised statistic entered')
#            return
#        statmed = np.nanmedian(binstat)
#        allmedstat = np.append(allmedstat, statmed)
#        
##    plt.plot(bins[0:np.size(bins)-1], allmedstat, 'k--')
#    return bins, allmedstat
#
#def my_chisquare_char(flux, char_var):
#    ''' Function that calculates the chi^2 of an object using the characteristic
#    variance for its average flux as the error '''
#    fluxn = normalise_flux(flux)
#    meanflux = np.nanmean(fluxn, axis=1)
#    top = np.square(fluxn-meanflux[:,None])
#    chi = np.nansum(top/char_var, axis=1)
#    return chi
#
#def my_chisquare_epoch(flux, sigsq):
#    ''' Function that calculates the chi^2 of an object using the sigma 
#    calculated from the spread of data points for its average flux as the error '''
#    sigsqn = np.tile(sigsq, [len(flux), 1])
#    fluxn, sigsqn = normalise_flux_and_errors(flux, sigsqn)
#    meanflux = np.nanmean(fluxn, axis=1)
#    top = np.square(fluxn-meanflux[:,None])
#    chi = np.nansum(top/(sigsqn**2), axis=1)
#    return chi


#def quadrants(initdata,sem):
#    
#    ira = initdata['X_IMAGE_'+sem]
#    idec = initdata['Y_IMAGE_'+sem]
#
#    ### define bounds of quadrants ###
#    midra = 12450
#    middec = 13310
#    
#    ### create masks for quadrant ###
#    mask1 = ira < midra
#    mask2 = idec < middec
#    quad1data = initdata[mask1*mask2]
#    
#    mask1 = ira >= midra
#    mask2 = idec < middec
#    quad2data = initdata[mask1*mask2]
#    
#    mask1 = ira < midra
#    mask2 = idec >= middec
#    quad3data = initdata[mask1*mask2]
#    
#    mask1 = ira >= midra
#    mask2 = idec >= middec
#    quad4data = initdata[mask1*mask2]
#    
#    return quad1data, quad2data, quad3data, quad4data
#
#def remove_edges(tbdata):
#    
#    ### Set X limits ###
#    x = tbdata['X_IMAGE_05B']
#    xmask1 = x > 1000
#    xmask2 = x < 24000
#    xmask = xmask1 * xmask2.astype(bool)
#    tbdata = tbdata[xmask]
#    
#    ### Set Y limits ###
#    y = tbdata['Y_IMAGE_05B']
#    ymask1 = y > 2000
#    ymask2 = y < 24500
#    ymask = ymask1 * ymask2.astype(bool)
#    tbdata = tbdata[ymask]
#    
#    return tbdata
#
#def flux_split(tbdata, half):
#    ''' Function to give either the sources greater that 1e4 or less than 1e4'''
#    flux = flux4_stacks(tbdata)
#    meanflux = np.nanmean(flux, axis=1)
#    if half == 'upper':
#        tbdata = tbdata[meanflux >= 8e3]
#    elif half == 'lower':
#        tbdata = tbdata[meanflux < 8e3]
#    else:
#        print('Invalid half string, should be upper or lower')
#    return tbdata


#
#def remove_low_flux(flux, tbdata):
#    ''' 
#        Function to remove objects that average below the detection limit
#    '''
#    avgflux = np.nanmean(flux, axis=1)
#    mask = avgflux >= 10**((30-25.3)/2.5)
#    flux = flux[mask]
#    tbdata = tbdata[mask]
#    return flux, tbdata