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

import xray_funcs #for calculating luminosities and alpha_ox

import cross_correlation #for running cross-correlation analysis

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

def no_ticks():
    ''' Function that removes all ticks from axes, useful on image plots '''
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    
def mag2flux(mag):
    ''' Function to convert DR11 magnitudes into fluxes
    Inputs:
        mag = AB magnitude from the DR11 catalogue
    Outputs:
        flux = flux of related to that magnitude 
    '''
    flux = 10**((mag-30)/(-2.5))
    return flux

def split_selection_band(varydata):
    ### Split table into J and K selected ###
    jdata = varydata[varydata['Chi_J'] > 32.08]
    kdata = varydata[varydata['Chi_K'] > 30]
    
    ### Find overlap ###
    jIDs = jdata['ID']
    kIDs = kdata['ID']
    jmask = np.isin(jIDs, kIDs)
    kmask = np.isin(kIDs, jIDs)
    
    ### Mask tables ###
    bothdata = jdata[jmask]
    jdata = jdata[~jmask]
    kdata = kdata[~kmask]
    
    return jdata, kdata, bothdata