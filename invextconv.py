#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:40:02 2018

Code to look at results from convolution

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs_no06 #my module to help run code neatly
from matplotlib.colors import LogNorm
plt.close('all') #close any open plots

combined = fits.open('stars_mag_flux_table.fits')
alldata = combined[1].data
combined = fits.open('stars_mag_flux_convS.fits')
alldataconv = combined[1].data
stars = fits.open('starsfwhm.fits')
sdata = stars[1].data

# remove saturated stars
sdata = sdata[alldata['MAG_APER_5_05B'] >= 12]
alldata = alldata[alldata['MAG_APER_5_05B'] >= 12]
alldataconv = alldataconv[alldataconv['MAG_APER_5_05B'] >= 12]

## Create flux stack
#allflux = vari_funcs_no06.flux5_stacks(alldata)
#allfluxconv = vari_funcs_no06.flux5_stacks(alldataconv)


# Create mag stack
#allflux = vari_funcs_no06.mag5_stacks(alldata)
#allflux, alldata2 = vari_funcs_no06.no99(allflux, alldata)
#allfluxconv = vari_funcs_no06.mag5_stacks(alldataconv)
#allfluxconv, alldataconv2 = vari_funcs_no06.no99(allfluxconv, alldataconv)

# Remove negative values
#allflux[allflux <= 0] = np.nan
#mask = ~np.isnan(allflux).any(axis=1)
#allflux = allflux[mask]
##allfluxconv[allfluxconv <= 0] = np.nan
##mask = ~np.isnan(allfluxconv).any(axis=1)
#allfluxconv = allfluxconv[mask]
#
## Normalise
#allflux = vari_funcs_no06.normalise_flux(allflux)
##allflux = vari_funcs_no06.normalise_mag(allflux)
##allfluxconv = vari_funcs_no06.psf_correct_flux(allflux, allflux, 'median')
#allfluxconv = vari_funcs_no06.normalise_flux(allfluxconv)
##allfluxconv = vari_funcs_no06.normalise_mag(allfluxconv)

avgflux = np.array([np.median(sdata['FWHM_05B']), 
                    np.median(sdata['FWHM_07B']), 
                    np.median(sdata['FWHM_08B']), 
                    np.median(sdata['FWHM_09B']), 
                    np.median(sdata['FWHM_10B']), 
                    np.median(sdata['FWHM_11B']), 
                    np.median(sdata['FWHM_12B'])]) *3600

avgfluxconv = np.array([np.median(alldataconv['FWHM_WORLD_05B']), 
                    np.median(alldataconv['FWHM_WORLD_07B']), 
                    np.median(alldataconv['FWHM_WORLD_08B']), 
                    np.median(alldataconv['FWHM_WORLD_09B']), 
                    np.median(alldataconv['FWHM_WORLD_10B']), 
                    np.median(alldataconv['FWHM_WORLD_11B']), 
                    np.median(alldataconv['FWHM_WORLD_12B'])]) *3600

## find and plot average
#avgflux = np.median(allflux, axis=0)
vari_funcs_no06.avg_lightcurve(avgflux)
#plt.title('Normalised Flux of Stars with Unconvolved')
#plt.ylim(0.9986, 1.0004)
#plt.ylim(21.36, 21.51)
plt.title('Median FWHM of stars before convolution')
plt.ylim(0.74, 0.88)
plt.ylabel('FWHM (arcsec)')
#plt.ylabel('Normalised Flux')
plt.savefig('plots/Lightcurves/FWHMbefore')

#avgfluxconv = np.median(allfluxconv, axis=0)
vari_funcs_no06.avg_lightcurve(avgfluxconv)
#plt.title('Normalised Flux of Stars curve with Convolved')
#plt.ylim(0.9986, 1.0004)
#plt.ylim(21.36, 21.51)
#plt.ylim(0.000205,0.000243)
plt.title('Median FWHM of stars after convolution')
plt.ylim(0.74, 0.88)
plt.ylabel('FWHM (arcsec)')
#plt.ylabel('Normalised Flux')
plt.savefig('plots/Lightcurves/FWHMafter')
