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
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from matplotlib.colors import LogNorm
plt.close('all') #close any open plots

def singlesigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        sig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux)
    N = np.size(flux)
    sig = [((flux - avgflux)**2 - (baseerr)**2)]# 
    sigsum = np.nansum(sig)
    normsig = sigsum/(N)
    return normsig

alldata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
alldataconv = fits.open('mag_flux_tables/stars_mag_flux_table_1519match.fits')[1].data
colname = 'FWHM_WORLD_'
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
goodstarIDs = np.load('PSF_IDs_original.npy')

#Extract the flux radii and remove negative values
avgfwhm = np.zeros(8)
avgfwhmconv = np.zeros(8)

#mask = np.isin(alldata['DR11_IDs'], goodstarIDs)
#alldata = alldata[mask]
#mask = np.isin(alldataconv['DR11_IDs'], goodstarIDs)
#alldataconv = alldataconv[mask]

for n, sem in enumerate(semesters):
    colnames = colname+sem
    
    # remove stars in old table
    mag = alldata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 19 #removes very faint stars
    mask = mask1 * mask2
    alldata = alldata[mask]
    avgfwhm[n] = np.median(alldata[colnames]) * 3600
    
    # remove stars in new table
    mag = alldataconv['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 19 #removes very faint stars
    mask = mask1 * mask2
    alldataconv = alldataconv[mask]
    avgfwhmconv[n] = np.median(alldataconv[colnames]) * 3600
    
# Create flux stack
allflux = vari_funcs.flux5_stacks(alldata)
allfluxconv = vari_funcs.flux5_stacks(alldataconv)

#allflux, alldata = vari_funcs.semfluxlim(allflux, alldata)
#allfluxconv, alldataconv = vari_funcs.semfluxlim(allfluxconv, alldataconv)

allflux, alldata = vari_funcs.noneg(allflux, alldata)
allfluxconv, alldataconv = vari_funcs.noneg(allfluxconv, alldataconv)

allfluxerr = vari_funcs.fluxerr1_stacks(alldata)
allfluxconverr = vari_funcs.fluxerr1_stacks(alldataconv)

#depths = np.load('fluxdepths.npy')
#allfluxerr = np.zeros(np.shape(allflux)) + depths[None,:]
#depthsconv = np.load('fluxdepthsconv_PSF.npy')
#allfluxconverr = np.zeros(np.shape(allfluxconv)) + depthsconv[None,:]


# Normalise
allflux, allfluxerr = vari_funcs.normalise_flux_and_errors(allflux, allfluxerr)
allfluxconv, allfluxconverr = vari_funcs.normalise_flux_and_errors(allfluxconv, allfluxconverr)

## Find FWHM values
#avgfwhm = np.array([np.median(alldata['FWHM_WORLD_05B']), 
#                    np.median(alldata['FWHM_WORLD_06B']), 
#                    np.median(alldata['FWHM_WORLD_07B']), 
#                    np.median(alldata['FWHM_WORLD_08B']), 
#                    np.median(alldata['FWHM_WORLD_09B']), 
#                    np.median(alldata['FWHM_WORLD_10B']), 
#                    np.median(alldata['FWHM_WORLD_11B']), 
#                    np.median(alldata['FWHM_WORLD_12B'])]) *3600
#
#avgfwhmconv = np.array([np.median(alldataconv['FWHM_WORLD_05B']), 
#                    np.median(alldataconv['FWHM_WORLD_06B']), 
#                    np.median(alldataconv['FWHM_WORLD_07B']), 
#                    np.median(alldataconv['FWHM_WORLD_08B']), 
#                    np.median(alldataconv['FWHM_WORLD_09B']), 
#                    np.median(alldataconv['FWHM_WORLD_10B']), 
#                    np.median(alldataconv['FWHM_WORLD_11B']), 
#                    np.median(alldataconv['FWHM_WORLD_12B'])]) *3600

### find and plot averages ###
#plot FWHM curve before
vari_funcs.avg_lightcurve(avgfwhm, shape='s', size=9)
#plt.title('Median FWHM of stars before convolution')
#plt.ylim(0.73, 0.9)
plt.ylabel('FWHM (arcsec)')
#plt.savefig('plots/Lightcurves/FWHMbefore')

#plot FWHM curve after
vari_funcs.avg_lightcurve(avgfwhmconv)
#plt.title('Median FWHM')
#plt.ylim(0.74, 1.5)
plt.ylabel('FWHM (arcsec)')
plt.xlabel('Semester')
plt.title('')
#plt.savefig('plots/Lightcurves/FWHMafter')

#fwhmmadbefore = median_absolute_deviation(avgfwhm)
#plt.text(0.7, 0.75, 'MAD before = '+str(round(fwhmmadbefore,4)))
#fwhmmadafter = median_absolute_deviation(avgfwhmconv)
#plt.text(0.7, 0.74, 'MAD after = '+str(round(fwhmmadafter,4)))

plt.figure()
#plot flux curve before
#plt.figure()
avgflux = np.nanmedian(allflux, axis=0)
avgerr = np.nanmedian(allfluxerr, axis=0)
vari_funcs.avg_lightcurve(avgflux, avgerr)
plt.title('Normalised Flux of Stars')
#plt.ylim(0.95, 1.06)
plt.ylabel('Normalised Flux')

#plot flux curve after
avgfluxconv = np.median(allfluxconv, axis=0)
avgerrconv = np.median(allfluxconverr, axis=0)
vari_funcs.avg_lightcurve(avgfluxconv, avgerrconv)
plt.title('Normalised Flux of Stars')
plt.ylabel('Normalised Flux')
#plt.ylim(0.95, 1.06)


#flaxvarbefore = singlesigmasq(avgflux, avgerr)
#plt.text(5, 0.989, 'Var before = %.2E' % flaxvarbefore)
#fluxvarafter = singlesigmasq(avgfluxconv, avgerrconv)
#plt.text(5, 0.986, 'Var after = %.2E' % fluxvarafter)