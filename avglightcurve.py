#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 11:28:23 2017

Code to examine average light curves of objects in certain flux bins.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('mag_flux_table_best.fits')
tbdata = combined[1].data
stars = fits.open('starsnotsatfwhm.fits')
sdata = stars[1].data
stars = fits.open('stars_mag_flux_table.fits')
sdata2 = stars[1].data

## Restrict objects to those in the Chandra field ###
tbdata = vari_funcs.chandra_only(tbdata)

## Create arrays of flux values ###
flux = vari_funcs.k_mag_flux.flux5_stacks(tbdata)
sfwhm = np.stack(([sdata['FWHM_05B'], sdata['FWHM_06B'],
                sdata['FWHM_07B'], sdata['FWHM_08B'],
                sdata['FWHM_09B'], sdata['FWHM_10B'], 
                sdata['FWHM_11B'], sdata['FWHM_12B']]), axis=1)
sflux = vari_funcs.flux5_stacks(sdata2)

### Calculate the average flux for each year within each flux bin ###
avgflux = np.mean(flux, axis=0)
savgfwhm = np.mean(sfwhm, axis=0)
savgflux = np.mean(sflux, axis=0)

### Plot light curves ###
vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux)
plt.title('Mean Lightcurve for All Objects')
vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux)
plt.title('Mean Lightcurve for Stars')
ax = vari_funcs.lightcurve_funcs.avg_lightcurve(savgfwhm)
plt.title('Mean FWHM')
plt.ylabel('FWHM')
ax.invert_yaxis()


### trying to investigate 06B ###
nonzero06fwhm = sdata['FWHM_06B'][np.nonzero(sdata['FWHM_06B'])]
avgfwhm06 = np.mean(nonzero06fwhm)
plt.figure()
plt.hist(nonzero06fwhm)

combined.close()
stars.close()
################################################################

#calculate average flux for each object
#avgflux = np.mean(flux, axis=1) #for UDS

### Define bins for histograms and average light curves ###
#fluxbin1 = vari_funcs.flux_funcs.fluxbin(avgflux, 0, 1e1, flux)
#fluxbin2 = vari_funcs.flux_funcs.fluxbin(avgflux, 1e1, 1e2, flux)
#fluxbin3 = vari_funcs.flux_funcs.fluxbin(avgflux, 1e2, 1e3, flux)
#fluxbin4 = vari_funcs.flux_funcs.fluxbin(avgflux, 1e3, 1e4, flux)
#fluxbin5 = vari_funcs.flux_funcs.fluxbin(avgflux, 1e4, 1e5, flux)
#fluxbin6 = vari_funcs.flux_funcs.fluxbin(avgflux, 1e6, 1e7, flux)

#avgflux1 = np.mean(fluxbin1, axis=0)
#avgflux2 = np.mean(fluxbin2, axis=0)
#avgflux3 = np.mean(fluxbin3, axis=0)
#avgflux4 = np.mean(fluxbin4, axis=0)
#avgflux5 = np.mean(fluxbin5, axis=0)
#avgflux6 = np.mean(fluxbin6, axis=0)

#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux1)
#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux2)
#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux3)
#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux4)
#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux5)
#vari_funcs.lightcurve_funcs.avg_lightcurve(avgflux6)
#
#### Create arrays of flux values ###
#sflux = vari_funcs.flux5_stacks(sdata)
#
##calculate average flux for each object
#savgflux = np.mean(sflux, axis=1) #for UDS
#
#### Define bins for histograms and average light curves ###
#sfluxbin1 = vari_funcs.flux_funcs.fluxbin(savgflux, 0, 1e1, sflux)
#sfluxbin2 = vari_funcs.flux_funcs.fluxbin(savgflux, 1e1, 1e2, sflux)
#sfluxbin3 = vari_funcs.flux_funcs.fluxbin(savgflux, 1e2, 1e3, sflux)
#sfluxbin4 = vari_funcs.flux_funcs.fluxbin(savgflux, 1e3, 1e4, sflux)
#sfluxbin5 = vari_funcs.flux_funcs.fluxbin(savgflux, 1e4, 1e5, sflux)
#sfluxbin6 = vari_funcs.flux_funcs.fluxbin(savgflux, 1e6, 1e7, sflux)
#
#### Calculate the average flux for each year within each flux bin ###
#savgflux = np.mean(sflux, axis=0)
#savgflux1 = np.mean(sfluxbin1, axis=0)
#savgflux2 = np.mean(sfluxbin2, axis=0)
#savgflux3 = np.mean(sfluxbin3, axis=0)
#savgflux4 = np.mean(sfluxbin4, axis=0)
#savgflux5 = np.mean(sfluxbin5, axis=0)
#savgflux6 = np.mean(sfluxbin6, axis=0)
#
#### Plot light curves ###
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux1)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux2)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux3)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux4)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux5)
#vari_funcs.lightcurve_funcs.avg_lightcurve(savgflux6)