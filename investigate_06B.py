#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:35:11 2018

Code to isolate variable sources using the bootstraped error bars and chi^2

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_cleaned.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_cleaned.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_cleaned.fits')[1].data
sigtb = Table.read('quad_epoch_sigma_table_cleaned.fits')

## Create arrays of flux values ###
flux = vari_funcs.flux5_stacks(tbdata)
fluxchan = vari_funcs.flux5_stacks(chandata) 
sflux = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
fluxchan, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata)
sflux, serr, sdata = vari_funcs.create_quad_error_array(sigtb, sdata)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Check chisq plot looks correct ###
fig, chisq = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux)

### find those with lightcurves that are > 3*MAD ###
#flux, fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
mad = median_absolute_deviation(flux, axis=1)
med = np.mean(flux, axis=1)
flux06 = flux[:,1]
diff06B = abs(flux06-med)

goodflux = flux[diff06B < 10*mad,:]
badflux = flux[diff06B > 10*mad,:]
bad06 = tbdata[diff06B > 10*mad]
badchi = chisq[diff06B > 10*mad]

plt.plot(np.mean(badflux, axis=1), badchi,'yo',mfc='None')

med = np.nanmedian(goodflux, axis=0)
plt.figure()
plt.plot(med)

med = np.nanmedian(badflux, axis=0)
plt.figure()
plt.plot(med)

#t = Table(bad06)
#t.write('06mad10table.fits')