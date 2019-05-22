#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:58:14 2019

Code to check the Chi-squared of full sample after most deviant point is removed

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra_clean_no06.fits')[1].data
#sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

### Remove edges ###
tbdata = vari_funcs.remove_edges(tbdata)
chandata = vari_funcs.remove_edges(chandata)
sdata = vari_funcs.remove_edges(sdata)

#### Limit to Chandra region ###
#tbdata = vari_funcs.chandra_only(tbdata)
#chandata = vari_funcs.chandra_only(chandata)
#sdata = vari_funcs.chandra_only(sdata)

## Create arrays of flux values ###
#flux = vari_funcs.flux5_stacks(tbdata)
#fluxchan = vari_funcs.flux5_stacks(chandata) 
#sflux = vari_funcs.flux5_stacks(sdata)
flux = vari_funcs.flux4_stacks(tbdata)
fluxchan = vari_funcs.flux4_stacks(chandata) 
sflux = vari_funcs.flux4_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Get error arrays ###
#flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
#fluxchan, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata)
#sflux, serr, sdata = vari_funcs.create_quad_error_array(sigtb, sdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
fluxchan, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata, aper=4)
sflux, serr, sdata = vari_funcs.create_quad_error_array(sigtb, sdata, aper=4)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Find difference from 1 ###
fluxn, fluxnerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
diff = abs(fluxn - 1)
fluxchann, channerr = vari_funcs.normalise_flux_and_errors(fluxchan, chanerr)
chandiff = abs(fluxchann - 1)
sfluxn, snerr = vari_funcs.normalise_flux_and_errors(sflux, serr)
sdiff = abs(sfluxn - 1)

### mask max difference ###
mask = np.argmax(diff, axis=1)
flux[np.arange(len(flux)),mask] = np.nan
fluxerr[np.arange(len(flux)),mask] = np.nan
chanmask = np.argmax(chandiff, axis=1)
fluxchan[np.arange(len(fluxchan)),chanmask] = np.nan
chanerr[np.arange(len(fluxchan)),chanmask] = np.nan
smask = np.argmax(sdiff, axis=1)
sflux[np.arange(len(sflux)),smask] = np.nan
serr[np.arange(len(sflux)),smask] = np.nan


### Check chisq plot looks correct ###
fig,chisq = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux)

plt.hlines(27.83, 8e1,1e7, label='Chi=27.83')
plt.hlines(20, 8e1,1e7, label='Chi=20')
plt.legend()

varydata = tbdata[chisq>27.83]
varydata = tbdata[chisq>20]