#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:04:52 2017

Variability using magnitudes instead of fluxes

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
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

## Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best.fits')[1].data
chandata = fits.open('mag_flux_tables/month_xray_mag_flux_table_best.fits')[1].data
sdata = fits.open('mag_flux_tables/month_stars_mag_flux_table.fits')[1].data

## Restrict objects to those in the Chandra field ###
#tbdata = vari_funcs.chandra_only(tbdata)
#chandata = vari_funcs.chandra_only(chandata)
#sdata = vari_funcs.chandra_only(sdata)

### Split chandra and xmm ###
#mask = np.isnan(chandata['RA'])
#xmmdata = chandata[mask]
#chandata = chandata[~mask]

### Create arrays of flux values ###
#fluxn = vari_funcs.mag5_stacks(tbdata)
#fluxchann = vari_funcs.mag5_stacks(chandata) # for chandra non-stellar objects
#sfluxn = vari_funcs.mag5_stacks(sdata)
#fluxxmm = vari_funcs.mag5_stacks(xmmdata)

# Create arrays of flux values ###
fluxn = vari_funcs.mag5_months(tbdata)
fluxchann = vari_funcs.mag5_months(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.mag5_months(sdata)

### get error arrays and correct them ###
fluxerrn = vari_funcs.magerr5_months(tbdata)
fluxerrchan = vari_funcs.magerr5_months(chandata)
sfluxerr = vari_funcs.magerr5_months(sdata)

## Change 99s to nans so they are ignored ###
mask = fluxerrn >= 99
fluxn[mask] = np.nan
fluxerrn[mask] = np.nan
mask = fluxerrchan >= 99
fluxchann[mask] = np.nan
fluxerrchan[mask] = np.nan
mask = sfluxerr >= 99
sfluxn[mask] = np.nan
sfluxerr[mask] = np.nan

## Remove rows where all are nans ###
mask = ~np.isnan(np.nanmean(fluxn, axis=0))
fluxn = fluxn[:,mask]
fluxerrn = fluxerrn[:,mask]
mask = ~np.isnan(np.nanmean(fluxchann, axis=0))
fluxchann = fluxchann[:,mask]
mask = ~np.isnan(np.nanmean(sfluxn, axis=0))
sfluxn = sfluxn[:,mask]
sfluxerr = sfluxerr[:,mask]

#### normalise to 1 ###
#fluxn = vari_funcs.normalise(fluxn)
#fluxchann = vari_funcs.normalise(fluxchann)
#sfluxn = vari_funcs.normalise(sfluxn)

#
### Multiply all flux values in a yearstack by the correct constant ###
#fluxcorrn = vari_funcs.psf_correct_mag(fluxn, fluxn, 'median') 
#fluxchancorrn = vari_funcs.psf_correct(fluxn, fluxchann, 'median') 
#sfluxcorrn = vari_funcs.psf_correct(fluxn, sfluxn, 'median') 



### Create Plot for non corrected ###
fig1 = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'mad',
                                            starflux=sfluxn, stars=True)
fig2 = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'excess',
                                       fluxerr = fluxerrn, 
                                       starfluxerr = sfluxerr,
                                            starflux=sfluxn, stars=True,
                                            chanerr = fluxerrchan)
fig2.canvas.mpl_connect('pick_event', vari_funcs.onpickmonth)
fig1.canvas.mpl_connect('pick_event', vari_funcs.onpickmonth)
#
