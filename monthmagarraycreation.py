#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:16:47 2018

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
tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best_err3.fits')[1].data
chandata = fits.open('mag_flux_tables/month_xray_mag_flux_table_best_err3.fits')[1].data
sdata = fits.open('mag_flux_tables/month_stars_mag_flux_table_err3.fits')[1].data

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
fluxerrn = vari_funcs.magerr5_months_errdiff(tbdata)
fluxerrchan = vari_funcs.magerr5_months_errdiff(chandata)
sfluxerr = vari_funcs.magerr5_months_errdiff(sdata)
#%% Save arrays
## Change 99s to nans so they are ignored ###
mask1 = fluxerrn >= 95
mask2 = fluxn == 99
mask = mask1 + mask2
fluxn[mask] = np.nan
fluxerrn[mask] = np.nan
mask1 = fluxerrchan >= 95
mask2 = fluxchann == 99
mask = mask1 + mask2
fluxchann[mask] = np.nan
fluxerrchan[mask] = np.nan
mask1 = sfluxerr >= 95
mask2 = sfluxn == 99
mask = mask1 + mask2
sfluxn[mask] = np.nan
sfluxerr[mask] = np.nan

## Remove rows where all are nans ###
mask = ~np.isnan(np.nanmean(fluxn, axis=0))
fluxn = fluxn[:,mask]
fluxerrn = fluxerrn[:,mask]
mask = ~np.isnan(np.nanmean(fluxchann, axis=0))
fluxchann = fluxchann[:,mask]
fluxerrchan = fluxerrchan[:,mask]
mask = ~np.isnan(np.nanmean(sfluxn, axis=0))
sfluxn = sfluxn[:,mask]
sfluxerr = sfluxerr[:,mask]


np.save('mag_arrays/monthmagarray', fluxn)
np.save('mag_arrays/monthmagerrdiffarray', fluxerrn)
np.save('mag_arrays/monthxraymagarray', fluxchann)
np.save('mag_arrays/monthxraymagerrdiffarray', fluxerrchan)
np.save('mag_arrays/monthstarmagarray', sfluxn)
np.save('mag_arrays/monthstarmagerrdiffarray', sfluxerr)