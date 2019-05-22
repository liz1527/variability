#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 13:51:08 2018

Code to test updates to sigmasq function

@author: ppxee
"""
### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
#from astropy.stats import median_absolute_deviation
#
#def sigmasq(flux, baseerr):
#    ''' Function that calculates the excess varience value for each row in an 
#    array 
#    Inputs:
#        flux = array of fluxes from objects in a number of epochs 
#        baseerr = array of errors that the mean error should be calculated from
#    Output:
#        sig = array of excess variance values for every object '''
#    avgflux = np.nanmean(flux, axis=0)
#    meanerr = np.nanmedian(baseerr, axis=0)
#    meanerrsq = meanerr**2
#    baseerrsq = baseerr**2
#    N = np.size(flux, axis=0)
#    numobs = np.size(flux, axis=1)
##    sig = [((flux[:, None, n]- avgflux[n])**2 - meanerrsq[n])/(N-1) for n in range(numobs)]# 
#    sig = [((flux[:, None, n]- avgflux[n])**2 - baseerrsq[:, None, n])/(N-1) for n in range(numobs)]# 
#    sig = np.array(sig).reshape(numobs, N)
#    sigsum = np.nansum(sig, axis=1)
#    return sigsum
def normsigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        sig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux, axis=0)
    baseerrsq = np.square(baseerr)
    N = np.size(flux, axis=0)
    numobs = np.size(flux, axis=1)
    sig = [((flux[:, n]- avgflux[n])**2 - baseerrsq[:, n]) for n in range(numobs)]# 
    sigsum = np.nansum(sig, axis=1)
    normsig = sigsum/(N*avgflux**2)
    return normsig
## Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best.fits')[1].data
chandata = fits.open('mag_flux_tables/month_xray_mag_flux_table_best.fits')[1].data
sdata = fits.open('mag_flux_tables/month_stars_mag_flux_table.fits')[1].data

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
fluxerrchan = fluxerrchan[:,mask]
mask = ~np.isnan(np.nanmean(sfluxn, axis=0))
sfluxn = sfluxn[:,mask]
sfluxerr = sfluxerr[:,mask]

var = normsigmasq(fluxn, fluxerrn)
plt.figure()
plt.hist(var)