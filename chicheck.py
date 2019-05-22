#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 10:35:16 2018

code to test chisquare

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_corr.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_corr.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_corr.fits')[1].data

## Create arrays of flux values ###
fluxn = vari_funcs.flux5_stacks(tbdata)
fluxchann = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
fluxn, tbdata = vari_funcs.noneg(fluxn, tbdata)
fluxchann, chandata = vari_funcs.noneg(fluxchann, chandata)
sfluxn, sdata = vari_funcs.noneg(sfluxn, sdata)

### get error arrays and correct them ###
fluxerrn = vari_funcs.fluxerr5_stacks_corr(tbdata)
fluxerrchann = vari_funcs.fluxerr5_stacks_corr(chandata)
sfluxerrn = vari_funcs.fluxerr5_stacks_corr(sdata)

fluxerr = vari_funcs.fluxerr5_stacks(tbdata)
fluxerrchan = vari_funcs.fluxerr5_stacks(chandata)
sfluxerr = vari_funcs.fluxerr5_stacks(sdata)

def single_min_chi_sq(mag, magerr):
    avgmag = np.nanmedian(mag) #use median mags as a start for the expected model
    testchange = np.linspace(-0.1,0.1,101)
    testexpect = testchange + avgmag
    chisq = np.array([(np.square(mag-testexpect[n]))/np.square(magerr) for n in range(101)])
    sumchisq = np.nansum(chisq, axis=1)
    minchisq = np.min(sumchisq)
    expect = testexpect[sumchisq == minchisq]
    # plot the light curve and best fit line to check chisq visually
    plt.figure()
    t = np.arange(1,9)
    plt.errorbar(t, mag, magerr, fmt='x')
    plt.hlines(expect, 1, 8)
    return minchisq 

meansmagnew = np.mean(sfluxn, axis=1)
meanmagnew = np.mean(fluxn, axis=1)
meanchanmagnew = np.mean(fluxchann, axis=1)

sfluxnorm, sfluxerrnorm = vari_funcs.normalise_flux_and_errors(sfluxn, sfluxerrn)
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(fluxn, fluxerrn)
fluxchannorm, fluxerrchannorm = vari_funcs.normalise_flux_and_errors(fluxchann, fluxerrchann)

schisq = single_min_chi_sq(sfluxnorm[53], sfluxerrnorm[53])
