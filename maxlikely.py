#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields
plt.style.use('dark_background')

## Open the fits files and get data ###
#tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

# Extract magnitude table and error table
flux = vari_funcs.flux5_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)


fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#%% define functions

def maximum_likelihood_fig(testmag, testmagerr, meanmag, posvar):

    # Calculate likelihood curve
    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
#    Lstd = np.std(L)
#    idx = (np.abs(L+Lstd).argmin())
#    sigstd = posvar[idx]
#    err = np.abs(sig+sigstd)
    err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
    plt.figure()
    plt.plot(posvar, L)
    plt.vlines(sig, np.min(L), np.max(L), color='w')
#    plt.vlines(sig+err, np.min(L), np.max(L))
#    plt.vlines(sig-err, np.min(L), np.max(L))
    plt.ylim(ymin=0)
    plt.xlabel(r'$\sigma$')
    plt.ylabel('L')
    return sig, err

def maximum_likelihood(testmag, testmagerr, meanmag, posvar):

    # Calculate likelihood curve
    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
    if np.sum(L) == 0:
        return sig, np.nan
    else:
        err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
        return sig, err

#%% Single curve
#get just one object to try with first
testmag = fluxnorm[300,:]
testmagerr = fluxerrnorm[300,:]

# set up array of test variabilities
posvar = np.arange(0,2,0.01)

### get mean magnitude
meanmag = np.nanmean(testmag)

sig2, err2 = maximum_likelihood_fig(testmag, testmagerr, meanmag, posvar)
##%% All points
#posvar = np.arange(0,2,0.01)
##start = time.time()
#
#numobs = np.shape(fluxnorm)[0]
#meanflux = np.nanmean(fluxnorm, axis=1)
#out = np.array([maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], meanflux[n], posvar) for n in range(numobs)])
#
##%% Plots
#out2 = out
##out2[out2[:,0]==0] = np.nan
#meanflux = np.nanmean(flux, axis=1)
#
#plt.figure()
#plt.errorbar(meanflux, out2[:,0], out2[:,1], fmt='x', zorder=0)
#plt.ylabel('Maximum Likelihood')
#plt.xlabel('Mean Magnitude')
#plt.xscale('log')
##plt.figure()
#plt.scatter(meanflux, out2[:,0],c='k', marker='x', zorder=2)
#plt.ylabel('Maximum Likelihood')
#plt.xlabel('Mean Magnitude')

end = time.time()
print(end-start)