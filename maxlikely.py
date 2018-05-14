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
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

## Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best_err3.fits')[1].data

# Extract magnitude table and error table
mag = vari_funcs.mag5_months(tbdata)
magerr = vari_funcs.magerr5_months_errdiff(tbdata)

## Change 99s to nans so they are ignored ###
mask = magerr >= 90
mag[mask] = np.nan
magerr[mask] = np.nan
mask = ~np.isnan(np.nanmean(mag, axis=0))
mag = mag[:,mask]
magerr = magerr[:,mask]

magnorm = vari_funcs.normalise_mag(mag)
magerrnorm = vari_funcs.err_correct(mag, magerr, magnorm)

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
#    plt.hlines(Lstd, 0, 3.5)
    plt.vlines(sig+err, np.min(L), np.max(L))
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
##get just one object to try with first
#testmag = magnorm[:,1]
#testmagerr = magerrnorm[:,1]
#
## set up array of test variabilities
#posvar = np.arange(0,3.5,0.01)
#
#### get mean magnitude
#meanmag = np.nanmean(testmag)
#
#sig2, err2 = maximum_likelihood_fig(testmag, testmagerr, meanmag, posvar)
#%% All points
posvar = np.arange(0,3.5,0.01)
#start = time.time()

numobs = np.shape(mag)[1]
meanmag = np.nanmean(magnorm, axis=0)
out = np.array([maximum_likelihood(magnorm[:,n], magerrnorm[:,n], meanmag[n], posvar) for n in range(numobs)])

#%% Plots
out2 = out
out2[out2[:,0]==0] = np.nan
meanmag = np.nanmean(mag, axis=0)

plt.figure()
plt.errorbar(meanmag, out2[:,0], out2[:,1], fmt='x', zorder=0)
plt.yscale('log')
#plt.xlim(xmin=9, xmax=26)
#plt.ylim(ymin=2e-6, ymax=3e0)
plt.ylabel('Maximum Likelihood')
plt.xlabel('Mean Magnitude')
#plt.figure()
plt.scatter(meanmag, out2[:,0],c='k', marker='x', zorder=2)
plt.yscale('log')
plt.xlim(xmin=12, xmax=26)
plt.ylim(ymin=4e-3, ymax=4e0)
plt.ylabel('Maximum Likelihood')
plt.xlabel('Mean Magnitude')

end = time.time()
print(end-start)