#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:19:08 2018

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
fluxerrchan = vari_funcs.fluxerr5_stacks_corr(chandata)
sfluxerr = vari_funcs.fluxerr5_stacks_corr(sdata)

#fig = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'mad',
#                                            starflux=sfluxn, stars=True)

fig = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'excess',
                                       fluxerr = fluxerrn, starfluxerr = sfluxerr,
                                            starflux=sfluxn, stars=True, 
                                            chanerr = fluxerrchan,
                                            normalised=True)

#fig.canvas.mpl_connect('pick_event', vari_funcs.onpick)
#
bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24])

bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

#mag, bindata = vari_funcs.fluxbin(bins[25], bins[26], allmag, tbdata)
#magerr = vari_funcs.magerr5_stacks(bindata)
#errchange, newmagerr = add_chi_err(mag, magerr)

### Bin data ###
allerrchange = np.array([])
allmedexcess = np.array([])
for n, binedge in enumerate(bins):
#    print(binedge)
    if n==np.size(bins)-1:
        break
    mag, bindata = vari_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
    magerr = vari_funcs.fluxerr5_stacks_corr(bindata) #make error array
    nmag, nmagerr = vari_funcs.normalise_flux_and_errors(mag, magerr)
    excess = vari_funcs.normsigmasq(nmag, nmagerr)
    medexcess = np.nanmedian(excess)
    allmedexcess = np.append(allmedexcess, medexcess)
    
plt.plot(bins[0:42], allmedexcess, 'k')

#%% Want to see what chi-square plot looks like with new errors
def min_chi_sq(mag, magerr):
    avgmag = np.nanmedian(mag, axis=1) #use median mags as a start for the expected model
    # create array that checks model around the median value
    testexpect = np.tile(np.array([avgmag]).transpose(), [1,50])
    testchanges = np.linspace(-0.1, 0.1) 
    testexpect += testchanges
    # Find chi-squared values with all possible models
    chisq = np.array([(np.square(mag-testexpect[:,n,None]))/np.square(magerr) for n in range(50)])
    chisq = np.nansum(chisq, axis=2) #sum to complete chi squared calculation
    chisqmin = np.nanmin(chisq, axis=0) #find the minimum chi-squared value
    return chisqmin

meansmagnew = np.mean(sfluxn, axis=1)
meanmagnew = np.mean(fluxn, axis=1)
meanchanmagnew = np.mean(fluxchann, axis=1)
schisq = min_chi_sq(sfluxn, sfluxerr)
chisq = min_chi_sq(fluxn, fluxerrn)
chanchisq = min_chi_sq(fluxchann, fluxerrchan)

plt.figure()
plt.plot(meansmagnew, schisq, 'm*', markersize=10, mfc='None')
plt.plot(meanmagnew, chisq, 'b+')
plt.plot(meanchanmagnew, chanchisq, 'ro', markersize=10, mfc='None')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Mean Flux')
plt.ylabel('Chi squared value')