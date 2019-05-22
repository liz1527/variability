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
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_1519match.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_1519match.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_1519match.fits')[1].data

## Create arrays of flux values ###
fluxn = vari_funcs.flux5_stacks(tbdata)
fluxchann = vari_funcs.flux5_stacks(chandata) 
sfluxn = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
fluxn, tbdata = vari_funcs.noneg(fluxn, tbdata)
fluxchann, chandata = vari_funcs.noneg(fluxchann, chandata)
sfluxn, sdata = vari_funcs.noneg(sfluxn, sdata)

fluxerr = vari_funcs.fluxerr5_stacks(tbdata)
fluxerrchan = vari_funcs.fluxerr5_stacks(chandata)
sfluxerr = vari_funcs.fluxerr5_stacks(sdata)


fig = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'mad',
                                            starflux=sfluxn, stars=True,
                                            normalised=True, scale='log')

#fig = plt.figure(figsize=[9,6])
#avgflux = np.nanmean(fluxn, axis=1)
#fluxnorm = vari_funcs.normalise_flux(fluxn)
#sigma = np.std(fluxnorm, axis=1)
#plt.scatter(avgflux, sigma, marker='+')
#plt.xscale('log')
##plt.yscale('log')
#plt.xlabel('Mean Flux')
#plt.ylabel('sigma')
fig, _ = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'excess',
                                       fluxerr = fluxerr,
                                       starfluxerr = sfluxerr,
                                        starflux=sfluxn, stars=True, 
                                        chanerr = fluxerrchan,
                                        normalised=True, scale='symlog')
##plt.ylim(ymin=-0.05, ymax=0.05)
fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux)

bins, medexcess = vari_funcs.plot_median_line(fluxn, tbdata, statistic='excess')


##%% Want to see what chi-square plot looks like with new errors
#def min_chi_sq(mag, magerr):
#    avgmag = np.nanmedian(mag, axis=1) #use median mags as a start for the expected model
#    # create array that checks model around the median value
#    testexpect = np.tile(np.array([avgmag]).transpose(), [1,50])
##    testrange = 0.1*avgmag
#    for n in range(len(avgmag)):
#        values = np.linspace(-0.1, 0.1) #np.linspace(-testrange[n], testrange[n]) 
#        if n==0:
#            testchanges = values
#        else:
#            testchanges = np.vstack([testchanges, values])
#    testexpect += testchanges
#    # Find chi-squared values with all possible models
#    chisq = np.array([(np.square(mag-testexpect[:,n,None]))/np.square(magerr) for n in range(50)])
#    chisq = np.nansum(chisq, axis=2) #sum to complete chi squared calculation
#    chisqmin = np.nanmin(chisq, axis=0) #find the minimum chi-squared value
#    return chisqmin
#
#def chi_sq(mag, magerr):
#    avgmag = np.nanmedian(mag, axis=1) #use median mags as a start for the expected model
#    chisq = np.square(mag-avgmag[:,None])/np.square(magerr)
#    sumchisq = np.nansum(chisq, axis=1)
#    return sumchisq
#
#meansmagnew = np.mean(sfluxn, axis=1)
#meanmagnew = np.mean(fluxn, axis=1)
#meanchanmagnew = np.mean(fluxchann, axis=1)
#
#sfluxnorm, sfluxerrnorm = vari_funcs.normalise_flux_and_errors(sfluxn, sfluxerrn)
#fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(fluxn, fluxerrn)
#fluxchannorm, fluxerrchannorm = vari_funcs.normalise_flux_and_errors(fluxchann, fluxerrchann)
#
#schisq = min_chi_sq(sfluxnorm, sfluxerrnorm)
#chisq = min_chi_sq(fluxnorm, fluxerrnorm)
#chanchisq = min_chi_sq(fluxchannorm, fluxerrchannorm)
#
##schisq = chi_sq(sfluxnorm, sfluxerrnorm)
##chisq = chi_sq(fluxnorm, fluxerrnorm)
##chanchisq = chi_sq(fluxchannorm, fluxerrchannorm)
#
#plt.figure()
#plt.plot(meansmagnew, schisq, 'm*', markersize=10, mfc='None')
#plt.plot(meanmagnew, chisq, 'b+')
#plt.plot(meanchanmagnew, chanchisq, 'ro', markersize=10, mfc='None')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Mean Flux')
#plt.ylabel('Chi squared value')
#
