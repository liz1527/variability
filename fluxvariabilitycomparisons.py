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
fluxchann = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
fluxn, tbdata = vari_funcs.semfluxlim(fluxn, tbdata)
fluxchann, chandata = vari_funcs.semfluxlim(fluxchann, chandata)
sfluxn, sdata = vari_funcs.semfluxlim(sfluxn, sdata)
#fluxn, tbdata = vari_funcs.noneg(fluxn, tbdata)
#fluxchann, chandata = vari_funcs.noneg(fluxchann, chandata)
#sfluxn, sdata = vari_funcs.noneg(sfluxn, sdata)

### get error arrays and correct them ###
#fluxerrn = vari_funcs.fluxerr5_stacks_corr(tbdata)
#fluxerrchann = vari_funcs.fluxerr5_stacks_corr(chandata)
#sfluxerrn = vari_funcs.fluxerr5_stacks_corr(sdata)

fluxerrn = vari_funcs.fluxerr5_stacks(tbdata)
fluxerrchann = vari_funcs.fluxerr5_stacks(chandata)
sfluxerrn = vari_funcs.fluxerr5_stacks(sdata)
#
#### create an error array that is just the flux depths ###
## read in depths array
#depths = np.load('fluxdepthsconv_PSF.npy')
## create empty arrays the correct shape
#fluxerrn = np.zeros(np.shape(fluxn)) + depths[None,:]
#fluxerrchann = np.zeros(np.shape(fluxchann)) + depths[None,:]
#sfluxerrn = np.zeros(np.shape(sfluxn)) + depths[None,:]

### Open the fits files and get data ###
tbdata2 = fits.open('mag_flux_tables/mag_flux_table_best.fits')[1].data
chandata2 = fits.open('mag_flux_tables/xray_mag_flux_table_best.fits')[1].data
sdata2 = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data

## Create arrays of flux values ###
flux = vari_funcs.flux5_stacks(tbdata2)
fluxchan = vari_funcs.flux5_stacks(chandata2) # for chandra non-stellar objects
sflux = vari_funcs.flux5_stacks(sdata2)

### remove values that are negative ###
flux, tbdata2 = vari_funcs.semfluxlim(flux, tbdata2)
fluxchan, chandata2 = vari_funcs.semfluxlim(fluxchan, chandata2)
sflux, sdata2 = vari_funcs.semfluxlim(sflux, sdata2)

### get error arrays and correct them ###
#fluxerr = vari_funcs.fluxerr5_stacks_corr(tbdata2)
#fluxerrchan = vari_funcs.fluxerr5_stacks_corr(chandata2)
#sfluxerr = vari_funcs.fluxerr5_stacks_corr(sdata2)

fluxerr = vari_funcs.fluxerr5_stacks(tbdata2)
fluxerrchan = vari_funcs.fluxerr5_stacks(chandata2)
sfluxerr = vari_funcs.fluxerr5_stacks(sdata2)

fig, _ = vari_funcs.flux_variability_plot(flux, fluxchan, 'excess',
                                       fluxerr = fluxerr, flux2=fluxn,
                                       fluxerr2 = fluxerrn, fluxchan2=fluxchann,
                                       chanerr2=fluxerrchann, comparison=True,
                                       starfluxerr = sfluxerr,
                                            starflux=sflux, stars=True, 
                                            chanerr = fluxerrchan,
                                            normalised=True)

bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24])

bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)
bins = bins[16:44] #because of flux limit
### Bin data ###
allmedexcess = np.array([])
allmedexcesscorr = np.array([])
for n, binedge in enumerate(bins):
#    print(binedge)
    if n==np.size(bins)-1:
        break
    mag, bindata = vari_funcs.fluxbin(binedge, bins[n+1], flux, tbdata2) #bindata
    magcorr, bincorr = vari_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
#    magerrcorr = vari_funcs.fluxerr5_stacks_corr(bincorr) #make error array
    magerrcorr = vari_funcs.fluxerr5_stacks(bincorr) #make error array
#    magerr = vari_funcs.fluxerr5_stacks_corr(bindata) #make error array
    magerr = vari_funcs.fluxerr5_stacks(bindata) #make error array
    nmag, nmagerr = vari_funcs.normalise_flux_and_errors(mag, magerr)
    nmagcorr, nmagerrcorr = vari_funcs.normalise_flux_and_errors(magcorr, magerrcorr)
    binexess = vari_funcs.normsigmasq(nmag, nmagerr)
    binexesscorr = vari_funcs.normsigmasq(nmagcorr, nmagerrcorr)
    medexcess = np.nanmedian(binexess)
    allmedexcess = np.append(allmedexcess, medexcess)
    medexcesscorr = np.nanmedian(binexesscorr)
    allmedexcesscorr = np.append(allmedexcesscorr, medexcesscorr)
    
plt.plot(bins[0:26], allmedexcess, 'k--')
plt.plot(bins[0:26], allmedexcesscorr, 'k')

#%% Want to see what chi-square plot looks like with new errors
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
