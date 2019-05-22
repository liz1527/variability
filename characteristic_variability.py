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
#from scipy.stats import chisquare
plt.close('all') #close any open plots

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


vari_funcs.flux_variability_plot(fluxn, fluxchann, 'var', starflux=sfluxn, 
                                 stars=True, normalised=True, scale='symlog')
#bins, medvar = vari_funcs.plot_median_line(fluxn, tbdata, statistic='var')

bins, medvar = vari_funcs.plot_median_line_stars(fluxn, tbdata, sfluxn, sdata, statistic='var')

#
vari_funcs.flux_variability_plot(fluxn, fluxchann, 'chisq',
                                   fluxerr = fluxerr,
                                   starfluxerr = sfluxerr,
                                    starflux=sfluxn, stars=True,
                                    chanerr = fluxerrchan,
                                    normalised=True, scale='symlog')
plt.text(5e2, 5e4, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{i}^{2}}}$')

def my_chisquare(flux, char_var):
    fluxn = vari_funcs.normalise_flux(flux)
    meanflux = np.nanmean(fluxn, axis=1)
    top = np.square(fluxn-meanflux[:,None])
    chi = np.nansum(top/char_var, axis=1)
    return chi

def my_chisquare_charanderr(flux, fluxerr, char_var):
    fluxn, fluxerrn = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    meanflux = np.nanmean(fluxn, axis=1)
    top = np.square(fluxn-meanflux[:,None])
    bot = np.square(fluxerrn) + char_var
    chi = np.nansum(top/bot, axis=1)
    return chi

### Put characteristic variance into chisquare ###
newmedvar = np.empty(np.shape(medvar))
finalmed = np.empty(np.shape(medvar))
for n, binedge in enumerate(bins[0:np.size(bins)-1]):
    flux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
    fluxchan, binchan = vari_funcs.fluxbin(binedge, bins[n+1], fluxchann, chandata) #bindata
    sflux, sbindata = vari_funcs.fluxbin(binedge, bins[n+1], sfluxn, sdata) #bindata
    # get errors
    fluxerr = vari_funcs.fluxerr5_stacks(bindata)
    fluxchanerr = vari_funcs.fluxerr5_stacks(binchan)
    sfluxerr = vari_funcs.fluxerr5_stacks(sbindata)
    
    meanflux = np.nanmean(flux, axis=1)
    meanchan = np.nanmean(fluxchan, axis=1)
    meansflux = np.nanmean(sflux, axis=1)
    chisq = my_chisquare(flux, medvar[n])
    chisqchan = my_chisquare(fluxchan, medvar[n])
    schisq = my_chisquare(sflux, medvar[n])
    
    ### plot ###
    plt.figure(3,figsize=[8,8])
    plt.plot(meanflux, chisq, 'b+',zorder=2)
    plt.plot(meanchan, chisqchan, 'ro', zorder=3, mfc='None', markersize=10)
    plt.plot(meansflux, schisq, 'm*', zorder=1, mfc='None', markersize=10)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Chi Squared')
    plt.xlabel('Mean Flux')
    plt.title('1st iteration')
    plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
    
    ### remove any with chisq > 70 ###
    newgflux = flux[chisq<50,:]
    newsflux = sflux[schisq<50,:]
    newflux = np.vstack((newgflux, newsflux))
    newfluxn = vari_funcs.normalise_flux(newflux)
    vary = np.var(newfluxn, axis=1, ddof=1)
    newmedvar[n] = np.nanmedian(vary)

    newchisq = my_chisquare(flux, newmedvar[n])
    newchisqchan = my_chisquare(fluxchan, newmedvar[n])
    newschisq = my_chisquare(sflux, newmedvar[n])
    
    ### plot new ###
    plt.figure(4, figsize=[8,8])
    plt.plot(meanflux, newchisq, 'b+',zorder=2)
    plt.plot(meanchan, newchisqchan, 'ro', zorder=3, mfc='None', markersize=10)
    plt.plot(meansflux, newschisq, 'm*', zorder=1, mfc='None', markersize=10)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Chi Squared')
    plt.xlabel('Mean Flux')
    plt.title('2nd iteration')
    plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
    
    ### plot new variance ###
    plt.figure(5, figsize=[8,8])
    #get errors
    binfluxn, binfluxerrn = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    binfluxchann, binfluxchanerrn = vari_funcs.normalise_flux_and_errors(fluxchan, fluxchanerr)
    binsfluxn, binsfluxerrn = vari_funcs.normalise_flux_and_errors(sflux, sfluxerr)
    #get variance
    sig = np.var(binfluxn, axis=1, ddof=1)
    sigchan = np.var(binfluxchann, axis=1, ddof=1)
    ssig = np.var(binsfluxn, axis=1, ddof=1)
    sigreal = sig - newmedvar[n]
    sigrealchan = sigchan - newmedvar[n]
    ssigreal = ssig - newmedvar[n]
    #plot
    plt.plot(meanflux, sigreal, 'b+', zorder=2)
    plt.plot(meanchan, sigrealchan, 'ro', zorder=3, mfc='None', markersize=10)
    plt.plot(meansflux, ssigreal, 'm*', zorder=1, mfc='None', markersize=10)
    plt.yscale('symlog', linthreshy=0.0001)
    plt.xscale('log')
    plt.ylabel(r'$\sigma^{2}_{real}$')
    plt.xlabel('Mean Flux')
    plt.text(1e5, 1e0, r'$\sigma^{2}_{real} = \sigma^{2} - \sigma_{noise}^{2}$')
    plt.tight_layout()
    finalmed[n] = np.nanmedian(sigreal)

#plt.plot(bins[0:42], finalmed, 'k--',zorder=3)
















