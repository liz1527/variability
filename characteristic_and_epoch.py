#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:19:08 2018

OLD CODE - NOT UPDATED WITH RESTRUCTURING 
Part of the build up to uncertainty calculations

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots

def my_chisquare_char(flux, char_var):
    fluxn = vari_funcs.normalise_flux(flux)
    meanflux = np.nanmean(fluxn, axis=1)
    top = np.square(fluxn-meanflux[:,None])
    chi = np.nansum(top/char_var, axis=1)
    return chi

def my_chisquare_epoch(flux, sigsq):
    sigsqn = np.tile(sigsq, [len(flux), 1])
    fluxn, sigsqn = vari_funcs.normalise_flux_and_errors(flux, sigsqn)
#    fluxn = vari_funcs.normalise_flux(flux)
    meanflux = np.nanmean(fluxn, axis=1)
    top = np.square(fluxn-meanflux[:,None])
    chi = np.nansum(top/(sigsqn**2), axis=1)
    return chi

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra_clean.fits')[1].data

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
bins, medvar = vari_funcs.plot_median_line(fluxn, tbdata, statistic='var')

bins, medvar = vari_funcs.plot_median_line_stars(fluxn, tbdata, sfluxn, sdata, statistic='var')


vari_funcs.flux_variability_plot(fluxn, fluxchann, 'chisq',
                                   fluxerr = fluxerr,
                                   starfluxerr = sfluxerr,
                                    starflux=sfluxn, stars=True,
                                    chanerr = fluxerrchan,
                                    normalised=True, scale='symlog')
plt.text(5e2, 5e4, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{i}^{2}}}$')
plt.hlines(24.322,2e2,1e7, label='99.9% confidence level', zorder=4)
plt.legend()

### Put characteristic variance into chisquare ###
newmedvar = np.empty(np.shape(medvar))
finalmed = np.empty(np.shape(medvar))
sigsqdict = {}
sigdict = {}
finalchisq = np.array([])
galchisq = np.array([])
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
    chisq = my_chisquare_char(flux, medvar[n])
    chisqchan = my_chisquare_char(fluxchan, medvar[n])
    schisq = my_chisquare_char(sflux, medvar[n])
    
    ### plot ###
    plt.figure(3,figsize=[8,8])    
    if n == np.size(bins)-2:
        plt.plot(meanflux, chisq, 'b+',zorder=2, label='Galaxy')
        plt.plot(meanchan, chisqchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
        plt.plot(meansflux, schisq, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin=2e2, xmax=1e7)
        plt.ylim(ymin=3e-2, ymax=4e4)
        plt.ylabel('Chi Squared')
        plt.xlabel('Mean Flux')
        plt.title('1st iteration')
        plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
        plt.hlines(24.322,2e2,1e7, label='99.9% confidence level', zorder=4)
        plt.legend()
    else:
        plt.plot(meanflux, chisq, 'b+',zorder=2)
        plt.plot(meanchan, chisqchan, 'ro', zorder=3, mfc='None', markersize=10)
        plt.plot(meansflux, schisq, 'm*', zorder=1, mfc='None', markersize=10)

    
    ### remove any with chisq > 70 ###
    newgflux = flux[chisq<50,:]
    newsflux = sflux[schisq<50,:]
    newflux = np.vstack((newgflux, newsflux))
    newfluxn = vari_funcs.normalise_flux(newflux)
    vary = np.var(newfluxn, axis=1, ddof=1)
    newmedvar[n] = np.nanmedian(vary)

    newchisq = my_chisquare_char(flux, newmedvar[n])
    newchisqchan = my_chisquare_char(fluxchan, newmedvar[n])
    newschisq = my_chisquare_char(sflux, newmedvar[n])
    
    ## plot new ###
    plt.figure(4, figsize=[8,8])
    if n == np.size(bins)-2:
        plt.plot(meanflux, newchisq, 'b+',zorder=2, label='Galaxy')
        plt.plot(meanchan, newchisqchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
        plt.plot(meansflux, newschisq, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin=2e2, xmax=1e7)
        plt.ylim(ymin=3e-2, ymax=4e4)
        plt.ylabel('Chi Squared')
        plt.xlabel('Mean Flux')
        plt.title('2nd iteration')
        plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
        plt.hlines(24.322,2e2,1e7, label='99.9% confidence level', zorder=4)
        plt.legend()
    else:
        plt.plot(meanflux, newchisq, 'b+',zorder=2)
        plt.plot(meanchan, newchisqchan, 'ro', zorder=3, mfc='None', markersize=10)
        plt.plot(meansflux, newschisq, 'm*', zorder=1, mfc='None', markersize=10)

        

    ### plot new variance ###
#    plt.figure(5, figsize=[8,8])
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
#    #plot
#    plt.plot(meanflux, sigreal, 'b+', zorder=2)
#    plt.plot(meanchan, sigrealchan, 'ro', zorder=3, mfc='None', markersize=10)
#    plt.plot(meansflux, ssigreal, 'm*', zorder=1, mfc='None', markersize=10)
#    plt.yscale('symlog', linthreshy=0.0001)
#    plt.xscale('log')
#    plt.ylabel(r'$\sigma^{2}_{real}$')
#    plt.xlabel('Mean Flux')
#    plt.text(1e5, 1e0, r'$\sigma^{2}_{real} = \sigma^{2} - \sigma_{noise}^{2}$')
#    plt.tight_layout()
#    finalmed[n] = np.nanmedian(sigreal)

    ### Calculate sigma for each epoch in flux bin ###
    
    avgflux = np.nanmean(newflux, axis=1)
    ### Find sigma^2 ###
    diff = newflux - avgflux[:,None]
    top = diff**2
    bot = len(avgflux)
    sigsqdict[str(int(binedge))] = np.nansum(top/bot, axis=0)
    sigdict[str(int(binedge))] = np.sqrt(sigsqdict[str(int(binedge))])
    
    ### Get another new chi-squared with epoch error ###
    newchisq2 = my_chisquare_epoch(flux, sigdict[str(int(binedge))])
    newchisqchan2 = my_chisquare_epoch(fluxchan, sigdict[str(int(binedge))])
    newschisq2 = my_chisquare_epoch(sflux, sigdict[str(int(binedge))])
    
    ### plot new ###
    plt.figure(6, figsize=[8,8])    
    if n == np.size(bins)-2:
        plt.plot(meanflux, newchisq2, 'b+',zorder=2, label='Galaxy')
        plt.plot(meanchan, newchisqchan2, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
        plt.plot(meansflux, newschisq2, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin=2e2, xmax=1e7)
        plt.ylim(ymin=3e-2, ymax=4e4)
        plt.ylabel('Chi Squared')
        plt.xlabel('Mean Flux')
        plt.title('With epoch flux bin errors')
        plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{epoch}^{2}}}$')
        plt.hlines(24.322,2e2,1e7, label='99.9% confidence level', zorder=4)
        plt.legend()
    else:
        plt.plot(meanflux, newchisq2, 'b+',zorder=2)
        plt.plot(meanchan, newchisqchan2, 'ro', zorder=3, mfc='None', markersize=10)
        plt.plot(meansflux, newschisq2, 'm*', zorder=1, mfc='None', markersize=10)


    ### Save final chisq for stars and gals in 1 large array ###
    tempchi = np.hstack([newchisq2, newschisq2])
    finalchisq = np.append(finalchisq, tempchi)
    galchisq = np.append(galchisq, newchisq2)
    
#plt.plot(bins[0:42], finalmed, 'k--',zorder=3)

histbins = np.logspace(-2.3,4.4,200)
plt.figure()
plt.hist(finalchisq, histbins, normed=True, label=r'Histogram of $\chi^{2}$ for all galaxies')
plt.xscale('log')
#plt.yscale('symlog')
plt.ylabel('Counts')
plt.xlabel('Chi Squared')

x = np.logspace(-2.3,4.4,500)
#x = np.linspace(3e-2,4e4,5000)
y = stats.chi2.pdf(x,7) #7 dof as 8 variables
plt.plot(x,y, label=r'Model $\chi^{2}$ with dof=7')
#plt.vlines(7, 0, 0.12)
plt.legend()

varychi = galchisq[galchisq > 24.322]

#### Turn dictionary into astropy table ###
#t = Table(sigdict)
#t.write('epoch_sigma_table.fits')






