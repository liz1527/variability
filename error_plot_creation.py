#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:19:08 2018

Code to find sigma in quadrant, epoch and flux bins

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
from scipy.stats import chisquare
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
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('mag_flux_tables/K/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
sdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_extra_clean_no06.fits')[1].data

### Split the data into the 4 quadrants ###
quaddata = vari_funcs.quadrants(tbdata, '05B')
chanquaddata = vari_funcs.quadrants(chandata, '05B')
squaddata = vari_funcs.quadrants(sdata, '05B')

sigsqdict = {}
sigdict = {}
finalchisq = np.array([])
galchisq = np.array([])
oldchisq = np.array([])
for m, qtbdata in enumerate(quaddata):
    ### get quad arrays for chan and stars ###
    qchandata = chanquaddata[m]
    qsdata = squaddata[m]
    
    ## Create arrays of flux values ###
    fluxn = vari_funcs.flux4_stacks(qtbdata)
    fluxchann = vari_funcs.flux4_stacks(qchandata) 
    sfluxn = vari_funcs.flux4_stacks(qsdata)
    
    ### remove values that are negative ###
    fluxn, qtbdata = vari_funcs.noneg(fluxn, qtbdata)
    fluxchann, qchandata = vari_funcs.noneg(fluxchann, qchandata)
    sfluxn, qsdata = vari_funcs.noneg(sfluxn, qsdata)
    
    fluxerr = vari_funcs.fluxerr4_stacks(qtbdata)
    fluxerrchan = vari_funcs.fluxerr4_stacks(qchandata)
    sfluxerr = vari_funcs.fluxerr4_stacks(qsdata)
    
    print(len(fluxn)+len(sfluxn))
    
    if m == 3:
        ### plot flux vs err ###
        plt.figure(1)
#        plt.subplot(2,2,m+1)
        plt.plot(fluxn[:,6], fluxerr[:,6], 'b+')
        plt.plot(sfluxn[:,6], sfluxerr[:,6], 'm*')
    
    ### get original chi sq ###
    chisq = vari_funcs.my_chisquare_err(fluxn, fluxerr)
    schisq = vari_funcs.my_chisquare_err(sfluxn, sfluxerr)
    
    allchi = np.append(chisq, schisq)
    oldchisq = np.append(oldchisq, allchi)
    
#    vari_funcs.flux_variability_plot(fluxn, fluxchann, 'var', starflux=sfluxn, 
#                                     stars=True, normalised=True, scale='symlog')
    #bins, medvar = vari_funcs.plot_median_line(fluxn, tbdata, statistic='var')
    
    bins, medvar = vari_funcs.plot_median_line_stars(fluxn, qtbdata, sfluxn, 
                                                     qsdata, statistic='var',
                                                     createplot=False)
#    plt.close(2)
    
#    vari_funcs.flux_variability_plot(fluxn, fluxchann, 'chisq',
#                                       fluxerr = fluxerr,
#                                       starfluxerr = sfluxerr,
#                                        starflux=sfluxn, stars=True,
#                                        chanerr = fluxerrchan,
#                                        normalised=True, scale='symlog')
#    plt.text(5e2, 5e4, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{i}^{2}}}$')
#    plt.hlines(24.322,2e2,1e7, label='99.9% confidence level', zorder=4)
#    plt.legend()
    
    ### Put characteristic variance into chisquare ###
    newmedvar = np.empty(np.shape(medvar))
    finalmed = np.empty(np.shape(medvar))
    for n, binedge in enumerate(bins[0:np.size(bins)-1]):
        flux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], fluxn, qtbdata) #bindata
        fluxchan, binchan = vari_funcs.fluxbin(binedge, bins[n+1], fluxchann, qchandata) #bindata
        sflux, sbindata = vari_funcs.fluxbin(binedge, bins[n+1], sfluxn, qsdata) #bindata
        # get errors
        fluxerr = vari_funcs.fluxerr4_stacks(bindata)
        fluxchanerr = vari_funcs.fluxerr4_stacks(binchan)
        sfluxerr = vari_funcs.fluxerr4_stacks(sbindata)
        
        medfluxerr = np.nanmedian(np.append(fluxerr[:,6], sfluxerr[:,6],))
        
#        print(len(flux))
        meanflux = np.nanmean(flux, axis=1)
        meanchan = np.nanmean(fluxchan, axis=1)
        meansflux = np.nanmean(sflux, axis=1)
        chisq = my_chisquare_char(flux, medvar[n])
        chisqchan = my_chisquare_char(fluxchan, medvar[n])
        schisq = my_chisquare_char(sflux, medvar[n])
        
        
        ### remove any with chisq > 50 ###
        newgflux = flux[chisq<50,:]
        newsflux = sflux[schisq<50,:]
        newflux = np.vstack((newgflux, newsflux))
        newfluxn = vari_funcs.normalise_flux(newflux)
        vary = np.var(newfluxn, axis=1, ddof=1)
        newmedvar[n] = np.nanmedian(vary)
    
#        print(len(newflux))
        
        newchisq = my_chisquare_char(flux, newmedvar[n])
        newchisqchan = my_chisquare_char(fluxchan, newmedvar[n])
        newschisq = my_chisquare_char(sflux, newmedvar[n])
        
        
    
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
        
        ### Calculate sigma for each epoch in flux bin ###

        avgflux = np.nanmean(newflux, axis=1)
        
        ### Find sigma^2 ###
        diff = newflux - avgflux[:,None]
        top = diff**2
        bot = len(avgflux)
        sig = np.nansum(top/bot, axis=0)
        dictkey = str(m+1)+' '+str(int(binedge))
        sigsqdict[dictkey] = np.nansum(top/bot, axis=0)
        sigdict[dictkey] = np.sqrt(sigsqdict[dictkey])
        
        if m == 3:
            ### Add point to plot ###
            plt.figure(1)
#            plt.subplot(2,2,m+1)
            plt.plot(int(binedge), sig[6], 'go', markersize=7)
            plt.yscale('log')
            plt.xscale('log')
            plt.xlabel('Flux')
            plt.ylabel(r'$\sigma_{F}$')
        
        ### Get another new chi-squared with epoch error ###
        newchisq2 = my_chisquare_epoch(flux, sigdict[dictkey])
        newchisqchan2 = my_chisquare_epoch(fluxchan, sigdict[dictkey])
        newschisq2 = my_chisquare_epoch(sflux, sigdict[dictkey])
    
        ### Save final chisq for stars and gals in 1 large array ###
        tempchi = np.hstack([newchisq2, newschisq2])
        finalchisq = np.append(finalchisq, tempchi)
        galchisq = np.append(galchisq, newchisq2)
    
        if m == 3 and n == 10:
            ### Plot distribution of flux values ###
            plt.figure(2)
            plt.hist(newflux[:,6], bins=30, density=True, 
                     label='Actual values', color='k')
            plt.xlabel('Flux')
            plt.ylabel('Normalised Frequency')
            
            print(len(newflux)) # print how many sources included
            
            print(np.nanmedian(newflux[:,6]))
            
            ### Plot guassian over distribution for our sig ###
            gaussig = np.sqrt(sig[6])
            f = np.linspace(np.nanmin(newflux[:,6]), np.nanmax(newflux[:,6]), 100)
            gaus = stats.norm.pdf(f, np.nanmedian(newflux[:,6]), gaussig)
            print(gaussig)
#            gaus = gaus/np.nansum(gaus)
            plt.plot(f, gaus, 'r-', label='Self-calibrated')
            
            ### Plot guassian over distribution for median SE sig ###
            gausSE = stats.norm.pdf(f, np.nanmedian(newflux[:,6]), medfluxerr)
            print(medfluxerr)
#            gausSE = gausSE/np.nansum(gausSE)
            plt.plot(f, gausSE, 'b--', label='SExtractor')
            
            ### Add text saying what the two sigmas are (maybe just put in caption ###
#            plt.text(100, 0.006, r'Median SExtractor $\sigma_{F}$')
#            plt.text(100, 0.0057, '= '+str(medfluxerr))
#            plt.text(100, 0.005, r'Self-calibrated $\sigma_{F}$')
#            plt.text(100, 0.0047, '= '+str(gaussig))
            
#            plt.xlim(xmax=1400)
            plt.legend(loc='upper right')
            plt.tight_layout()
#plt.plot(bins[0:42], finalmed, 'k--',zorder=3)
#%%
histbins = np.logspace(-2.3,4.4,100)
plt.figure()
plt.hist(finalchisq, histbins, label='Self-calibrated uncertainties')
plt.hist(oldchisq, histbins, label='SExtractor Uncertainties')
plt.xscale('log')
#plt.yscale('symlog')
plt.ylabel('Normalised Frequency')
plt.xlabel(r'$\chi^{2}$ ')
plt.tight_layout()

x = np.logspace(-2.3,4.4,500)
#x = np.linspace(3e-2,4e4,5000)
y = stats.chi2.pdf(x,6) #6 dof as 7 variables
#plt.plot(x,y, label=r'Model')
#plt.vlines(7, 0, 0.12)
plt.legend()

varychi = galchisq[galchisq > 24.322]

### Turn dictionary into astropy table ###
t = Table(sigdict)
#t.write('quad_epoch_sigma_table_extra_clean_no06_2arcsec_neg.fits')




