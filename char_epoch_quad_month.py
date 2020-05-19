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
#from scipy.stats import chisquare
plt.close('all') #close any open plots

import time
start = time.time()
print(start)

def compute_p(flux, chisq):
    p = np.empty(len(chisq))
    for n, chi in enumerate(chisq):
        countarr = np.ones(len(flux[n,:]))
        countarr[np.isnan(flux[n,:])] = 0 # set nans to 1
        dof = np.sum(countarr) #should count number of not nan values
        dof -= 1
        p[n] = 1 - stats.chi2.cdf(chi, dof)
    return p

def get_quaddata(filename):
    ### Open the fits files and get data ###
#    tbdata = fits.open(filename)[1].data
    tbdata = Table.read(filename)
    print('Removing edges')
    ### Remove edges ###
    tbdata = vari_funcs.field_funcs.remove_edges(tbdata, 'sep05')
    print('Removed edges')
    ### Split the data into the 4 quadrants ###
    quaddata = vari_funcs.field_funcs.quadrants(tbdata, 'sep05')
    
    ### Delete full table to save space ###
    del tbdata
    return quaddata

filename = 'mag_flux_tables/K/month/month_mag_flux_table_best_K_extra_quad_clean_38.fits'
chanfilename = 'mag_flux_tables/K/month/month_xray_mag_flux_table_best_K_extra_quad_clean_38.fits'
sfilename = 'mag_flux_tables/K/month/month_stars_mag_flux_table_K_extra_quad_clean_38.fits'

quaddata = get_quaddata(filename)
chanquaddata = get_quaddata(chanfilename)
squaddata = get_quaddata(sfilename)

ap = 4 # set up aperture size

#### Remove edges ###
#tbdata = vari_funcs.field_funcs.remove_edges(tbdata, 'sep05')
#chandata = vari_funcs.field_funcs.remove_edges(chandata, 'sep05')
#sdata = vari_funcs.field_funcs.remove_edges(sdata, 'sep05')
#
#### Split the data into the 4 quadrants ###
#quaddata = vari_funcs.field_funcs.quadrants(tbdata, 'sep05')
#chanquaddata = vari_funcs.field_funcs.quadrants(chandata, 'sep05')
#squaddata = vari_funcs.field_funcs.quadrants(sdata, 'sep05')
#
#### delete full tables to save memory ###
#del tbdata, chandata, sdata

sigsqdict = {}
sigdict = {}
finalchisq = np.array([])
galchisq = np.array([])
oldchisq = np.array([])
lost=0
for m, qtbdata in enumerate(quaddata):
    ### get quad arrays for chan and stars ###
    qchandata = chanquaddata[m]
    qsdata = squaddata[m]
    print(len(qtbdata)+len(qsdata))
    
    ## Create arrays of flux values ###
    fluxn = vari_funcs.k_mag_flux.month_flux_stacks(qtbdata, aper=ap)
    fluxchann = vari_funcs.k_mag_flux.month_flux_stacks(qchandata, aper=ap) 
    sfluxn = vari_funcs.k_mag_flux.month_flux_stacks(qsdata, aper=ap)
    
#    ### remove values that are negative ###
#    fluxn, qtbdata = vari_funcs.flux_funcs.noneg(fluxn, qtbdata)
#    fluxchann, qchandata = vari_funcs.flux_funcs.noneg(fluxchann, qchandata)
#    sfluxn, qsdata = vari_funcs.flux_funcs.noneg(sfluxn, qsdata)
    
    fluxerr = vari_funcs.k_mag_flux.month_fluxerr_stacks(qtbdata, aper=ap)
    fluxerrchan = vari_funcs.k_mag_flux.month_fluxerr_stacks(qchandata, aper=ap)
    sfluxerr = vari_funcs.k_mag_flux.month_fluxerr_stacks(qsdata, aper=ap)
    
#    ## nan values that are negative ###
#    fluxn, fluxerr, qtbdata = vari_funcs.flux_funcs.nanneg(fluxn, fluxerr, qtbdata)
#    fluxchann, fluxerrchan, qchandata = vari_funcs.flux_funcs.nanneg(fluxchann, 
#                                                                     fluxerrchan,
#                                                                     qchandata)
#    sfluxn, sfluxerr, qsdata = vari_funcs.flux_funcs.nanneg(sfluxn, sfluxerr, qsdata)
    
    
#    ### plot flux vs err ###
#    plt.figure(11)
#    plt.subplot(2,2,m+1)
#    plt.plot(fluxn[:,0], fluxerr[:,0], 'b+')
#    plt.plot(sfluxn[:,0], sfluxerr[:,0], 'm*')
    
#    vari_funcs.selection_plot_funcs.flux_variability_plot(fluxn, fluxchann, 'var', starflux=sfluxn, 
#                                     stars=True, normalised=True, scale='symlog')
    #bins, medvar = vari_funcs.selection_plot_funcs.plot_median_line(fluxn, tbdata, statistic='var')
    
    bins, medvar = vari_funcs.selection_plot_funcs.plot_median_line_stars(fluxn, 
                                                                          qtbdata, 
                                                                          sfluxn, 
                                                                          qsdata, 
                                                                          statistic='var',
                                                                          createplot=False)
    
    
#    vari_funcs.selection_plot_funcs.flux_variability_plot(fluxn, fluxchann, 'chisq',
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
        flux, bindata = vari_funcs.flux_funcs.fluxbin(binedge, bins[n+1], fluxn, qtbdata) #bindata
        fluxchan, binchan = vari_funcs.flux_funcs.fluxbin(binedge, bins[n+1], fluxchann, qchandata) #bindata
        sflux, sbindata = vari_funcs.flux_funcs.fluxbin(binedge, bins[n+1], sfluxn, qsdata) #bindata
        # get errors
        fluxerr = vari_funcs.k_mag_flux.month_fluxerr_stacks(bindata, aper=ap)
        fluxchanerr = vari_funcs.k_mag_flux.month_fluxerr_stacks(binchan, aper=ap)
        sfluxerr = vari_funcs.k_mag_flux.month_fluxerr_stacks(sbindata, aper=ap)
        
#        ### nan values that are negative ###
#        flux, fluxerr, bindata = vari_funcs.flux_funcs.nanneg(flux, fluxerr, 
#                                                              bindata)
#        fluxchan, fluxchanerr, binchan = vari_funcs.flux_funcs.nanneg(fluxchan, 
#                                                                      fluxchanerr, 
#                                                                      binchan)
#        sflux, sfluxerr, sbindata = vari_funcs.flux_funcs.nanneg(sflux, sfluxerr, 
#                                                                 sbindata)
        
#        print(len(flux)+len(sflux))
        meanflux = np.nanmean(flux, axis=1)
        meanchan = np.nanmean(fluxchan, axis=1)
        meansflux = np.nanmean(sflux, axis=1)
        chisq = vari_funcs.vary_stats.my_chisquare_char(flux, medvar[n])
        chisqchan = vari_funcs.vary_stats.my_chisquare_char(fluxchan, medvar[n])
        schisq = vari_funcs.vary_stats.my_chisquare_char(sflux, medvar[n])
        p = compute_p(flux, chisq)
        pchan = compute_p(fluxchan, chisqchan)
        sp = compute_p(sflux, schisq)
        
        ### plot ###
#        plt.figure(3,figsize=[8,8])    
#        if n == np.size(bins)-2 and m==3:
#            plt.plot(meanflux, chisq, 'b+',zorder=2, label='Galaxy')
#            plt.plot(meanchan, chisqchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
#            plt.plot(meansflux, schisq, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
#            plt.yscale('log')
#            plt.xscale('log')
#            plt.xlim(xmin=8e1, xmax=1e7)
##            plt.ylim(ymin=3e-2, ymax=4e4)
#            plt.ylabel('Chi Squared')
#            plt.xlabel('Mean Flux')
#            plt.title('1st iteration')
##            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
#            plt.hlines(0.00003931,2e2,1e7, label='99.9% confidence level', zorder=4)
#            plt.legend()
#        else:
#            plt.plot(meanflux, chisq, 'b+',zorder=2)
#            plt.plot(meanchan, chisqchan, 'ro', zorder=3, mfc='None', markersize=10)
#            plt.plot(meansflux, schisq, 'm*', zorder=1, mfc='None', markersize=10)
    
#        plt.figure(3,figsize=[8,8])    
#        if n == np.size(bins)-2 and m==3:
#            plt.plot(meanflux, p, 'b+',zorder=2, label='Galaxy')
#            plt.plot(meanchan, pchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
#            plt.plot(meansflux, sp, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
#            plt.yscale('log')
#            plt.xscale('log')
#            plt.xlim(xmin=8e1, xmax=1e7)
##            plt.ylim(ymin=3e-2, ymax=4e4)
#            plt.ylabel('Chi Squared p value')
#            plt.xlabel('Mean Flux')
#            plt.title('1st iteration')
##            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
#            plt.hlines(0.00003931,2e2,1e7, label='99.9% confidence level', zorder=4)
#            plt.legend()
#        else:
#            plt.plot(meanflux, p, 'b+',zorder=2)
#            plt.plot(meanchan, pchan, 'ro', zorder=3, mfc='None', markersize=10)
#            plt.plot(meansflux, sp, 'm*', zorder=1, mfc='None', markersize=10)
        
        ### remove any that are significantly variable ###
        newgflux = flux[chisq<84]#[p>0.00003931,:]
        newsflux = sflux[schisq<84]#[sp>0.00003931,:]
        newflux = np.vstack((newgflux, newsflux))
        print(len(newflux))
        lost += (len(flux)+len(sflux)) - len(newflux)
        newfluxn = vari_funcs.flux_funcs.normalise_flux(newflux)
        vary = np.nanvar(newfluxn, axis=1, ddof=1)
        newmedvar[n] = np.nanmedian(vary)
    
        newchisq = vari_funcs.vary_stats.my_chisquare_char(flux, newmedvar[n])
        newchisqchan = vari_funcs.vary_stats.my_chisquare_char(fluxchan, newmedvar[n])
        newschisq = vari_funcs.vary_stats.my_chisquare_char(sflux, newmedvar[n])
        newp = compute_p(flux, newchisq)
        newpchan = compute_p(flux, newchisqchan)
        newsp = compute_p(sflux, newschisq)
        
#        ## plot new ###
#        plt.figure(4, figsize=[8,8])
#        if n == np.size(bins)-2 and m==3:
#            plt.plot(meanflux, newchisq, 'b+',zorder=2, label='Galaxy')
#            plt.plot(meanchan, newchisqchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
#            plt.plot(meansflux, newschisq, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
#            plt.yscale('log')
#            plt.xscale('log')
#            plt.xlim(xmin=8e1, xmax=1e7)
#            plt.ylim(ymin=3e-2, ymax=4e4)
#            plt.ylabel('Chi Squared')
#            plt.xlabel('Mean Flux')
#            plt.title('2nd iteration')
#            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
#            plt.hlines(84,2e2,1e7, label='99.9% confidence level', zorder=4)
#            plt.legend()
#        else:
#            plt.plot(meanflux, newchisq, 'b+',zorder=2)
#            plt.plot(meanchan, newchisqchan, 'ro', zorder=3, mfc='None', markersize=10)
#            plt.plot(meansflux, newschisq, 'm*', zorder=1, mfc='None', markersize=10)
    
#        plt.figure(4, figsize=[8,8])
#        if n == np.size(bins)-2 and m==3:
#            plt.plot(meanflux, newp, 'b+',zorder=2, label='Galaxy')
#            plt.plot(meanchan, newpchan, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
#            plt.plot(meansflux, newsp, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
#            plt.yscale('log')
#            plt.xscale('log')
#            plt.xlim(xmin=8e1, xmax=1e7)
##            plt.ylim(ymin=3e-2, ymax=4e4)
#            plt.ylabel('Chi Squared P value')
#            plt.xlabel('Mean Flux')
#            plt.title('2nd iteration')
#            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{noise}^{2}}}$')
#            plt.hlines(0.00003931,2e2,1e7, label='99.9% confidence level', zorder=4)
#            plt.legend()
#        else:
#            plt.plot(meanflux, newp, 'b+',zorder=2)
#            plt.plot(meanchan, newpchan, 'ro', zorder=3, mfc='None', markersize=10)
#            plt.plot(meansflux, newsp, 'm*', zorder=1, mfc='None', markersize=10)
            
            
#        ### remove any that are significantly variable ###
#        newgflux = flux[newp>0.00003931,:]
#        newsflux = sflux[newsp>0.00003931,:]
#        newflux = np.vstack((newgflux, newsflux))
#        print(len(newflux))
#        lost += (len(flux)+len(sflux)) - len(newflux)
#        newfluxn = vari_funcs.flux_funcs.normalise_flux(newflux)
#        vary = np.nanvar(newfluxn, axis=1, ddof=1)
#        newmedvar[n] = np.nanmedian(vary)
#    
#    
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
        
#        ### Add point to plot ###
#        plt.figure(11)
#        plt.subplot(2,2,m+1)
#        plt.plot(int(binedge), sig[0], 'go', markersize=10)
#        plt.yscale('log')
#        plt.xscale('log')
#        
        ### Get another new chi-squared with epoch error ###
        newchisq2 = vari_funcs.vary_stats.my_chisquare_epoch(flux, sigdict[dictkey])
        newchisqchan2 = vari_funcs.vary_stats.my_chisquare_epoch(fluxchan, sigdict[dictkey])
        newschisq2 = vari_funcs.vary_stats.my_chisquare_epoch(sflux, sigdict[dictkey])
        newp2 = compute_p(flux, newchisq2)
        newpchan2 = compute_p(flux, newchisqchan2)
        newsp2 = compute_p(sflux, newschisq2)
        
#        ### plot new ###
        plt.figure(6, figsize=[8,8])    
        if n == np.size(bins)-2 and m==3:
            plt.plot(meanflux, newchisq2, 'b+',zorder=2, label='Galaxy')
            plt.plot(meanchan, newchisqchan2, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
            plt.plot(meansflux, newschisq2, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(3e-2,3e4)
            plt.xlim(8e1, 1e7)
#            plt.xlim(xmin=8e1, xmax=1e7)
#            plt.ylim(ymin=3e-2, ymax=4e4)
            plt.ylabel('Chi Squared')
            plt.xlabel('Mean Flux')
            plt.title('With quad epoch flux bin errors')
#            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{quad-epoch}^{2}}}$')
#            plt.hlines(22.458,4e-1,1e7, label='99.9% confidence level', zorder=4)
            plt.hlines(84,4e-1,1e7, label='Chi=30', zorder=4)
            plt.legend()
        else:
            plt.plot(meanflux, newchisq2, 'b+',zorder=2)
            plt.plot(meanchan, newchisqchan2, 'ro', zorder=3, mfc='None', markersize=10)
            plt.plot(meansflux, newschisq2, 'm*', zorder=1, mfc='None', markersize=10)
#    
        plt.figure(7, figsize=[8,8])    
        if n == np.size(bins)-2 and m==3:
            plt.plot(meanflux, newp2, 'b+',zorder=2, label='Galaxy')
            plt.plot(meanchan, newpchan2, 'ro', zorder=3, mfc='None', markersize=10, label='X-ray detected')
            plt.plot(meansflux, newsp2, 'm*', zorder=1, mfc='None', markersize=10, label='DR11 Star')
            plt.yscale('log')
            plt.xscale('log')
#            plt.ylim(3e-2,3e4)
            plt.xlim(8e1, 1e7)
#            plt.xlim(xmin=8e1, xmax=1e7)
#            plt.ylim(ymin=3e-2, ymax=4e4)
            plt.ylabel('Chi Squared P value')
            plt.xlabel('Mean Flux')
            plt.title('With quad epoch flux bin errors')
#            plt.text(5e2, 1e3, r'$\chi^{2} = \sum{\frac{( \,{x_{i} - \bar{x}})^{2} \,}{\sigma_{quad-epoch}^{2}}}$')
#            plt.hlines(22.458,4e-1,1e7, label='99.9% confidence level', zorder=4)
            plt.hlines(0.00003931,4e-1,1e7, label='Chi=30', zorder=4)
            plt.legend()
        else:
            plt.plot(meanflux, newp2, 'b+',zorder=2)
            plt.plot(meanchan, newpchan2, 'ro', zorder=3, mfc='None', markersize=10)
            plt.plot(meansflux, newsp2, 'm*', zorder=1, mfc='None', markersize=10)
#    
    
    
        ### Save final chisq for stars and gals in 1 large array ###
        tempchi = np.hstack([newchisq2, newschisq2])
        finalchisq = np.append(finalchisq, tempchi)
        galchisq = np.append(galchisq, newchisq2)
    
#plt.plot(bins[0:42], finalmed, 'k--',zorder=3)
#%%
#histbins = np.logspace(-2.3,4.4,100)
#plt.figure()
#plt.hist(finalchisq, histbins, label='Self-calibrated uncertainties')
#plt.hist(oldchisq, histbins, label='SExtractor Uncertainties')
#plt.xscale('log')
##plt.yscale('symlog')
#plt.ylabel('Normalised Counts')
#plt.xlabel(r'$\chi^{2}$ ')
#plt.tight_layout()
#
#x = np.logspace(-2.3,4.4,500)
##x = np.linspace(3e-2,4e4,5000)
#y = stats.chi2.pdf(x,6) #6 dof as 7 variables
##plt.plot(x,y, label=r'Model')
##plt.vlines(7, 0, 0.12)
#plt.legend()

#varychi = galchisq[galchisq > 84]

# Turn dictionary into astropy table ###
#t = Table(sigdict) allmonths.fits', overwrite=True)


end = time.time()
print(end-start)


