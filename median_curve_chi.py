#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:23:10 2018

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

def my_chisquare_quad_median_err(flux, fluxerr, mediancurve):
#    fluxn, fluxerrn = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    top = np.square(flux-mediancurve[None,:])
    bot = np.square(fluxerr)
    chi = np.nansum(top/bot, axis=1)
    return chi

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra_clean.fits')[1].data
sigtb = Table.read('quad_epoch_sigma_table_extra_clean.fits')

## Create arrays of flux values ###
flux = vari_funcs.flux5_stacks(tbdata)
fluxchan = vari_funcs.flux5_stacks(chandata) 
sflux = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Get error arrays ###
quadflux, quaderr, quaddata = vari_funcs.create_quad_error_array(sigtb, tbdata, quadoutput=True)
quadchanflux, quadchanerr, quadchandata = vari_funcs.create_quad_error_array(sigtb, chandata, quadoutput=True)
quadsflux, quadserr, quadsdata = vari_funcs.create_quad_error_array(sigtb, sdata, quadoutput=True)

### get flux bounds ###
bins = np.array(sigtb.colnames)
binarr = np.empty(55)
k = 0
### Create binedge array ###
for bin in bins:
    if bin[0] == '1':
        binarr[k] = int(bin[2:])
        k+=1
plt.figure(figsize=[8,8])      
### reset X-ray column as messed up by stacking ###
for key, qdata in quaddata.items():
    qdata['X-ray'][qdata['X-ray']==70] = False 
    qdata['X-ray'][qdata['X-ray']==84] = True
    
    qflux = quadflux[key]
    qerr = quaderr[key]
    qchanflux = quadchanflux[key]
    qchanerr = quadchanerr[key]
    qsflux = quadsflux[key]
    qserr = quadserr[key]
    
    for n, binedge in enumerate(binarr[0:-1]):
        #get chi for gal
        qbflux, qberr = vari_funcs.fluxbinerr(binedge, binarr[n+1], qflux, qerr)
        meanqb = np.nanmean(qbflux, axis=1)
        qbflux, qberr = vari_funcs.normalise_flux_and_errors(qbflux, qberr)
        mediancurve = np.nanmedian(qbflux, axis=0)
        print(mediancurve)
        
        chisq = my_chisquare_quad_median_err(qbflux, qberr, mediancurve)
        chisqold = vari_funcs.my_chisquare_err(qbflux, qberr)
#        plt.plot(meanqb, chisqold, 'mo', mfc='None', markersize=10)
        
        # get chi for chan
        qbchanflux, qbchanerr = vari_funcs.fluxbinerr(binedge, binarr[n+1], qchanflux, qchanerr)
        meanqbchan = np.nanmean(qbchanflux, axis=1)
        qbchanflux, qbchanerr = vari_funcs.normalise_flux_and_errors(qbchanflux, qbchanerr)
        
        chisqchan = my_chisquare_quad_median_err(qbchanflux, qbchanerr, mediancurve)
        chisqoldchan = vari_funcs.my_chisquare_err(qbchanflux, qbchanerr)
#        plt.plot(meanqb, chisqold, 'mo', mfc='None', markersize=10)
    
        # get chi for stars 
        qbsflux, qbserr = vari_funcs.fluxbinerr(binedge, binarr[n+1], qsflux, qserr)
        meanqbs = np.nanmean(qbsflux, axis=1)
        qbsflux, qbserr = vari_funcs.normalise_flux_and_errors(qbsflux, qbserr)
        mediancurves = np.nanmedian(qbsflux, axis=0)
        
#        if n==0:
        chisqs = my_chisquare_quad_median_err(qbsflux, qbserr, mediancurves)
        chisqolds = vari_funcs.my_chisquare_err(qbsflux, qbserr)
#        plt.plot(meanqb, chisqold, 'mo', mfc='None', markersize=10)
        
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim(ymin=3e-2, ymax=4e4)
        plt.xlim(xmin=2e2, xmax=1e7)
        plt.xlabel('Mean Flux')
        plt.ylabel(r'$\chi^{2}$')
        plt.title('With quad epoch flux bin errors and median curve')
        if key==0 and n==0:
            plt.plot(meanqb, chisq,'b+',zorder=1,label='UDS Source')
            plt.plot(meanqbchan, chisqchan,'ro', mfc='None', markersize=10,zorder=2,label='X-ray Source')
            plt.plot(meanqbs, chisqs,'m*', mfc='None', markersize=10,zorder=0,label='DR11 Star')
            plt.legend()
        else:
            plt.plot(meanqb, chisq,'b+',zorder=1)
            plt.plot(meanqbchan, chisqchan,'ro', mfc='None', markersize=10,zorder=2)
            plt.plot(meanqbs, chisqs,'m*', mfc='None', markersize=10,zorder=0)

            