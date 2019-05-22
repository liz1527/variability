#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 10:19:33 2018

Code to determine error on individual epochs as a function of flux.

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
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_1519match.fits')[1].data

### Get flux array and remove values that are negative ###
flux = vari_funcs.flux5_stacks(tbdata)
sflux = vari_funcs.flux5_stacks(sdata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)

### Everything needs to be done on flux bins so create bin array ###
#bins = np.array([13, 15])
#bins = np.append(bins, np.arange(16,22.4,0.2))
bins = np.arange(13,22.4,0.2)
bins = np.append(bins, [24])
bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

### set up time variable ###
t = np.arange(8)
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']

lastind = len(bins)-1
sigsq = np.empty([lastind,8])
sigsqdict = {}
sigdict = {}
oldnum = 0
for n, binedge in enumerate(bins[0:lastind]):
    sbinflux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], sflux, sdata) #bin data
    print(str(binedge) + ' '+ str(len(sbinflux)))
    
    ### find average flux of stars ###
    avgflux = np.nanmean(sbinflux, axis=1)
    
    ### Find sigma^2 ###
    diff = sbinflux - avgflux[:,None]
    top = diff**2
    bot = len(avgflux)
    sigsq[n,:] = np.nansum(top/bot, axis=0)
    sigsqdict[binedge] = np.nansum(top/bot, axis=0)
    sigdict[binedge] = np.sqrt(sigsqdict[binedge])
    
#    ### Plot average lightcurve for bin ###
#    plt.figure()
#    avgperbin = np.nanmean(sbinflux, axis=0)
#    plt.errorbar(t,avgperbin, sigdict[binedge], marker='o', color='k')
#    for m in range(len(sbinflux)):
#        plt.scatter(t, sbinflux[m,:])
#    plt.xticks(t, semesters)
#    plt.title(str(int(binedge))+' < F < '+str(int(bins[n+1])))
#    plt.ylabel('Flux')
#    plt.xlabel('Semester')
#    plt.tight_layout()
#    plt.savefig('Bin_Average_Curves/'+str(int(binedge))+'-'+str(int(bins[n+1]))+'_error_check.png')
#    plt.close()
    plt.figure()
    plt.hist(sbinflux[:,5])
    
#    binflux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], flux, tbdata) #bindata
#    meanbinflux = np.nanmean(binflux, axis=1)
#    fullsig = np.empty(np.shape(binflux))
#    num = len(binflux)
#    fullsig[0:num,:] = sigdict[binedge]
#    binfluxn, binerrn = vari_funcs.normalise_flux_and_errors(binflux, fullsig)
#    vary = vari_funcs.normsigmasq(binfluxn, binerrn)
#    plt.plot(meanbinflux, vary, 'b+', zorder=1)
##    plt.yscale('symlog', limthreshy=1e-7)
##    plt.yscale('log')
#    plt.xscale('log')