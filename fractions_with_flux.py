#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:56:43 2019

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots

    
### Get data ###
spuds = fits.open('UDS_catalogues/SpUDS_catalogue_DR11data_mag_flux.fits')[1].data
tbdata =fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
xdata =fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
varyspuds = fits.open('variable_tables/no06_variables_chi30_2arcsec_SpUDSmatch.fits')[1].data
xvaryspuds = fits.open('variable_tables/no06_variables_chi30_2arcsec_SpUDSmatch_xray.fits')[1].data
vary = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
xvary = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data



def find_number(tb):
    ### Get flux ###
    flux = vari_funcs.flux4_stacks(tb) 
    
    bins = np.array([13, 15])
    bins = np.append(bins, np.arange(16,24,0.2))
    bins = np.append(bins, [24])
    
    bins = 10**((30-bins)/2.5)
    bins = np.flip(bins, axis=0)
    #bins = bins[16:44] #because of flux limit
    
    num = np.zeros(len(bins))
    
    ### Bin data ###
    for n, binedge in enumerate(bins):
    #    print(binedge)
        if n==np.size(bins)-1:
            break
        mag, _ = vari_funcs.fluxbin(binedge, bins[n+1], flux, tb)
        num[n] = len(mag)
        
    return num, bins

#spuds = spuds[~np.isnan(spuds['mag_ap_24'])]
spudsnum, bins = find_number(spuds)
tbnum, bins = find_number(tbdata)
xnum, bins = find_number(xdata)
spudsvarynum, bins = find_number(varyspuds)
xvaryspudsnum, bins = find_number(xvaryspuds)
varynum, bins = find_number(vary)
xvarynum, bins = find_number(xvary)

#plt.figure(figsize=[8,6])
#plt.bar(bins, spudsnum, color='g', label='SpUDS variable')
#plt.bar(bins, xspudsnum, color='m', label='X-ray and SpUDS variable')
#plt.bar(bins, xvarynum, color='r', label='X-ray variable')
#plt.bar(bins, varynum, color='b', label='Variable')
#plt.xscale('log')
#plt.xlabel('Flux')
#plt.ylabel('Number')
#plt.legend()
#plt.tight_layout()

### Fractions of things ###
fracspuds = spudsnum/tbnum
fracspudsvary = spudsvarynum/varynum
plt.figure(figsize=[8,6])
plt.plot(bins, fracspuds, label='SpUDS/All UDS')
plt.plot(bins, fracspudsvary, label='Variable SpUDS/All Variable')
plt.xscale('log')
plt.xlabel('Flux')
plt.ylabel('Fraction')
plt.xlim(8e1, 1e7)
plt.legend()
plt.tight_layout()


#%% Fractions of things ###
fracxspuds = xvaryspudsnum/spudsvarynum
fracxvary = xvarynum/varynum
plt.figure(figsize=[8,6])
plt.plot(bins, fracxspuds, label='Variable SpUDS X-ray/Variable SpUDS')
plt.plot(bins, fracxvary, label='Variable X-ray/Variable ')
plt.xscale('log')
plt.xlabel('Flux')
plt.ylabel('Fraction')
plt.xlim(8e1, 1e7)
plt.ylim(ymax=1.2)
plt.legend(loc='upper left')
plt.tight_layout()