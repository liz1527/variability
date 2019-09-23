#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:25:12 2019

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
plt.close('all')

#%% histograms of z ###

varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data

varydatalow = vari_funcs.flux_split(varydata, 'lower')
varydataup = vari_funcs.flux_split(varydata, 'upper')

zlow = vari_funcs.get_z(varydatalow)
zup = vari_funcs.get_z(varydataup)

plt.hist([zlow, zup], color=['m','g'],histtype='step', 
         label=[r'$\bar{F} < 8 \times 10^{3}$',r'$\bar{F} \geq 8 \times 10^{3}$'],
         normed=True)
plt.legend()
plt.xlabel('z')
plt.ylabel('Number')

### with x/nox split ###
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data

noxz = vari_funcs.get_z(noxvarydata)
xz = vari_funcs.get_z(xvarydata)

plt.figure()
plt.hist([noxz, xz], color=['b','r'],histtype='step', 
         label=['Non X-ray',r'X-ray'])#,
#         normed=True)
plt.legend()
plt.xlabel('z')
plt.ylabel('Number')

#%% histograms of flux ###

### Get fits ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
dr11 = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best.fits')[1].data
### get fluxes from mag-flux ###
def get_mean_flux(tbdata):
    flux = vari_funcs.flux4_stacks(tbdata)
    meanflux = np.nanmean(flux, axis=1)
    return meanflux

def get_jansky_flux(tbdata):
    meanmag = tbdata['KMAG_20']
    meanflux = 10**(23-((meanmag+48.6)/2.5))
    return meanflux
    
meanflux = get_mean_flux(tbdata)
meannoxvary = get_mean_flux(noxvarydata)
meanxvary = get_mean_flux(xvarydata)
meanvary = get_mean_flux(varydata)

### get fluxs from DR11data ###
#meanflux = get_jansky_flux(dr11)
#dr11vary = get_jansky_flux(varydata)

### plot ###
#bins = np.logspace(-8,-3,50)
#bins = np.logspace(0,7,50)    
bins = np.arange(13,25,0.2)
bins = np.append(bins, [25])
bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

plt.figure()
plt.hist(meanflux,  bins, histtype='step', normed='True')
plt.hist(meannoxvary,  bins, histtype='step', color='b', normed='True')
plt.xscale('log')
plt.xlabel('2" K band flux')
plt.ylabel('Number')
plt.tight_layout()

plt.figure()
#plt.hist(meanvary,  bins, histtype='step', color='k')#, normed='True')
plt.hist(meannoxvary,  bins, histtype='step', color='b', label='Non X-ray')#, normed='True')
plt.hist(meanxvary,  bins, histtype='step', color='r', label='X-ray')#, normed='True')
plt.xscale('log')
plt.xlabel('2" K band flux')
plt.ylabel('Number')
plt.legend()
plt.tight_layout()