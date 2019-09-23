#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:42:06 2019

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

### Import tables ###
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
fullxraydata = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_xray_novarys.fits')
fullxraymag = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06_novarys.fits')[1].data
fullmag = fits.open('mag_flux_tables/mag_flux_table_extra_clean_no06.fits')[1].data
#fullmag = fits.open('DR11_raw_output.fits')[1].data

### with x/nox split ###
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data

### get mean fluxes in 2 apertures ###
def get_mean_mag_1(tbdata):
    mag = vari_funcs.mag1_stacks(tbdata)
    meanmag = np.nanmean(mag, axis=1)
    return meanmag

def get_mean_mag_4(tbdata):
    mag = vari_funcs.mag2_stacks(tbdata)
    meanmag = np.nanmean(mag, axis=1)
    return meanmag

meanmag1 = get_mean_mag_1(fullmag)
#meanmag1 = fullmag['MAG_APER'][:,0]
xmeanmag1 = get_mean_mag_1(fullxraymag)
xvarymeanmag1 = get_mean_mag_1(xvarydata)
noxvarymeanmag1 = get_mean_mag_1(noxvarydata)

meanmag4 = get_mean_mag_4(fullmag)
#meanmag4 = fullmag['MAG_APER'][:,1]
xmeanmag4 = get_mean_mag_4(fullxraymag)
xvarymeanmag4 = get_mean_mag_4(xvarydata)
noxvarymeanmag4 = get_mean_mag_4(noxvarydata)

### find differences ###
mag4_mag1 = meanmag4 - meanmag1
xmag4_mag1 = xmeanmag4 - xmeanmag1
xvarymag4_mag1 = xvarymeanmag4 - xvarymeanmag1
noxvarymag4_mag1 = noxvarymeanmag4 - noxvarymeanmag1


### Plot ###
plt.figure(figsize=[8,8])
plt.plot(meanmag4, mag4_mag1, '.', markersize=1, color='tab:grey', alpha=0.35, 
         label='UDS Source')
plt.plot(xmeanmag4, xmag4_mag1, 'k+', label='X-ray Non-Variable')
plt.plot(noxvarymeanmag4, noxvarymag4_mag1, 'bo', label='X-ray Variable')
plt.plot(xvarymeanmag4, xvarymag4_mag1, 'ro', label='Non-X-ray Variable')
#plt.xlim(xmin=15,xmax=25)
#plt.ylim(ymin=-2.2,ymax=-0.5)
plt.xlim(xmin=15,xmax=25)
plt.ylim(ymin=-1,ymax=-0.3)
plt.xlabel('K band 1" magnitude')
plt.ylabel('K band 1" magnitude - K band 0.7" magnitude')
plt.legend(loc='upper left')