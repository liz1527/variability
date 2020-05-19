#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

Code to identify how many objects are within a certain radius of the variables

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table, Column #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy import units as u
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

    
### Open the fits files and get data ###
tbdata1 = Table.read('variable_tables/matched_Xray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')
tbdata2 = Table.read('variable_tables/matched_notXray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')
tbdata3 = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')

### restrict redshift ###
tbdata1 = tbdata1[tbdata1['z']<1]
tbdata2 = tbdata2[tbdata2['z']<1]
tbdata3 = tbdata3[tbdata3['z']<1]

#%% match variables and fullcat catalogue ###

### Set match distance ###
arcsec_dist = 5

def find_nearby(varys, arcsec_dist):
    ''' find number of UDS objects within a certain distance from each variable '''
    dist = arcsec_dist * 7.45 # convert dist in pixels 
    num = np.zeros(len(varys))
    
    ### Open full catalogue ###
    fullcat = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019.fits')

    ### Get coordinates ###
    x_fullcat = fullcat['X_IMAGE']
    y_fullcat = fullcat['Y_IMAGE']
    
    for n, id in enumerate(varys['ID']):
        obdata = varys[varys['ID']==id]
        x_varys = obdata['X_IMAGE']
        y_varys = obdata['Y_IMAGE']
        
        ### define radial distance from variable to all x talks ###
        r = np.sqrt(np.square(x_varys-x_fullcat) + np.square(y_varys-y_fullcat))
        
        ### identify any fullcats within 5 arcsec ###
        mask = r <= dist
        near_vary = x_fullcat[mask]
        
        num[n] = len(near_vary)
    return num-1 #as will always find the variable object in the catalogue

num1 = find_nearby(tbdata1, arcsec_dist)
num2 = find_nearby(tbdata2, arcsec_dist)
num3 = find_nearby(tbdata3, arcsec_dist)

num1_mean = np.nanmean(num1)
num2_mean = np.nanmean(num2)
num1_std = np.nanstd(num1)
num2_std = np.nanstd(num2)

### Split by X-ray ###
num3_xray = num3[tbdata3['X-ray']==True]
num3_notxray = num3[tbdata3['X-ray']==False]

num3_xray_mean = np.nanmean(num3_xray)
num3_xray_std = np.nanstd(num3_xray)
num3_notxray_mean = np.nanmean(num3_notxray)
num3_notxray_std = np.nanstd(num3_notxray)

### Get luminosity values ###
tb3_lum = tbdata3['M_K_z_p']
tb3_xray_lum = tb3_lum[tbdata3['X-ray']==True]
tb3_notxray_lum = tb3_lum[tbdata3['X-ray']==False]

tb1_lum = tbdata1['M_K_z_p']
tb2_lum = tbdata2['M_K_z_p']

#%% Make plots ###
plt.figure()
plt.scatter(tb2_lum, num2, c='b', label='Not X-ray')
plt.scatter(tb1_lum, num1, c='r', label='X-ray')
plt.xlabel('K-Band Absolute Magitude')
plt.ylabel('Num within 5"')
plt.xlim(-27, -14)
plt.title('Matched Variable Sample')
plt.legend()
plt.tight_layout()


plt.figure()
plt.scatter(tb3_notxray_lum, num3_notxray, c='b', label='Not X-ray')
plt.scatter(tb3_xray_lum, num3_xray, c='r', label='X-ray')
plt.xlabel('K-Band Absolute Magitude')
plt.ylabel('Num within 5"')
plt.title('Full Variable Sample')
plt.xlim(-27, -14)
plt.legend()
plt.tight_layout()



















