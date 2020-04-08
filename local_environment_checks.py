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
tbdata1 = Table.read('variable_tables/K/variables_no06_chi30_neg_DR11data.fits')
tbdata2 = Table.read('variable_tables/K/variables_no06_chi30_1arcsec_neg_DR11data.fits')

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


num1_mean = np.nanmean(num1)
num2_mean = np.nanmean(num2)
num1_std = np.nanstd(num1)
num2_std = np.nanstd(num2)

### Split by X-ray ###
num1_xray = num1[tbdata1['X-ray']==True]
num1_notxray = num1[tbdata1['X-ray']==False]
num2_xray = num2[tbdata2['X-ray']==True]
num2_notxray = num2[tbdata2['X-ray']==False]

num1_xray_mean = np.nanmean(num1_xray)
num2_xray_mean = np.nanmean(num2_xray)
num1_xray_std = np.nanstd(num1_xray)
num2_xray_std = np.nanstd(num2_xray)
num1_notxray_mean = np.nanmean(num1_notxray)
num2_notxray_mean = np.nanmean(num2_notxray)
num1_notxray_std = np.nanstd(num1_notxray)
num2_notxray_std = np.nanstd(num2_notxray)

### KS tests ###
from scipy import stats
[D_1_2, p_1_2] = stats.ks_2samp(num1, num2)
[D_tb1_x_nox, p_tb1_x_nox] = stats.ks_2samp(num1_xray, num1_notxray)
[D_tb2_x_nox, p_tb2_x_nox] = stats.ks_2samp(num2_xray, num2_notxray)

### plot hists ###

f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, sharex=True,figsize=[7,7])
#
#plt.figure()
a0.hist([num1, num2], bins=np.linspace(0,8,8),
         label=['2 arcsec', '1 arcsec'], histtype='step',
         density=True)
#plt.hist([num1, num2, all_n250], bins=np.linspace(0,6,50),
#         label=['2 arcsec', '1 arcsec', 'full cat'], histtype='step',
#         density=True)
a0.text(4.4,0.25,'p-value = '+str(round(p_1_2,3)))
a3.set_xlabel('Num of sources within '+str(arcsec_dist)+' arcsec')
a0.legend()

a2.vlines(num1_mean, 0, 3.5, 'C0', linestyle='dashdot', label=r'Mean Num')
a3.vlines(num2_mean, 0, 3.5, 'C1', label=r'Mean Num')
a2.axvspan(num1_mean-num1_std, num1_mean+num1_std, alpha=0.5, color='C0')
a3.axvspan(num2_mean-num2_std, num2_mean+num2_std, alpha=0.5, color='C1')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/local_enviro_comp_apersize_neg.png', overwrite=True)

#plt.figure()
f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, 
   sharex=True,figsize=[7,7])
a0.hist([num1_xray, num1_notxray], bins=np.linspace(0,8,8),
         color=['r','b'], label=['X-ray', 'Not X-ray'], histtype='step',
         density=True)
a3.set_xlabel('Num of sources within '+str(arcsec_dist)+' arcsec')
#plt.vlines(np.nanmedian(num1_xray), 0, 1, color='C0')
#plt.vlines(np.nanmedian(num1_notxray), 0, 1, color='C1')
a0.legend()
a0.set_title('2 arcsec Selection')
a0.text(4.4,0.25,'p-value = '+str(round(p_tb1_x_nox,3)))
#plt.xscale('log')

a2.vlines(num1_xray_mean, 0, 3.5, 'r', linestyle='dashdot', label=r'Mean Num')
a3.vlines(num1_notxray_mean, 0, 3.5, 'b', label=r'Mean Num')
a2.axvspan(num1_xray_mean-num1_xray_std, 
           num1_xray_mean+num1_xray_std, alpha=0.5, color='r')
a3.axvspan(num1_notxray_mean-num1_notxray_std, 
           num1_notxray_mean+num1_notxray_std, alpha=0.5, color='b')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/local_enviro_comp_xray_2_neg.png', overwrite=True)

#plt.figure()
f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, 
   sharex=True,figsize=[7,7])
a0.hist([num2_xray, num2_notxray], bins=np.linspace(0,8,8),
         color=['r','b'], label=['X-ray', 'Not X-ray'], histtype='step',
         density=True)
a3.set_xlabel('Num of sources within '+str(arcsec_dist)+' arcsec')
#plt.vlines(np.nanmedian(num2_xray), 0, 1, color='C0')
#plt.vlines(np.nanmedian(num2_notxray), 0, 1, color='C1')
a0.legend()
a0.set_title('1 arcsec Selection')
a0.text(4.4,0.25,'p-value = '+str(round(p_tb2_x_nox,3)))
#plt.xscale('log')

a2.vlines(num2_xray_mean, 0, 3.5, 'r', linestyle='dashdot', label=r'Mean Num')
a3.vlines(num2_notxray_mean, 0, 3.5, 'b', label=r'Mean Num')
a2.axvspan(num2_xray_mean-num2_xray_std, 
           num2_xray_mean+num2_xray_std, alpha=0.5, color='r')
a3.axvspan(num2_notxray_mean-num2_notxray_std, 
           num2_notxray_mean+num2_notxray_std, alpha=0.5, color='b')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/local_enviro_comp_xray_1_neg.png', overwrite=True)































