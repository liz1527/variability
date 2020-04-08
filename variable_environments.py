#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:49:08 2020

Code to compare environments of variable samples

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table, Column, vstack #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots

### Open required tables ###
tbdata1 = Table.read('variable_tables/K/variables_no06_chi30_DR11data.fits')
tbdata2 = Table.read('variable_tables/K/variables_no06_chi30_1arcsec_noneg_DR11data.fits')

### Open environment tables ###
envirodata = Table.read('UDS_catalogues/DR11-aperture-densities-DTM-Jun30.fits')

### Get enviroment info for each table ###
def get_data(tbdata, envirodata):
    tbID = tbdata['ID']
    enID = envirodata['DR11-ID']
    mask = np.isin(enID, tbID)
    tb_enviro = envirodata[mask]
    enID2 = tb_enviro['DR11-ID']
    mask = np.isin(tbID, enID2)
    tbdatanew = tbdata[mask]
    tb_n250 = tb_enviro['n-250']
    return tb_enviro, tbdatanew, tb_n250

tb1_enviro, tbdata1, tb1_n250 = get_data(tbdata1, envirodata)
tb2_enviro, tbdata2, tb2_n250 = get_data(tbdata2, envirodata)
all_n250 = envirodata['n-250']

tb1_n250_mean = np.nanmean(tb1_n250)
tb2_n250_mean = np.nanmean(tb2_n250)
tb1_n250_std = np.nanstd(tb1_n250)
tb2_n250_std = np.nanstd(tb2_n250)

### Split by X-ray ###
tb1_n250_xray = tb1_n250[tbdata1['X-ray']==True]
tb1_n250_notxray = tb1_n250[tbdata1['X-ray']==False]
tb2_n250_xray = tb2_n250[tbdata2['X-ray']==True]
tb2_n250_notxray = tb2_n250[tbdata2['X-ray']==False]

tb1_n250_xray_mean = np.nanmean(tb1_n250_xray)
tb2_n250_xray_mean = np.nanmean(tb2_n250_xray)
tb1_n250_xray_std = np.nanstd(tb1_n250_xray)
tb2_n250_xray_std = np.nanstd(tb2_n250_xray)
tb1_n250_notxray_mean = np.nanmean(tb1_n250_notxray)
tb2_n250_notxray_mean = np.nanmean(tb2_n250_notxray)
tb1_n250_notxray_std = np.nanstd(tb1_n250_notxray)
tb2_n250_notxray_std = np.nanstd(tb2_n250_notxray)

### KS tests ###
from scipy import stats
[D_1_2, p_1_2] = stats.ks_2samp(tb1_n250, tb2_n250)
[D_tb1_x_nox, p_tb1_x_nox] = stats.ks_2samp(tb1_n250_xray, tb1_n250_notxray)
[D_tb2_x_nox, p_tb2_x_nox] = stats.ks_2samp(tb2_n250_xray, tb2_n250_notxray)

### plot hists ###

f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, sharex=True,figsize=[7,7])
#
#plt.figure()
a0.hist([tb1_n250, tb2_n250], bins=np.linspace(0,5,30),
         label=['2 arcsec', '1 arcsec'], histtype='step',
         density=True)
#plt.hist([tb1_n250, tb2_n250, all_n250], bins=np.linspace(0,6,50),
#         label=['2 arcsec', '1 arcsec', 'full cat'], histtype='step',
#         density=True)
a0.text(3,0.4,'p-value = '+str(round(p_1_2,3)))
a3.set_xlabel('Num density in 250 kpc aperture')
a0.legend()

a2.vlines(tb1_n250_mean, 0, 3.5, 'C0', linestyle='dashdot', label=r'Mean $n_{250}$')
a3.vlines(tb2_n250_mean, 0, 3.5, 'C1', label=r'Mean $n_{250}$')
a2.axvspan(tb1_n250_mean-tb1_n250_std, tb1_n250_mean+tb1_n250_std, alpha=0.5, color='C0')
a3.axvspan(tb2_n250_mean-tb2_n250_std, tb2_n250_mean+tb2_n250_std, alpha=0.5, color='C1')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/enviro_comp_apersize_noneg.png', overwrite=True)

#plt.figure()
f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, 
   sharex=True,figsize=[7,7])
a0.hist([tb1_n250_xray, tb1_n250_notxray], bins=np.linspace(0,5,30),
         color=['r','b'], label=['X-ray', 'Not X-ray'], histtype='step',
         density=True)
a3.set_xlabel('Num density in 250 kpc aperture')
#plt.vlines(np.nanmedian(tb1_n250_xray), 0, 1, color='C0')
#plt.vlines(np.nanmedian(tb1_n250_notxray), 0, 1, color='C1')
a0.legend()
a0.set_title('2 arcsec Selection')
a0.text(3,0.4,'p-value = '+str(round(p_tb1_x_nox,3)))
#plt.xscale('log')

a2.vlines(tb1_n250_xray_mean, 0, 3.5, 'r', linestyle='dashdot', label=r'Mean $n_{250}$')
a3.vlines(tb1_n250_notxray_mean, 0, 3.5, 'b', label=r'Mean $n_{250}$')
a2.axvspan(tb1_n250_xray_mean-tb1_n250_xray_std, 
           tb1_n250_xray_mean+tb1_n250_xray_std, alpha=0.5, color='r')
a3.axvspan(tb1_n250_notxray_mean-tb1_n250_notxray_std, 
           tb1_n250_notxray_mean+tb1_n250_notxray_std, alpha=0.5, color='b')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/enviro_comp_xray_2_noneg.png', overwrite=True)

#plt.figure()
f, (a0, a2, a3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [10, 1, 1]}, 
   sharex=True,figsize=[7,7])
a0.hist([tb2_n250_xray, tb2_n250_notxray], bins=np.linspace(0,5,30),
         color=['r','b'], label=['X-ray', 'Not X-ray'], histtype='step',
         density=True)
a3.set_xlabel('Num density in 250 kpc aperture')
#plt.vlines(np.nanmedian(tb2_n250_xray), 0, 1, color='C0')
#plt.vlines(np.nanmedian(tb2_n250_notxray), 0, 1, color='C1')
a0.legend()
a0.set_title('1 arcsec Selection')
a0.text(3,0.4,'p-value = '+str(round(p_tb2_x_nox,3)))
#plt.xscale('log')

a2.vlines(tb2_n250_xray_mean, 0, 3.5, 'r', linestyle='dashdot', label=r'Mean $n_{250}$')
a3.vlines(tb2_n250_notxray_mean, 0, 3.5, 'b', label=r'Mean $n_{250}$')
a2.axvspan(tb2_n250_xray_mean-tb2_n250_xray_std, 
           tb2_n250_xray_mean+tb2_n250_xray_std, alpha=0.5, color='r')
a3.axvspan(tb2_n250_notxray_mean-tb2_n250_notxray_std, 
           tb2_n250_notxray_mean+tb2_n250_notxray_std, alpha=0.5, color='b')

a2.set_yticks([])
a3.set_yticks([])

a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
f.savefig('plots/new_catalogue/environment/enviro_comp_xray_1_noneg.png', overwrite=True)

