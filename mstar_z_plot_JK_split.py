#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:23:28 2019

Code to plot M_star vs z

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
import matplotlib.colors as colors
plt.close('all')


### Get fits ###
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_K_extra_clean.fits')[1].data
dr11 = fits.open('UDS_catalogues/DR11-2arcsec-Jun-30-2019_best.fits')[1].data
varydata = fits.open('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')[1].data
fullxray = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019_best_chandra.fits')

### Split table into J and K selected ###
jdata, kdata, bothdata = vari_funcs.split_selection_band(varydata)

def ignore_zeros(m):
    m[m==0] = np.nan
    mask = ~np.isnan(m)
    m = m[mask]
    return m, mask

def prep_variables(tbdata):
    z = vari_funcs.get_z(tbdata) # get redshift
    m = tbdata['Mstar_z_p'] #get mass
    m, mask = ignore_zeros(m) #mask those with null mass
    z = z[mask] #apply mask to z array
    return z, m


z, m = prep_variables(dr11)
Jz, Jm = prep_variables(jdata)
Kz, Km = prep_variables(kdata)
bothz, bothm = prep_variables(bothdata)
allxz, allxm = prep_variables(fullxray)

#### KS test ###
#from scipy import stats
#[Dx_nox, px_nox] = stats.ks_2samp(noxm, xm)
#[Dallx_vary, pallx_vary] = stats.ks_2samp(varym, allxm)
##[Ddev_vary, pdev_vary] = stats.ks_2samp(varym, devm)
#[Dx_allx, px_allx] = stats.ks_2samp(xm, allxm)
#[Dnox_allx, pnox_allx] = stats.ks_2samp(noxm, allxm)

plt.figure()
plt.hist([Jm, Km, bothm], bins=np.logspace(4.5,12), color=['green','purple','y'], 
         histtype='step', label=['J Selected', 'K Selected', 'J and K Selected'])
plt.xlabel('$M_{star}$')
plt.ylabel('Number')
plt.xscale('log')
plt.tight_layout()

plt.figure(figsize=[10,7])
plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.2, label='Galaxy', rasterized=True)
plt.plot(allxz, allxm, 'ks',markersize=5, label='Non-Variable X-ray AGN')
plt.plot(Jz, Jm, 'o', color='green', label='Variable J selected AGN')
plt.plot(Kz, Km, 'o', color='purple', label='Variable K selected AGN')
plt.plot(bothz, bothm, 'o', color='y', label='Variable J and K selected AGN')

plt.yscale('log')
plt.xlim(xmin=-0.1, xmax=4.5)
plt.ylim(ymin=2e6, ymax=5e11)
plt.legend(loc='lower right')
plt.xlabel('Redshift')
plt.ylabel('Stellar Mass ($M_{\odot}$)')
plt.tight_layout()
