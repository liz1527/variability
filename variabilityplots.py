#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:10:23 2018

Investigate properties of the variable sources

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from matplotlib.mlab import griddata
plt.close('all') #close any open plots

# Import fits table
varydata = fits.open('variable_tables/variable_mag_flux_table_months_no99_DR11_modz.fits')[1].data

# Extract magnitude table and error table
mag = vari_funcs.mag5_months(varydata)
magerr = vari_funcs.magerr5_months(varydata)

## Change 99s to nans so they are ignored ###
mask = magerr >= 99
mag[mask] = np.nan
magerr[mask] = np.nan
varydata['KMAG_20'][varydata['KMAG_20']==99] = np.nan

## Remove rows where all are nans ###
mask = ~np.isnan(np.nanmean(mag, axis=0))
mag = mag[:,mask]
magerr = magerr[:,mask]
varydata = varydata[mask]

# Calculate excess variance
excess = vari_funcs.normsigmasq(mag, magerr)
mad = median_absolute_deviation(mag, axis=0, ignore_nan=True)

# remove anonamolous excess
excess[excess==np.nanmin(excess)] = np.nan

plt.figure()
plt.subplot(231)
plt.scatter(varydata['z_m2'], excess, c=varydata['X-ray'])
plt.xlabel('photometric redshift')
plt.ylabel('excess variance')
#plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)

plt.subplot(232)
plt.scatter(varydata['z_m2'], mad,c=varydata['X-ray'])
plt.xlabel('photometric redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)

plt.subplot(234)
plt.scatter(varydata['z_spec'], excess,c=varydata['X-ray'])
plt.xlabel('spectroscopic redshift')
plt.ylabel('excess variance')
#plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)

plt.subplot(235)
plt.scatter(varydata['z_spec'], mad,c=varydata['X-ray'])
plt.xlabel('photometric redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)

plt.subplot(233)
plt.scatter(varydata['z_m2'], varydata['mod_z_score'], c=varydata['X-ray'])
plt.xlabel('photometric redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)

plt.subplot(236)
plt.scatter(varydata['z_spec'], varydata['mod_z_score'],c=varydata['X-ray'])
plt.xlabel('spectroscopic redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)

plt.figure()
plt.subplot(231)
plt.scatter(varydata['z_m2'], excess, c=varydata['KMAG_20'])
plt.xlabel('photometric redshift')
plt.ylabel('excess variance')
#plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(232)
plt.scatter(varydata['z_m2'], mad,c=varydata['KMAG_20'])
plt.xlabel('photometric redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(234)
plt.scatter(varydata['z_spec'], excess,c=varydata['KMAG_20'])
plt.xlabel('spectroscopic redshift')
plt.ylabel('excess variance')
#plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(235)
plt.scatter(varydata['z_spec'], mad,c=varydata['KMAG_20'])
plt.xlabel('photometric redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(233)
plt.scatter(varydata['z_m2'], varydata['mod_z_score'], c=varydata['KMAG_20'])
plt.xlabel('photometric redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(236)
plt.scatter(varydata['z_spec'], varydata['mod_z_score'],c=varydata['KMAG_20'])
plt.xlabel('spectroscopic redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)
plt.colorbar()