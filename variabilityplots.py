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
varydata = fits.open('variable_tables/variable_month_3.5_alldetails.fits')[1].data

# Extract magnitude table and error table
mag = vari_funcs.mag5_months(varydata)
magerr = vari_funcs.magerr5_months(varydata)

## Change 99s to nans so they are ignored ###
mask = magerr >= 99
mag[mask] = np.nan
magerr[mask] = np.nan
varydata['KMAG_20'][varydata['KMAG_20']==99] = np.nan

### Remove rows where all are nans ###
mask = ~np.isnan(np.nanmean(mag, axis=0))
mag = mag[:,mask]
magerr = magerr[:,mask]
#varydata = varydata[mask]

# Calculate variance measures
excess = vari_funcs.normsigmasq(mag, magerr)
mad = varydata['MAD'][mask]
modz = varydata['mod_z_score'][mask]

# get redshifts
specz = varydata['z_spec'][mask]
photz = varydata['z_m2'][mask]

# get colours
c1 = varydata['X-ray'][mask]
c2 = varydata['KMAG_20'][mask]


plt.figure()
plt.subplot(231)
plt.scatter(photz, excess, c=c2)
plt.plot(photz[c1], excess[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('photometric redshift')
plt.ylabel('excess variance')
#plt.yscale('symlog')
#plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(232)
plt.scatter(photz, mad, c=c2)
plt.plot(photz[c1], mad[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('photometric redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(234)
plt.scatter(specz, excess, c=c2)
plt.plot(specz[c1], excess[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('spectroscopic redshift')
plt.ylabel('excess variance')
plt.ylim(ymin = np.nanmin(excess), ymax=np.nanmax(excess))
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(235)
plt.scatter(specz, mad, c=c2)
plt.plot(specz[c1], mad[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('spectroscopic redshift')
plt.ylabel('MAD')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(233)
plt.scatter(photz, modz, c=c2)
plt.plot(photz[c1], modz[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('photometric redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)
plt.colorbar()

plt.subplot(236)
plt.scatter(specz, modz, c=c2)
plt.plot(specz[c1], modz[c1],'ro', mfc = 'none', markersize = 10)
plt.xlabel('spectroscopic redshift')
plt.ylabel('Modified z_score')
plt.xlim(xmin=0)
plt.colorbar()