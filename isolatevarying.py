#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:04:36 2017

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

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best.fits')[1].data
chandata = fits.open('mag_flux_tables/month_xray_mag_flux_table_best.fits')[1].data
sdata = fits.open('mag_flux_tables/month_stars_mag_flux_table.fits')[1].data

### Restrict objects to those in the Chandra field ###
#tbdata = vari_funcs.chandra_only(tbdata)

## Create arrays of flux values but without 06B ###
flux = vari_funcs.mag5_months(tbdata)
fluxchan = vari_funcs.mag5_months(chandata) # for chandra non-stellar objects
sflux = vari_funcs.mag5_months(sdata)

### remove values that are +/-99 ###
#flux, tbdata = vari_funcs.no99(flux, tbdata)
#fluxchan, chandata = vari_funcs.no99(fluxchan, chandata)
#sflux, sdata = vari_funcs.no99(sflux, sdata)

### Change 99s to nans so they are ignored ###
flux[flux == 99] = np.nan
fluxchan[fluxchan == 99] = np.nan
sflux[sflux == 99] = np.nan

### Remove rows where all are nans ###
mask = ~np.isnan(np.nanmean(flux, axis=0))
flux = flux[:,mask]
tbdata = tbdata[mask]
mask = ~np.isnan(np.nanmean(fluxchan, axis=0))
fluxchan = fluxchan[:,mask]
chandata = chandata[mask]
mask = ~np.isnan(np.nanmean(sflux, axis=0))
sflux = sflux[:,mask]
sdata = sdata[mask]

print('Producing plot')
fig = vari_funcs.flux_variability_plot(flux, fluxchan, 'mad', starflux=sflux,
                                      stars=True)

#fig.canvas.mpl_connect('pick_event', vari_funcs.onpick)
## Find outliers ###
bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24, 25, 26])

##outliers, tb, modz = vari_funcs.find_outliers(fluxn, tbdata, bins)
##varys = tb[outliers]
##varyflux = vari_funcs.mag5_stacks(varys)
##varymad = median_absolute_deviation(varyflux, axis=1) 
##varymean = np.mean(varyflux, axis=1)
##plt.figure(1)
##plt.plot(varymean, varymad, 'md', mfc='none', markersize=10)
#

print('Identifying outliers')
outliers2, tb2, modz2 = vari_funcs.find_outliers(flux, tbdata, bins, threshold=5)
tb2['X-ray'][tb2['X-ray']==70] = False 
tb2['X-ray'][tb2['X-ray']==84] = True
### Find plotting values for the new flux table ###
flux2 = vari_funcs.mag5_months(tb2) 
flux2[flux2 == 99] = np.nan
avgfluxperob2 = np.nanmean(flux2, axis=0) 
vary2 = median_absolute_deviation(flux2, axis=0, ignore_nan=True)

print('Tabulating outliers')
### create table of varying ###
magmask = avgfluxperob2 < 21
outliers2 = outliers2*magmask
varydata = tb2[outliers2]
cols = fits.ColDefs(varydata)
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('variable_tables/variable_mag_flux_table_months_no99.fits')

varyfluxcorr = flux2[:,outliers2]
varymadcorr = vary2[outliers2]
varymeancorr = avgfluxperob2[outliers2]
varymodz = modz2[outliers2]
#plt.plot(varymeancorr, varymadcorr, 'kd', mfc='none', markersize=10)

#%% Find the ones that used to be variable but not aren't
#oldvary = fits.open('variable_mag_flux_table.fits')[1].data
#nums = ~np.isin(oldvary['NUMBER_05B'], varydata['NUMBER_05B'])
#oldvary = oldvary[nums]
#ind = np.isin(tbdata['NUMBER_05B'], oldvary['NUMBER_05B'])
#oldvary = tbdata[ind]

#varyflux = vari_funcs.mag5_stacks(oldvary)
#varymean = np.mean(varyflux, axis=1)
#varymad = median_absolute_deviation(varyflux, axis=1)
#plt.plot(varymean, varymad, 'gd', mfc='none', markersize=10)

#plt.figure()
#plt.scatter(varymeancorr, varymadcorr, c=varymodz, marker='+')
#plt.colorbar()

### Plot mod z in scatter ###
print('Plotting mod-z colours')
plt.figure()
plt.scatter(avgfluxperob2, vary2, c=modz2, marker='+', vmax=6)
cbar = plt.colorbar()
cbar.set_label('Modified z-score')
plt.yscale('log')
plt.ylim(1e-4, 1e1)
plt.xlabel('Mean Magnitude')
plt.ylabel('MAD Magnitude')
#
##plt.plot(varymeancorr, varymadcorr, 'md', mfc='none', markersize=10)

### plot mod z score as contours ###
print('Plotting contours')
xi = np.linspace(min(avgfluxperob2), max(avgfluxperob2), 1000)
yi = np.logspace(-4, 1, 1000)


zi = griddata(avgfluxperob2, vary2, modz2, xi, yi, interp='linear')
#plt.plot(avgfluxperob, vary, 'b+')
plt.figure(1)
plt.contour(xi, yi, zi, [0,4,5,6,8,10], zorder=3)
cbar = plt.colorbar()
cbar.set_label('Modified z-score')

#interest = tb2['NUMBER_05B'] == 1
#xint = avgfluxperob2[interest]
#yint = vary2[interest]
#plt.plot(xint, yint, 'ks', markersize=10)